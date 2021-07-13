#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021, Yalin Li <yalinli2@illinois.edu>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

'''
Unit construction and functions for creating wastewwater treatment system.

References
----------
[1] Humbird et al., Process Design and Economics for Biochemical Conversion of
Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic
Hydrolysis of Corn Stover; Technical Report NREL/TP-5100-47764;
National Renewable Energy Lab (NREL), 2011.
https://www.nrel.gov/docs/fy11osti/47764.pdf

[2] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic
Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update;
NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018.
https://doi.org/10.2172/1483234

TODO:
    AerobicDigestion and MembraneBioreactor will be replaced by
    an anaerobic membrane bioresactor (AnMBR).
'''


# %%

import thermosteam as tmo
import biosteam as bst
from biosteam.units.decorators import cost
from thermosteam import Stream, separations

#!!! Need to enable relative importing
from _chemicals import default_insolubles, get_insoluble_IDs, get_soluble_IDs
from utils import (
    auom,
    get_BD_dct,
    compute_stream_COD,
    get_digestion_rxns,
    get_MB_split,
    )
from _settings import new_price
from _internal_circulation_rx import InternalCirculationRx

_MGD_2_m3hr = auom('gallon').conversion_factor('m3')*1e6/24
_GPM_2_m3hr = auom('gallon').conversion_factor('m3')*60
_Gcal_2_kJ = auom('kcal').conversion_factor('kJ')*1e6 # (also MMkcal/hr)

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction
CEPCI = bst.units.design_tools.CEPCI_by_year
Unit = bst.Unit

__all__ = ('create_wastewater_treatment_system',)


# %%

@cost(basis='COD flow', ID='Ammonia addition', units='kg-O2/hr',
      kW=13.4226, cost=195200, S=5600, CE=CEPCI[2012], n=0.6, BM=1.5)
@cost(basis='COD flow', ID='Caustic feed', units='kg-O2/hr',
      kW=4.4742, cost=20000, S=5600, CE=CEPCI[2012], n=0.6, BM=3)
@cost(basis='Volumetric flow', ID='Aerobic basin', units='m3/hr',
      # power usage including polymer addition and feed pumping
      kW=1.4914+134.226,
      # 2.7 in million gallons per day (MGD)
      cost=4804854, S=2.7*_MGD_2_m3hr, CE=CEPCI[2012], n=0.6, BM=2.1)
@cost(basis='COD flow', ID='Blowers', units='kg-O2/hr',
      kW=6711.3, cost=2070000, S=5600, CE=CEPCI[2012], n=0.6, BM=2)
class AerobicDigestion(Unit):
    _N_ins = 6
    _N_outs = 2
    _units= {'COD flow': 'kg-O2/hr',
             'Volumetric flow': 'm3/hr'}

    # 4350 and 356069 are water flows from streams 622 and 611 in ref [1]
    evaporation = 4350 / 356069

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 caustic_mass=0, need_ammonia=False):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.caustic_mass = caustic_mass
        self.need_ammonia = need_ammonia

        # Reactions from auto-populated combustion reactions.
        # Based on P49 in ref [1], 96% of remaining soluble organic matter
        # is removed after aerobic digestion, of which 74% is converted to
        # water and CO2 and 22% to cell mass
        chems = self.chemicals
        growth_rxns = get_digestion_rxns(self.ins[0],
                                         BD=get_BD_dct(chems),
                                         X_biogas=0., X_growth=0.22,
                                         biomass_ID='WWTsludge')

        combustion_rxns = chems.get_combustion_reactions()
        self.digestion_rxns = ParallelRxn([*[i*0.74 for i in combustion_rxns
                                             if i.reactant in growth_rxns.reactants],
                                           *growth_rxns])

        #                                 Reaction definition       Reactant Conversion
        self.nitrification_rxn = Rxn('NH4OH + 2 O2 -> HNO3 + 2 H2O', 'NH4OH',  1)

        self.neutralization_rxns = ParallelRxn([
            #              Reaction definition       Reactant Conversion
            Rxn('H2SO4 + 2 NaOH -> Na2SO4 + 2 H2O',  'H2SO4',    1),
            Rxn('HNO3 + NaOH -> NaNO3 + H2O',        'HNO3',     1)
            ])

    def _run(self):
        influent, recycle, caustic, ammonia, polymer, air = self.ins
        vent, effluent = self.outs
        vent.phase = 'g'

        caustic.imass['NaOH'] = self.caustic_mass
        # Ammonia as a nutrient
        if self.need_ammonia is True:
            # Based on Table 33 on Page 73 of ref [2], originally as NH3
            ammonia.imass['NH4OH'] = 36 * 35.046/17.031
        else: ammonia.empty()

        effluent.mix_from(self.ins[0:5])

        self.design_results['Volumetric flow'] = influent.F_vol
        self.nitrification_rxn.force_reaction(effluent.mol)
        self.neutralization_rxns.force_reaction(effluent.mol)
        if effluent.imass['NaOH'] < 0:
            caustic.imass['NaOH'] += -effluent.imol['NaOH']
            effluent.imol['NaOH'] = 0

        # Ratio based on equipment sizing in ref [2]
        ratio = self.design_results['Volumetric flow'] / (2*_MGD_2_m3hr)
        polymer.imass['Polymer'] = 2 * ratio

        # 4693, 54718 and 180206 from stream 601 in ref [2]
        air.imass['O2'] = 54718 * ratio
        air.imass['N2'] = 180206 * ratio
        effluent.mix_from([effluent, air])

        self.design_results['COD flow'] = compute_stream_COD(effluent)
        self.digestion_rxns(effluent.mol)
        vent.copy_flow(effluent, ('CO2', 'O2', 'N2'), remove=True)
        vent.imol['Water'] = influent.imol['Water'] * self.evaporation
        effluent.imol['Water'] -= vent.imol['Water']
        vent.T = effluent.T

    def _design(self):
        self._decorated_cost()
        if self.need_ammonia == False:
            self.cost_items['Ammonia addition'].cost = 0


@cost(basis='Volumetric flow', ID='Reactor', units='m3/hr',
      # power usage including pumps
      kW=63.3845+715.872+14.914,
      # 2.7 in million gallons per day (MGD)
      cost=4898500, S=2.7*_MGD_2_m3hr, CE=CEPCI[2012], n=1, BM=1.6)
@cost(basis='COD flow', ID='Conveyor', units='kg-O2/hr',
      kW=7.457, cost=7000, S=5600, CE=CEPCI[2012], n=0.6, BM=2.9)
class MembraneBioreactor(bst.Splitter):
    _N_ins = 1
    _N_outs = 2
    _units= {'Volumetric flow': 'm3/hr',
             'COD flow': 'kg-O2/hr'}

    def __init__(self, ID='', ins=None, outs=(), thermo=None, split=None):
        self._load_thermo(thermo)
        if split is None:
            split = get_MB_split(self.chemicals, split)
        bst.Splitter.__init__(self, ID, ins, outs, thermo, split=split, order=None)

    def _run(self):
        mixture = self.ins[0].copy()
        mixture.mix_from(self.ins)
        separations.split(mixture, *self.outs, self.split)

    def _design(self):
        self.design_results['Volumetric flow'] = self.outs[0].F_vol
        self.design_results['COD flow'] = compute_stream_COD(self.ins[0])


@cost(basis='COD flow', ID='Thickeners', units='kg-O2/hr',
      kW=107.3808, cost=750000, S=5600, CE=CEPCI[2012], n=0.6, BM=1.6)
class BeltThickener(Unit):
    _ins_size_is_fixed = False
    _N_outs = 2
    _units= {'COD flow': 'kg-O2/hr'}

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 insolubles=default_insolubles):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.insolubles = get_insoluble_IDs(self.chemicals, insolubles)

    def _run(self):
        centrate, solids = self.outs
        insolubles = self.insolubles
        solubles = get_soluble_IDs(self.chemicals, insolubles)

        influent = self.ins[0].copy()
        influent.mix_from(self.ins)
        solids.copy_flow(influent, insolubles)
        # Concentrate sludge to 4% solids
        solids.imass['Water'] = 0.96/0.04 * influent.imass[insolubles].sum()
        if solids.imass['Water'] < influent.imass['Water']:
            ratio = solids.imass['Water'] / influent.imass['Water']
            solids.imass[solubles] = ratio * influent.imass[solubles]
            solids.T = influent.T

            centrate.mol = influent.mol - solids.mol
            centrate.T = influent.T
        else:
            centrate.empty()
            solids.copy_like(influent)

    def _design(self):
        self.design_results['COD flow'] = compute_stream_COD(self.ins[0])


@cost(basis='COD flow', ID='Centrifuge', units='kg-O2/hr',
      # power usage includes feed pumping and centrifuge
      kW=22.371+123.0405, cost=686800, S=5600, CE=CEPCI[2012], n=0.6, BM=2.7)
class SludgeCentrifuge(Unit):
    _N_ins = 1
    _N_outs = 2
    _units= {'COD flow': 'kg-O2/hr'}

    __init__ = BeltThickener.__init__

    def _run(self):
        influent = self.ins[0]
        centrate, solids = self.outs
        centrate.T = solids.T = influent.T
        insolubles = self.insolubles
        solubles = get_soluble_IDs(self.chemicals, insolubles)

        # Centrifuge captures 95% of the solids at 20% solids
        solids.imass[insolubles] = 0.95 * influent.imass[insolubles]
        solids.imass['Water'] = 0.8/0.2 * (influent.imass[insolubles].sum())
        if solids.imass['Water'] < influent.imass['Water']:
            ratio = solids.imass['Water'] / influent.imass['Water']
            solids.imass[solubles] = ratio * influent.imass[solubles]

            centrate.mol = influent.mol - solids.mol
        else:
            centrate.empty()
            solids.copy_like(influent)

    _design = BeltThickener._design


@cost(basis='Volumetric flow', ID='Reactor', units='m3/hr',
      # 2.7 in million gallons per day (MGD)
      cost=2450000, S=2.7*_MGD_2_m3hr, CE=CEPCI[2012], n=1, BM=1.8)
@cost(basis='Volumetric flow', ID='Evaporator', units='m3/hr',
      # 2.7 in million gallons per day (MGD)
      kW=1103.636, cost=5000000, S=2.7*_MGD_2_m3hr, CE=CEPCI[2012], n=0.6, BM=1.6)
class ReverseOsmosis(Unit):
    _N_ins = 1
    _N_outs = 2
    _units = {'Volumetric flow': 'm3/hr'}

    def _run(self):
        influent = self.ins[0]
        water, brine = self.outs

        self.design_results['Volumetric flow'] = self.F_vol_in

        # Based on stream 626 and 627 in ref [1]
        water.imass['Water'] = 376324/(376324+4967) * influent.imass['Water']
        brine.mol = influent.mol - water.mol
        water.T = brine.T = influent.T


def create_wastewater_treatment_units(ins, outs, IC_method,
                                      dry_flow_tpd, need_ammonia):
    wwt_streams = ins
    biogas, vent_R602, S604_CHP, recycled_water, brine = outs
    vent_R602.phase = 'g'

    ammonia_R602 = Stream('ammonia_R602', units='kg/hr')
    caustic_R602 = Stream('caustic_R602', units='kg/hr', price=new_price['Caustics'])
    polymer_R602 = Stream('polymer_R602', units='kg/hr', price=new_price['Polymer'])
    air_R602 = Stream('air_R602', phase='g', units='kg/hr')

    ######################## Units ########################
    # Mix waste liquids for treatment
    M601 = bst.units.Mixer('M601', ins=wwt_streams)

    R601 = InternalCirculationRx('R601', ins=M601-0,
                                 outs=(biogas, 'IC_eff', 'IC_sludge'),
                                 method=IC_method)

    R602 = AerobicDigestion('R602',
                            ins=(R601-1, '', caustic_R602,
                                 ammonia_R602, polymer_R602, air_R602),
                                 outs=(vent_R602, 'aerobic_treated_water'),
                                 caustic_mass=2252*dry_flow_tpd/2205,
                                 need_ammonia=need_ammonia)

    # Membrane bioreactor to split treated wastewater from R602
    S601 = MembraneBioreactor('S601', ins=R602-1,
                              outs=('membrane_treated_water', 'membrane_sludge'))

    # Recycled sludge stream of memberane bioreactor, the majority of it (96%)
    # goes to aerobic digestion
    S602 = bst.units.Splitter('S602', ins=S601-1, outs=('to_aerobic_digestion', ''),
                              split=0.96)

    # S603 = BeltThickener('S603', ins=(R601-2, S602-1),
    S603 = BeltThickener('S603', ins=(R601-2, S602-1),
                         outs=('S603_centrate', 'S603_solids'))

    # Sludge centrifuge to separate water (centrate) from sludge
    # Ref [1] included polymer addition in process flow diagram, but did not include
    # in the variable operating cost, thus followed ref [4] to add polymer in AerobicDigestion
    S604 = SludgeCentrifuge('S604', ins=S603-1,
                            outs=('S604_centrate', S604_CHP))

    # Mix recycles to aerobic digestion
    bst.units.Mixer('M602', ins=(S602-0, S603-0, S604-0), outs=1-R602)

    # Reverse osmosis to treat membrane separated water
    ReverseOsmosis('S605', ins=S601-0, outs=(recycled_water, brine))


create_wastewater_treatment_system = bst.SystemFactory(
    f=create_wastewater_treatment_units,
    ID='wastewater_treatment_sys',
    outs=[dict(ID='biogas'),
          dict(ID='vent_R602'),
          dict(ID='S604_CHP'),
          dict(ID='recycled_water'),
          dict(ID='brine')],
    fixed_ins_size=False,
)