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
'''


# %%

import thermosteam as tmo
import biosteam as bst
from biosteam.units.decorators import cost

# from biorefineries.wwt import (
#     default_insolubles, get_insoluble_IDs, get_soluble_IDs,
#     InternalCirculationRx,
#     FilterTank,
#     AnMBR
#     )
# from biorefineries.wwt.utils import auom, compute_stream_COD
from _chemicals import default_insolubles, get_insoluble_IDs, get_soluble_IDs
from utils import auom, compute_stream_COD
from _internal_circulation_rx import InternalCirculationRx
from _filter_tank import FilterTank
from _membrane_bioreactor import AnMBR

_mgd_to_cmh = auom('gallon').conversion_factor('m3')*1e6/24
_gpm_to_cmh = auom('gallon').conversion_factor('m3')*60
_Gcal_to_kJ = auom('kcal').conversion_factor('kJ')*1e6 # (also MMkcal/hr)

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction
CEPCI = bst.units.design_tools.CEPCI_by_year
Unit = bst.Unit

__all__ = ('create_wastewater_treatment_system',)


# %%

# =============================================================================
# Other units
# =============================================================================
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
      cost=2450000, S=2.7*_mgd_to_cmh, CE=CEPCI[2012], n=1, BM=1.8)
@cost(basis='Volumetric flow', ID='Evaporator', units='m3/hr',
      # 2.7 in million gallons per day (MGD)
      kW=1103.636, cost=5000000, S=2.7*_mgd_to_cmh, CE=CEPCI[2012], n=0.6, BM=1.6)
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


# %%

# =============================================================================
# System function
# =============================================================================
def create_wastewater_treatment_units(ins, outs, IC_method):
    wwt_streams = ins
    biogas, S604_CHP, recycled_water, brine = outs

    ######################## Units ########################
    # Mix waste liquids for treatment
    M601 = bst.units.Mixer('M601', ins=wwt_streams)

    R601 = InternalCirculationRx('R601', ins=M601-0,
                                 outs=('biogas_R601', 'IC_eff', 'IC_sludge'),
                                 method=IC_method)

    R602 = AnMBR('R602', ins=(R601-1, '', 'naocl_R602', 'citric_R602',
                              'bisulfite', 'air_R602'),
                 outs=('biogas_R602', 'permeate_R602', 'sludge_R602', 'vent_R602'),
                 reactor_type='CSTR',
                 membrane_configuration='cross-flow',
                 membrane_type='multi-tube',
                 membrane_material='ceramic',
                 include_aerobic_filter=False,
                 add_GAC=False,
                 include_degassing_membrane=True)

    R603 = FilterTank('R603', ins=(R602-1, '', 'air_R603'),
                      outs=('biogas_R603', 'treated_R603', 'sludge_R603', 'vent_R603'),
                      filter_type='aerobic')

    bst.units.Mixer('M602', ins=(R601-0, R602-0, R603-0), outs=biogas)

    # Recycled sludge stream of memberane bioreactor, the majority of it (96%)
    # goes back to the aerobic filter
    S602 = bst.units.Splitter('S602', ins=R603-2, outs=('recycled_S602', 'wasted_S602'),
                              split=0.96)

    S603 = BeltThickener('S603', ins=(R601-2, R602-2, S602-1),
                         outs=('centrate_S603', 'solids_S603'))

    # Sludge centrifuge to separate water (centrate) from sludge
    # Ref [1] included polymer addition in process flow diagram, but did not include
    # in the variable operating cost, thus followed ref [4] to add polymer in AerobicDigestion
    S604 = SludgeCentrifuge('S604', ins=S603-1,
                            outs=('centrate_S604', S604_CHP))

    # Mix recycles to aerobic digestion
    bst.units.Mixer('M603', ins=(S602-0, S603-0, S604-0), outs=1-R603)

    # Reverse osmosis to treat aerobically polished water
    ReverseOsmosis('S605', ins=R603-1, outs=(recycled_water, brine))


create_wastewater_treatment_system = bst.SystemFactory(
    f=create_wastewater_treatment_units,
    ID='wastewater_treatment_sys',
    outs=[dict(ID='biogas'),
          dict(ID='S604_CHP'),
          dict(ID='recycled_water'),
          dict(ID='brine')],
    fixed_ins_size=False,
)