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
TODO:
	Need to consider the endogenous decay coefficient b
	Need to consider the C/N ratio, George's analysis shows it's very low
	Need to consider sulfate and sulfide
	Check with Brian's AnMBR paper and see the COD<1300 mg/L not preferable thing
'''

import numpy as np
import biosteam as bst
import thermosteam as tmo
from collections import defaultdict
from chemicals.elements import molecular_weight
from biorefineries import cornstover as cs

Rxn = tmo.reaction.Reaction
PRxn = tmo.reaction.ParallelReaction

s1 = cs.R601.ins[0] # influent into AD
# from .design_tools.tank_design import (
#     TankPurchaseCostAlgorithm,
#     compute_number_of_tanks_and_total_purchase_cost,
#     storage_tank_purchase_cost_algorithms,
#     mix_tank_purchase_cost_algorithms)

__all__ = ('IC',)


# %%

# =============================================================================
# Util functions
# =============================================================================

def get_CHONSP(chemical):
    organic = True
    atoms = chemical.atoms

    CHONSP = []
    for atom in ('C', 'H', 'O', 'N', 'S', 'P'):
        CHONSP.append(atoms.get(atom) or 0.,)

    if CHONSP[0] <= 0 or CHONSP[1] <= 0: # does not have C or H
        if not (len(atoms) == 1 and CHONSP[1] == 2): # count H2 as organic
            organic = False

    if sum(v for v in atoms.values()) != sum(CHONSP): # contains other elements
        organic = False

    return CHONSP if organic else [0.]*6


def get_COD_stoichiometry(chemical):
    r'''
    Get the molar stoichiometry for the theoretical
    chemical oxygen demand (COD) of a given chemical as in:

    .. math::
        C_nH_aO_bN_cS_dP_e + \frac{2n+0.5a-b-1.5c+3d+2.5e}{2}O_2
        -> nCO_2 + \frac{a-3c}{2}H_2O + cNH_3 + dSO_3 + \frac{e}{2}P_2O_5
    '''
    Xs = nC, nH, nO, nN, nS, nP = get_CHONSP(chemical)

    dct = {
        chemical.formula: -1. if sum(Xs)!=0 else 0.,
        'O2': -(nC+0.25*nH-0.5*nO-0.75*nN+1.5*nS+1.25*nP),
        'CO2': nC,
        'H2O': 0.5*nH-1.5*nN,
        'NH3': nN,
        'SO3': nS,
        'P2O5': 0.5*nP
        }

    return dct


def get_BMP_stoichiometry(chemical):
    r'''
    Compute the theoretical biochemical methane potential (BMP) in
    mol :math:`CH_4`/mol chemical of a given chemical using:

    .. math::
        C_vH_wO_xN_yS_z + \frac{4v-w-2x+3y+2z}{2}H2O ->
        \frac{4v+w-2x-3y-2z}{8}CH4 + \frac{(4v-w+2x+3y+2z}{8}CO2 + yNH_3 + zH_2S
    '''
    Xs = nC, nH, nO, nN, nS, nP = get_CHONSP(chemical)

    dct = {
        chemical.formula: -1. if sum(Xs)!=0 else 0.,
        'H2O': -(nC-0.25*nH-0.5*nO+0.75*nN+0.5*nS),
        'CH4': 0.5*nC+0.125*nH-0.25*nO-0.375*nN-0.25*nS,
        'CO2': 0.5*nC-0.125*nH+0.25*nO+0.375*nN+0.25*nS,
        'NH3': nN,
        'H2S': nS
        }

    return dct


# Biodegradability, 0.87 from glucose (treated as the maximum value)
#!!! Add other sugars, etc
def get_BD_dct(default=0.87, BD=None):
    BD_dct = defaultdict(lambda: default)

    # Based on Kontos thesis
    BD_dct['AceticAcid'] = 0.87
    BD_dct['Arabinose'] = 0.2
    BD_dct['Glucose'] = 0.87
    BD_dct['GlucoseOligomer'] = BD_dct['Glucan'] = 0.81
    BD_dct['HMF'] = 0.85
    BD_dct['LacticAcid'] = 0.85
    BD_dct['Lignin'] = BD_dct['SolubleLignin'] = 0.001
    BD_dct['Tar'] = 0.
    BD_dct['Xylose'] = 0.8
    BD_dct['XyloseOligomer'] = BD_dct['Xylan'] = 0.75

    # Other assumptions
    C6_oligomer_to_monomer = BD_dct['GlucoseOligomer'] / BD_dct['Glucose']
    C5_oligomer_to_monomer = BD_dct['XyloseOligomer'] / BD_dct['Xylose']

    BD_dct['GalactoseOligomer'] = BD_dct['Galactan'] = \
        C6_oligomer_to_monomer * BD_dct['Galactose']
    BD_dct['MannoseOligomer'] = BD_dct['Mannan'] = \
        C6_oligomer_to_monomer * BD_dct['Mannose']
    BD_dct['ArabinoseOligomer'] = BD_dct['Arabinan'] = \
        C5_oligomer_to_monomer * BD_dct['Arabinose']

    if BD:
        BD_dct.update(BD)

    return BD_dct


def compute_stream_COD(stream, BD_dct=None):
    r'''
    Compute the chemical oxygen demand (COD) of a given stream in kg/m3
    by summing the COD of each chemical in the stream using:

    .. math::
        COD [\frac{kg}{m^3}] = mol_{chemical} [\frac{kmol}{m^3}] * \frac{g O_2}{mol chemical}
    '''
    chems = stream.chemicals
    imol = stream.z_mol
    iCOD = np.array([-get_COD_stoichiometry(i)['O2'] for i in chems])

    # Consider biodegradability
    if BD_dct:
        iCOD *= np.array([BD_dct[i.ID] for i in chems])

    COD = (imol*iCOD).sum() / stream.F_vol * molecular_weight({'O': 2})
    return COD


# The final conversion should be the conversion of the biodegradable part
def get_AD_rxns(stream, BD_dct=defaultdict(lambda: 0.87), conversion=1.):
    chems = stream.chemicals
    iBD = np.array([BD_dct[i.ID] for i in chems])
    iX = conversion * iBD

    # for i in chems:



    return iBD


def get_growth_rxns(stream, BD_dct=defaultdict(lambda: 0.87), conversion=1.):
    '''NOT READY YET.'''



# %%

class IC(bst.MixTank):
    '''
    Internal circulation (IC) reactor for anaerobic digestion,
    including a high-rate bottom reactor for rapid organic removal and
    a low-rate top reactor for polishing.
    Both reactors are similar to upflow anaerobic blanket reactor (UASB).

    Design of the reactor follows steps described in [1]_
    (assuming steady state and pseudo-zeroth-order kinetics),
    where two methods are used based on Irizar et al.[2]_ and
    Tchobanoglous et al.[3]_.

    Parameters
    ----------
    method : str
        Either "separate" to design the bottom and top reactors separately as in [2]_,
        or "lumped" to design the entire IC reactor as a black box following [3]_.
    OLR_overall : float
        Overall organic loading rate, [kg COD/m3/hr].
    COD_rm : float
        COD removal.
    qw : float
        Ratio between the bottom reactor waste flow and the influent.
    Y : float
        Biomass yield, [kg biomass/kg substrate COD].
    mu_max : float
        Maximum specific growth rate, [/hr].
    b : float
        Specific endogenous decay coefficient, [/hr].

    tau : float
        Residence time [hr].
    V_wf : float
        Fraction of working volume over total volume.
    vessel_material : str
        Vessel material.
    kW_per_m3 : float
        Electricity requirement per unit volume, [kW/m^3].
        Default to 0 as IC reactors realizes mixing through internal circulation
        caused by the rising force of the generated biogas.
    kwargs : dict
        Other keyword arguments (e.g., F_Xb, F_Xt).

    References
    ----------
    .. [1] Kontos, G. A. Advanced Anaerobic Treatment for Energy Recovery and
    Improved Process Economics in the Management of Biorefinery Wastewaters,
    University of Illinois at Urbana-Champaign, Champaign, IL, 2021.

    .. [2] Irizar et al., Model-Based Design of a Software Sensor for Real-Time
    Diagnosis of the Stability Conditions in High-Rate Anaerobic Reactors –
    Full-Scale Application to Internal Circulation Technology.
    Water Research 2018, 143, 479–491.
    `<https://doi.org/10.1016/j.watres.2018.06.055>`_.

    .. [3] Tchobanoglous et al., Wastewater Engineering: Treatment and Resource Recovery,
    5th ed.; McGraw-Hill Education: New York, 2013.


    '''

    _N_ins = 1
    _N_outs = 3
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True

    # Design and operating assumptions
    _default_kW_per_m3 = 0.
    _F_Xb = 0.0032
    _F_Xt = 0.0281

    # Utils
    _get_COD_stoichiometry = get_COD_stoichiometry
    _get_BMP_stoichiometry = get_BMP_stoichiometry

    _get_BD_dct = get_BD_dct
    _BD_dct = get_BD_dct()

    _compute_stream_COD = compute_stream_COD

    #!!! Need to figure this out
    # endogenous_rxn = Rxn('WWTsludge -> CO2 + XX ', 'WWTsludge', 1.)


    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 method='separate', OLR_overall=1.25, COD_rm=0.8, qw=0.05,
                 Y=0.07, mu_max=0.01, b=0.00083,
                 tau=None, V_wf=0.8, vessel_material='Stainless steel',
                 kW_per_m3=0., **kwargs):
        bst.MixTank.__init__(self, ID, ins, outs, thermo,
                             tau=tau, V_wf=V_wf,
                             vessel_material=vessel_material,
                             kW_per_m3=kW_per_m3)
        self.method = method
        self.OLR_overall = OLR_overall
        self.COD_rm = COD_rm
        self.qw = qw
        self.Y = Y
        self.mu_max = mu_max
        self.b = b

        for k, v in kwargs.items():
            setattr(self, k, v)


    def _run(self):
        method = self.method.lower()
        if method == 'separate':
            self._run_separate()
        elif method == 'lumped':
            self._run_lumped()


    def _run_separate(self):
        inf, = self.ins
        eff, waste, biogas = self.outs
        cal_COD = self._compute_stream_COD
        BD_dct = self.BD_dct
        OLR_overall = self.OLR_overall # overall organic loading rate, kg COD/m3/hr
        Xinf = inf.imass['WWTsludge'] # influent biomass concentration, kg/hr

        Qi= inf.F_vol # influent volumetric flow rate, m3/hr
        Si = cal_COD(inf, BD_dct) # influent COD, kg/m3
        Vtot = Qi * Si / OLR_overall # total volume (bottom and top rxs), m3
        Qw = self.qw * Qi # waste volumetric flow rate, m3/hr
        Qe = Qb = Qt = Qi - Qw # effluent volumetric flow rate, m3/hr
        Se = St = Qi*Si*(1-self.COD_rm) / Qe # effluent COD, kg/m3






    def _run_lumped(self):
        '''NOT READY YET.'''


    @property
    def OLR_overall(self):
        '''
        [float] Overall organic loading rate, [kg COD/m3/hr].
        '''
        return self._OLR_overall
    @OLR_overall.setter
    def OLR_overall(self, i):
        if i < 0:
            raise ValueError('`OLR_overall` should be >=0, '
                             f'the input value {i} is outside the range.')
        self._OLR_overall = i

    @property
    def COD_rm(self):
        '''
        [float] COD removal.
        '''
        return self._COD_rm
    @COD_rm.setter
    def COD_rm(self, i):
        if not 0<=i<=1:
            raise ValueError('`COD_rm` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._COD_rm = i

    @property
    def qw(self):
        '''[float] Ratio between the bottom reactor waste flow and the influent.'''
        return self._qw
    @qw.setter
    def qw(self, i):
        if not 0<=i<=1:
            raise ValueError('`qw` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._qw = i

    @property
    def Y(self):
        '''[float] Biomass yield, [kg biomass/kg substrate COD].'''
        return self._Y
    @Y.setter
    def Y(self, i):
        if i < 0:
            raise ValueError('`Y` should be >= 0, '
                             f'the input value {i} is outside the range.')
        self._Y = i

    @property
    def mu_max(self):
        '''[float] Maximum specific growth rate, [/hr].'''
        return self._mu_max
    @mu_max.setter
    def mu_max(self, i):
        if i < 0:
            raise ValueError('`mu_max` should be >= 0, '
                             f'the input value {i} is outside the range.')
        self._mu_max = i

    @property
    def b(self):
        '''[float] Specific endogenous decay coefficient, [/hr].'''
        return self._b
    @b.setter
    def b(self, i):
        if i < 0:
            raise ValueError('`b` should be >= 0, '
                             f'the input value {i} is outside the range.')
        self._b = i

    @property
    def F_Xb(self):
        '''
        [float] Biomass transfer ratio from the bottom reacor to the top reactor,
        should be within [0, 1] (ideal to no retention).
        '''
        return self._F_Xb
    @F_Xb.setter
    def F_Xb(self, i):
        if not 0<=i<=1:
            raise ValueError('`F_Xb` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._F_Xb = i

    @property
    def F_Xt(self):
        '''
        [float] Biomass transfer ratio from the top reacor to the effluent,
        should be within [0, 1] (ideal to no retention).
        '''
        return self._F_Xt
    @F_Xt.setter
    def F_Xt(self, i):
        if not 0<=i<=1:
            raise ValueError('`F_Xt` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._F_Xt = i

    @property
    def BD_dct(self):
        '''
        [dict] Biodegradability of chemicals.
        '''
        return self._BD_dct
    @BD_dct.setter
    def BD_dct(self, i):
        for k, v in i.items():
            if v < 0:
                raise ValueError('Biodegradability should be >=0, '
                                 f'the input value for chemical "{k}" is '
                                 'outside the range.')
        self._BD_dct = self._get_BD_dct(i)