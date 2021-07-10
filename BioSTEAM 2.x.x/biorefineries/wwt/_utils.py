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
Util functions
'''

import numpy as np
from chemicals.elements import molecular_weight
from thermosteam.reaction import (
    Reaction as Rxn,
    ParallelReaction as PRxn
    )
from biosteam.utils import (
    ExponentialFunctor,
    remove_undefined_chemicals,
    default_chemical_dict
    )
from biosteam.units.design_tools.tank_design import (
    mix_tank_purchase_cost_algorithms,
    TankPurchaseCostAlgorithm
    )
from _chemicals import default_insolubles

__all__ = (
    'get_BD_dct',
    'compute_stream_COD',
    'get_digestion_rxns',
    'IC_purchase_cost_algorithms',
    'get_MB_split',
    # 'remove_undefined_chemicals',
    )


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
        -> nCO_2 + \frac{a-3c-2d}{2}H_2O + cNH_3 + dH_2SO_4 + \frac{e}{4}P_4O_10
    '''
    Xs = nC, nH, nO, nN, nS, nP = get_CHONSP(chemical)

    excluded = ('O2', 'CO2', 'H2O', 'NH3', 'H2SO4', 'P4O10',
                *default_insolubles)
    if chemical.ID in excluded:
        return dict.fromkeys(excluded, 0)

    dct = {
        chemical.ID: -1. if sum(Xs)!=0 else 0.,
        'O2': -(nC+0.25*nH-0.5*nO-0.75*nN+1.5*nS+1.25*nP),
        'CO2': nC,
        'H2O': 0.5*nH-1.5*nN-nS, # assume one water reacts with SO3 to H2SO4
        'NH3': nN,
        'H2SO4': nS,
        'P4O10': 0.25*nP
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

    default_IDs = ('H2O', 'CH4', 'CO2', 'NH3', 'H2S')
    if chemical.ID in default_IDs:
        return dict.fromkeys(default_IDs, 0)

    dct = {
        chemical.ID: -1. if sum(Xs)!=0 else 0.,
        'H2O': -(nC-0.25*nH-0.5*nO+0.75*nN+0.5*nS),
        'CH4': 0.5*nC+0.125*nH-0.25*nO-0.375*nN-0.25*nS,
        'CO2': 0.5*nC-0.125*nH+0.25*nO+0.375*nN+0.25*nS,
        'NH3': nN,
        'H2S': nS
        }

    return dct


# Biodegradability, 0.87 from glucose (treated as the maximum value)
def get_BD_dct(chemicals, default_BD=0.87, **kwargs):
    BD_dct = dict.fromkeys([i.ID for i in chemicals], default_BD)

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

    # Assume biodegradabilities based on glucose and xylose
    BD_dct['Galactose'] = BD_dct['Mannose'] = \
        BD_dct['Arabinose']/BD_dct['Xylose'] * BD_dct['Glucose']
    BD_dct['GalactoseOligomer'] = BD_dct['Galactan'] = \
        BD_dct['MannoseOligomer'] = BD_dct['Mannan'] = \
             BD_dct['Glucan']

    BD_dct['ArabinoseOligomer'] = BD_dct['Arabinan'] = BD_dct['Xylan']

    if kwargs: # other input biodegradabilities
        BD_dct.update(kwargs)

    return BD_dct


def compute_stream_COD(stream):
    r'''
    Compute the chemical oxygen demand (COD) of a given stream in kg-O2/m3
    by summing the COD of each chemical in the stream using:

    .. math::
        COD [\frac{kg}{m^3}] = mol_{chemical} [\frac{kmol}{m^3}] * \frac{g O_2}{mol chemical}
    '''
    chems = stream.chemicals
    mol = stream.mol
    iCOD = np.array([-get_COD_stoichiometry(i)['O2'] for i in chems])

    COD = (mol*iCOD).sum() / stream.F_vol * molecular_weight({'O': 2})
    return COD


def get_digestion_rxns(stream, BD, X_biogas, X_growth, biomass_ID):
    biomass_MW = getattr(stream.chemicals, biomass_ID).MW
    chems = [i for i in stream.chemicals if i.ID!=biomass_ID]
    if isinstance(BD, float):
        BD = dict.fromkeys([i.ID for i in chems], BD)

    if X_biogas+X_growth > 1:
        raise ValueError(f'Sum of `X_biogas` and `X_biogas` is {X_biogas+X_growth}, '
                         'larger than 100%.')

    biogas_rxns = []
    growth_rxns = []
    for i in chems:
        X = BD.get(i.ID)
        if not X:
            continue # assume no entry means not biodegradable

        biogas_stoyk = get_BMP_stoichiometry(i)
        if biogas_stoyk[i.ID] == 0: # no conversion of this chemical
            continue

        iX_biogas = X * X_biogas # the amount of chemical used for biogas production
        iX_growth = X * X_growth # the amount of chemical used for cell growth

        if iX_biogas:
            biogas_rxn = Rxn(reaction=biogas_stoyk, reactant=i.ID, X=iX_biogas,
                             check_atomic_balance=True)
            biogas_rxns.append(biogas_rxn)

        if iX_growth:
        # Cannot check atom balance since the substrate may not have the atom
        #!!! Maybe balance this with CSL, DAP, and some other chemicals
            growth_rxn = Rxn(f'{i.ID} -> {i.MW/biomass_MW}{biomass_ID}',
                             reactant=i.ID, X=iX_growth,
                             check_atomic_balance=False)


            growth_rxns.append(growth_rxn)

    if len(biogas_rxns)+len(growth_rxns)>1:
        return PRxn(biogas_rxns+growth_rxns)

    return []


# Tank cost algorithms
IC_purchase_cost_algorithms = mix_tank_purchase_cost_algorithms.copy()
conventional = IC_purchase_cost_algorithms['Conventional']
#!!! Need to check if the cost correlation still holds for the ranges beyond
ic = TankPurchaseCostAlgorithm(
    ExponentialFunctor(A=conventional.f_Cp.A,
                       n=conventional.f_Cp.n),
    V_min=np.pi/4*1.5**2*16, # 1.5 and 16 are the lower bounds of the width and height ranges in ref [1]
    V_max=np.pi/4*12**2*25, # 12 and 25 are the lower bounds of the width and height ranges in ref [1]
    V_units='m^3',
    CE=conventional.CE,
    material='Stainless steel')

IC_purchase_cost_algorithms['IC'] = ic


# %%

# Split for the membrane bioreactor
def get_MB_split(chemicals, split_dct=None):
    # Copied from the cornstover biorefinery,
    # which is based on the 2011 NREL report (Humbird et al.)
    split = dict(
        Water=0.1454,
        Glycerol=0.125,
        LacticAcid=0.145,
        SuccinicAcid=0.125,
        HNO3=0.1454,
        Denaturant=0.125,
        DAP=0.1454,
        AmmoniumAcetate=0.145,
        AmmoniumSulfate=0.1454,
        H2SO4=0.1454,
        NaNO3=0.1454,
        Oil=0.125,
        N2=0.1351,
        NH3=0.1579,
        O2=0.15,
        CO2=0.1364,
        Xylose=0.25,
        Sucrose=0.125,
        Mannose=0.125,
        Galactose=0.125,
        Arabinose=0.125,
        Extract=0.145,
        NaOH=0.1454,
        SolubleLignin=0.145,
        GlucoseOligomer=0.1429,
        GalactoseOligomer=0.1429,
        MannoseOligomer=0.1429,
        XyloseOligomer=0.1429,
        ArabinoseOligomer=0.1429,
        Xylitol=0.125,
        Cellobiose=0.125,
        Cellulase=0.145
        )
    remove_undefined_chemicals(split, chemicals)
    default_chemical_dict(split, chemicals, 0.15, 0.125, 0.145)

    if split_dct is not None:
        split.update(split_dct)

    return split