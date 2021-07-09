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
from collections import defaultdict
from chemicals.elements import molecular_weight
from thermosteam.reaction import (
    Reaction as Rxn,
    ParallelReaction as PRxn
    )
from biosteam.utils import ExponentialFunctor
from biosteam.units.design_tools.tank_design import (
    mix_tank_purchase_cost_algorithms,
    TankPurchaseCostAlgorithm
    )

__all__ = (
    'get_BD_dct',
    'compute_stream_COD',
    'get_AD_rxns',
    'IC_purchase_cost_algorithms'
    )


# %%

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

    default_IDs = ('O2', 'CO2', 'H2O', 'NH3', 'H2SO4', 'P4O10')
    if chemical.ID in default_IDs:
        return dict.fromkeys(default_IDs, 0)

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
def get_BD_dct(default_BD=0.87, **kwargs):
    BD_dct = defaultdict(lambda: default_BD)

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

    if kwargs: # if user set those biodegradability
        for i in ('Galactose', 'Mannose', 'Arabinose'):
            BD_dct[i] = kwargs.get(i) or BD_dct[i]

    BD_dct['GalactoseOligomer'] = BD_dct['Galactan'] = \
        C6_oligomer_to_monomer * BD_dct['Galactose']
    BD_dct['MannoseOligomer'] = BD_dct['Mannan'] = \
        C6_oligomer_to_monomer * BD_dct['Mannose']
    BD_dct['ArabinoseOligomer'] = BD_dct['Arabinan'] = \
        C5_oligomer_to_monomer * BD_dct['Arabinose']

    if kwargs: # other input biodegradabilities
        BD_dct.update(kwargs)

    return BD_dct


def compute_stream_COD(stream, BD_dct=None):
    r'''
    Compute the chemical oxygen demand (COD) of a given stream in kg-O2/m3
    by summing the COD of each chemical in the stream using:

    .. math::
        COD [\frac{kg}{m^3}] = mol_{chemical} [\frac{kmol}{m^3}] * \frac{g O_2}{mol chemical}

    Can be limited to only biodegradable COD if `BD_dct` is provided.
    '''
    chems = stream.chemicals
    mol = stream.mol
    iCOD = np.array([-get_COD_stoichiometry(i)['O2'] for i in chems])

    # Consider biodegradability
    if BD_dct:
        iCOD *= np.array([BD_dct[i.ID] for i in chems])

    COD = (mol*iCOD).sum() / stream.F_vol * molecular_weight({'O': 2})
    return COD


def get_AD_rxns(stream, BD_dct, X_biogas, X_growth, biomass_ID):
    biomass_MW = getattr(stream.chemicals, biomass_ID).MW
    chems = [i for i in stream.chemicals if i.ID!=biomass_ID]
    if X_biogas+X_growth > 1:
        raise ValueError(f'Sum of `X_biogas` and `X_biogas` is {X_biogas+X_growth}, '
                         'larger than 100%.')

    biogas_rxns = []
    growth_rxns = []
    for i in chems:
        X = BD_dct.get(i.ID)
        if not X: # no entry for the chemical
            if isinstance(BD_dct, defaultdict): # check if have default value
                X = BD_dct[i.ID]
            else:
                continue # assume is not biodegradable

        iX_biogas = X * X_biogas # the amount of chemical used for biogas production
        iX_growth = X * X_growth # the amount of chemical used for cell growth

        biogas_stoyk = get_BMP_stoichiometry(i)
        if biogas_stoyk[i.ID] == 0: # no conversion of this chemical
            continue

        biogas_rxn = Rxn(reaction=biogas_stoyk, reactant=i.ID, X=iX_biogas,
                         check_atomic_balance=True)

        # Cannot check atom balance since the substrate may not have the atom
        growth_rxn = Rxn(f'{i.ID} -> {i.MW/biomass_MW}{biomass_ID}',
                         reactant=i.ID, X=iX_growth, check_atomic_balance=False)

        biogas_rxns.append(biogas_rxn)
        growth_rxns.append(growth_rxn)

    if len(biogas_rxns)>1:
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