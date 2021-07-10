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
Set properties of the chemicals used in the biorefineries.

Chemical data from the lactic acid biorefinery:
https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/blob/master/BioSTEAM%202.x.x/biorefineries/lactic/_chemicals.py

TODO:
    Make separate chems for cornstover, sugarcane, and lipidcane biorefineries
'''

import thermosteam as tmo
from thermosteam.units_of_measure import AbsoluteUnitsOfMeasure as auom

__all__ = (
    'insolubles',
    'get_soluble_ID',
    'create_cs_chemicals',
    )

insolubles = ('Tar', 'Lime', 'CaSO4', 'Ash', 'Lignin', 'Z_mobilis', 'T_reesei',
              'Cellulose', 'Protein', 'Enzyme', 'DenaturedEnzyme', 'WWTsludge')

def get_soluble_ID(chemicals, insolubles):
    return tuple(i.ID for i in chemicals if not i.ID in insolubles)


_cal2joule = auom('cal').conversion_factor('J')

def add_wwt_chemicals(chemicals):
    chems = tmo.Chemicals([i for i in chemicals])
    exist_IDs = [i.ID for i in chems]
    def chemical_database(ID, exist_IDs, phase=None, **data):
        if not ID in exist_IDs:
            chemical = tmo.Chemical(ID, **data)
            if phase:
                chemical.at_state(phase)
                chemical.phase_ref = phase
            chems.append(chemical)
            return chemical

    def chemical_defined(ID, exist_IDs, **data):
        if not ID in exist_IDs:
            chemical = tmo.Chemical.blank(ID, **data)
            chems.append(chemical)
            return chemical

    # CSL stream is modeled as 50% water, 25% protein, and 25% lactic acid,
    # its formula was obtained using the following codes
    # get_atom = lambda chemical, element: chemical.atoms.get(element) or 0.
    # CSL_atoms = {}
    # for i in ('C', 'H', 'O', 'N', 'S'):
    #     CSL_atoms[i] = 0.5*get_atom(chems.Water, i)+\
    #         0.25*get_atom(chems.Protein, i)+0.25*get_atom(chems.LacticAcid, i)
    chems.CSL.formula = 'CH2.8925O1.3275N0.0725S0.00175'

    NaNO3 = chemical_database('NaNO3', exist_IDs, phase='l', Hf=-118756*_cal2joule)
    Na2SO4 = chemical_database('Na2SO4', exist_IDs, phase='l', Hf=-1356380)
    Polymer = chemical_defined('Polymer', exist_IDs,
                               phase='s', MW=1, Hf=0, HHV=0, LHV=0)
    Polymer.Cn.add_model(evaluate=0, name='Constant')

    for i in chems:
        i.default()

    return chems



# %%

# For the cornstover biorefinery
def create_cs_chemicals():
    '''
    Create compiled chemicals for the cornstover biorefinery with the new
    wastewater treatment process.
    '''
    from biorefineries.cornstover import create_chemicals
    cs_chems = create_chemicals()
    new_chems = add_wwt_chemicals(cs_chems)

    return new_chems