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
'''

import thermosteam as tmo
from thermosteam.units_of_measure import AbsoluteUnitsOfMeasure as auom

__all__ = (
    'default_insolubles',
    'get_insoluble_IDs',
    'get_soluble_IDs',
    'create_cs_chemicals',
    )

default_insolubles = (
    # Cornstover biorefinery
    'Lime', 'CaSO4', 'Ash', 'P4O10',
    'Tar', 'Lignin', 'Cellulose', 'Hemicellulose',
    'Protein', 'Enzyme', 'DenaturedEnzyme', 'Z_mobilis', 'T_reesei', 'WWTsludge',
    # Sugarcane biorefinery
    'Yeast', 'CaO', 'Solids', 'Flocculant'
    )

def get_insoluble_IDs(chemicals, insolubles):
    chem_IDs = set([i.ID for i in chemicals])
    new_insolubles = set(insolubles).intersection(chem_IDs)
    return tuple(new_insolubles)

def get_soluble_IDs(chemicals, insolubles):
    return tuple(i.ID for i in chemicals if not i.ID in insolubles)


_cal2joule = auom('cal').conversion_factor('J')

def add_wwt_chemicals(chemicals):
    chems = tmo.Chemicals([i for i in chemicals])
    exist_IDs = [i.ID for i in chems]
    def chemical_database(ID, phase=None, **data):
        if not ID in exist_IDs:
            chemical = tmo.Chemical(ID, **data)
            if phase:
                chemical.at_state(phase)
                chemical.phase_ref = phase
            chems.append(chemical)
            return chemical

    def chemical_defined(ID, **data):
        if not ID in exist_IDs:
            chemical = tmo.Chemical.blank(ID, **data)
            chems.append(chemical)
            return chemical

    chemical_database('NH3', phase='g', Hf=-10963*_cal2joule)
    chemical_database('H2S', phase='g', Hf=-4927*_cal2joule)
    chemical_database('SO2', phase='g')
    chemical_database('NH4OH', search_ID='AmmoniumHydroxide', phase='l', Hf=-336719)
    chemical_database('H2SO4', phase='l')
    chemical_database('HNO3', phase='l', Hf=-41406*_cal2joule)
    chemical_database('NaOH', phase='l')
    chemical_database('NaNO3', phase='l', Hf=-118756*_cal2joule)
    chemical_database('Na2SO4', phase='l', Hf=-1356380)
    chemical_database('CaSO4', phase='s', Hf=-342531*_cal2joule)
    chems.CaSO4.Cn.move_up_model_priority('Lastovka solid', 0)

    chemical_defined('WWTsludge', phase='s',
                     formula='CH1.64O0.39N0.23S0.0035', Hf=-23200.01*_cal2joule)
    chemical_defined('Polymer', phase='s', MW=1, Hf=0, HHV=0, LHV=0)
    chems.Polymer.Cn.add_model(evaluate=0, name='Constant')

    for i in chems:
        i.default()

    return chems


# For the cornstover biorefinery
def create_cs_chemicals():
    '''
    Create compiled chemicals for the cornstover biorefinery with the new
    wastewater treatment process.
    '''
    from biorefineries.cornstover import create_chemicals
    cs_chems = create_chemicals()
    # CSL stream is modeled as 50% water, 25% protein, and 25% lactic acid,
    # its formula was obtained using the following codes
    # get_atom = lambda chemical, element: chemical.atoms.get(element) or 0.
    # CSL_atoms = {}
    # for i in ('C', 'H', 'O', 'N', 'S'):
    #     CSL_atoms[i] = 0.5*get_atom(chems.Water, i)+\
    #         0.25*get_atom(chems.Protein, i)+0.25*get_atom(chems.LacticAcid, i)
    cs_chems.CSL.formula = 'CH2.8925O1.3275N0.0725S0.00175'

    new_chems = add_wwt_chemicals(cs_chems)
    new_chems.compile()

    return new_chems


# For the sugarcane biorefinery
def create_sc_chemicals():
    '''
    Create compiled chemicals for the sugarcane biorefinery with the new
    wastewater treatment process.
    '''
    from biorefineries.sugarcane import create_chemicals
    sc_chems = create_chemicals()
    new_chems = add_wwt_chemicals(sc_chems)
    new_chems.compile()
    return new_chems