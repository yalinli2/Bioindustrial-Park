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

TODO:
    Make separate chems for cornstover, sugarcane, and lipidcane biorefineries
'''

# %%

# =============================================================================
# For the cornstover biorefinery
# =============================================================================

import thermosteam as tmo
from biorefineries.cornstover import create_chemicals

def create_cs_chemicals():
    '''
    Create compiled chemicals for the cornstover biorefinery with the new
    wastewater treatment process.
    '''
    cs_chems = create_chemicals()
    chems = tmo.Chemicals([i for i in cs_chems])