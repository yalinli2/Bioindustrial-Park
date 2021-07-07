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

from biorefineries import cornstover as cs


chems = cs.cornstover.chemicals
# CSL stream is modeled as 50% water, 25% protein, and 25% lactic acid
# get_atom = lambda chemical, element: chemical.atoms.get(element) or 0.

# CSL_atoms = {}
# for i in ('C', 'H', 'O', 'N', 'S'):
#     CSL_atoms[i] = 0.5*get_atom(chems.Water, i)+\
#         0.25*get_atom(chems.Protein, i)+0.25*get_atom(chems.LacticAcid, i)
chems.CSL.formula = 'CH2.8925O1.3275N0.0725S0.00175'