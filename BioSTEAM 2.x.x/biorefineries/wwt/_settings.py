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
Settings for cornstover (cs), sugarcane (sc), and lipidcane (lc) biorefineries.
'''

__all__ = (
    'cs_price', 'load_cs_settings',
    'load_sc_settings'
    )

from biorefineries.cornstover._process_settings import (
    price as cs_price,
    load_process_settings as load_cs_settings,
    )

from biorefineries.sugarcane._process_settings import \
    load_process_settings as load_sc_settings