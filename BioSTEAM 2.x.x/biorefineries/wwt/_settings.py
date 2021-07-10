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

References
----------
[1] Hossain et al. Techno-Economic Evaluation of Heat Integrated
Second Generation Bioethanol and Furfural Coproduction.
Biochemical Engineering Journal 2019, 144, 89–103.
https://doi.org/10.1016/j.bej.2019.01.017.

'''

__all__ = (
    'price', 'load_cs_settings',
    'load_sc_settings'
    )

from biorefineries import cornstover as cs
from biorefineries.cornstover._process_settings import (
    price as cs_price,
    load_process_settings as load_cs_settings,
    )

from biorefineries.sugarcane._process_settings import \
    load_process_settings as load_sc_settings

price = cs_price.copy()
price['Caustics'] = cs.caustic.price
price['Wastewater'] = -0.03 # ref [1], negative value for cost from product