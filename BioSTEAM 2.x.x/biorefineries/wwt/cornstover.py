#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021, Yalin Li <yalinli2@illinois.edu>
#
# Part of this module is based on the cornstover biorefinery:
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/cornstover
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


import biosteam as bst
from biosteam import main_flowsheet as F
from biorefineries import (
    cornstover as cs,
    sugarcane as sc
    )

# from biorefineries.wwt import (
#     create_cs_chemicals,
#     cs_price, load_cs_settings,
#     create_wastewater_treatment_system
#     )
# from biorefineries.utils import get_MESP
from _chemicals import create_cs_chemicals
from _settings import cs_price, load_cs_settings
from _wwt_sys import create_wastewater_treatment_system
from utils import get_MESP


# %%

# =============================================================================
# Function to make the system
# =============================================================================

load_cs_settings()
chems = create_cs_chemicals()
bst.settings.set_thermo(chems)


@bst.SystemFactory(
    ID='cornstover_sys',
    ins=[*cs.create_dilute_acid_pretreatment_system.ins,
          dict(ID='denaturant',
               Octane=1,
               price=cs_price['Denaturant'])],
    outs=[dict(ID='ethanol', price=cs_price['Ethanol'])],
)
def create_cs_system(ins, outs, include_blowdown_recycle=False):
    feedstock, denaturant = ins
    ethanol, = outs
    f = bst.main_flowsheet
    s = f.stream
    u = f.unit
    U101 = cs.FeedStockHandling('U101', feedstock)
    U101.cost_items['System'].cost = 0.

    pretreatment_sys = cs.create_dilute_acid_pretreatment_system(
        ins=U101-0,
        mockup=True
    )

    fermentation_sys = cs.create_cellulosic_fermentation_system(
        ins=pretreatment_sys-0,
        mockup=True,
    )

    ethanol_purification_sys = sc.create_ethanol_purification_system(
        ins=[fermentation_sys-1, denaturant],
        outs=[ethanol],
        IDs={'Beer pump': 'P401',
             'Beer column heat exchange': 'H401',
             'Beer column': 'D402',
             'Beer column bottoms product pump': 'P402',
             'Distillation': 'D403',
             'Distillation bottoms product pump': 'P403',
             'Ethanol-denaturant mixer': 'M701',
             'Recycle mixer': 'M402',
             'Heat exchanger to superheat vapor to molecular sieves': 'H402',
             'Molecular sieves': 'U401',
             'Ethanol condenser': 'H403',
             'Ethanol day tank': 'T701',
             'Ethanol day tank pump': 'P701',
             'Denaturant storage': 'T702',
             'Denaturant pump': 'P702',
             'Product tank': 'T703'},
        mockup=True,
    )
    ethanol, stillage, stripper_bottoms_product = ethanol_purification_sys.outs
    recycled_water = bst.Stream(Water=1,
                                T=47+273.15,
                                P=3.9*101325,
                                units='kg/hr')
    S401 = bst.PressureFilter('S401', (stillage, recycled_water))
    if include_blowdown_recycle:
        blowdown_to_wastewater = bst.Stream('blowdown_to_wastewater')
    else:
        blowdown_to_wastewater = None

    create_wastewater_treatment_system(
        ins=[S401-1, pretreatment_sys-1, blowdown_to_wastewater],
        outs=['biogas', 'S603_CHP', 'recycled_water', 'brine'],
        mockup=True,
        IC_method='lumped',
    )

    M501 = bst.Mixer('M501', (u.S603-1, S401-0))

    cs.create_facilities(
        solids_to_boiler=M501-0,
        gas_to_boiler=s.biogas,
        process_water_streams=(s.stripping_water,
                               s.warm_process_water_1,
                               s.warm_process_water_2,
                               s.pretreatment_steam,
                               s.saccharification_water),
        feedstock=feedstock,
        RO_water=u.S604-0,
        recycle_process_water=stripper_bottoms_product,
        blowdown_to_wastewater=blowdown_to_wastewater,
    )


# %%

# =============================================================================
# Create the system
# =============================================================================

flowsheet = bst.Flowsheet('wwt_cornstover')
F.set_flowsheet(flowsheet)
cornstover_sys = create_cs_system(include_blowdown_recycle=False)

u = F.unit
wwt_units = [i for i in u if i.ID[1:3]=='60']
OSBL_units = (*wwt_units, u.CWP, u.CT, u.PWC, u.ADP,
              u.T701, u.T702, u.P701, u.P702, u.M701, u.FT,
              u.CSL_storage, u.DAP_storage, u.BT)
cornstover_tea = cs.create_tea(cornstover_sys, OSBL_units, [u.U101])
ethanol = F.stream.ethanol

# Compare MESP
assert(cornstover_tea.IRR==cs.cornstover_tea.IRR==0.1)
print(f'\n\nIRR = {cornstover_tea.IRR:.0%}')
MESP_old = get_MESP(cs.ethanol, cs.cornstover_tea, 'old cs sys')
MESP_new = get_MESP(ethanol, cornstover_tea, 'new cs sys')


# %%

# =============================================================================
# Util functions for result comparison
# =============================================================================

# Old system
cs_wwt_units = [i for i in cs.cornstover_sys.units if i.ID[1:3]=='60']
old_capex = cs.WWTC.installed_cost / 1e6

old_power_wwt = cs.WWTC.power_utility.rate
old_power_tot = sum(i.power_utility.consumption for i in cs.cornstover_sys.units)
old_power_ratio = old_power_wwt / old_power_tot

old_power_net = sum(i.power_utility.rate for i in cs.cornstover_sys.units)

# New system
new_capexes = {i.ID: i.installed_cost/1e6 for i in wwt_units}
new_capex = sum(i for i in new_capexes.values())

new_powers_wwt = {i.ID: i.power_utility.rate/1e3 for i in wwt_units}
new_power_wwt = sum(i for i in new_powers_wwt.values())
new_power_tot = sum(i.power_utility.consumption for i in cornstover_sys.units)/1e3
new_power_ratio = new_power_wwt / new_power_tot

new_power_net = sum(i.power_utility.rate for i in cornstover_sys.units)/1e3

s = F.stream
# Hf in kJ/hr
net_e = s.ethanol.Hf/3600/1e3 + u.BT.power_utility.rate/1e3
net_e_ratio = net_e/(s.cornstover.Hf/3600/1e3)


# %%

# =============================================================================
# Backup codes for checks
# =============================================================================

# # C/N ratio
# from _utils import get_CN_ratio
# R601_CN = get_CN_ratio(u.R601._inf)

# # Methane production
# from utils import get_digestable_chemicals
# BD_dct = {k.ID: 1. for k in get_digestable_chemicals(chems)}
# u.R601.biodegradability = BD_dct
# u.R601._refresh_rxns(X_biogas=0.86, X_growth=0.05)
# # About 5481, vs. 5681 from Humbird, about 3-4% less,
# # might be due to the different CH4/CO2 ratios in the biogas production reaction,
# # Humbird assumed 51%:49% (1.04) CH4:CO2 on a molar basis, here is about 1
# print(u.R601.outs[0].imass['CH4'])
# get_MESP(ethanol, cornstover_tea, '')