#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021, Yalin Li <yalinli2@illinois.edu>
#
# Part of this module is based on the sugarcane biorefinery:
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/sugarcane
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


import biosteam  as bst
from biosteam import main_flowsheet as F
from biorefineries import sugarcane as sc

#!!! Need to enable relative importing
from _chemicals import create_sc_chemicals
from _utils import kph_to_tpd, get_MESP
from _settings import new_price, load_sc_settings
from _wwt_sys import create_wastewater_treatment_system


# %%

# =============================================================================
# Function to make the system
# =============================================================================

load_sc_settings()
chems = create_sc_chemicals()
bst.settings.set_thermo(chems)

@bst.SystemFactory(
    ID='sugarcane_sys',
    ins=[*sc.create_juicing_system_with_fiber_screener.ins,
         sc.create_ethanol_purification_system.ins[1]], # denaturant
    outs=[sc.create_ethanol_purification_system.outs[0], # ethanol
          dict(ID='vent_R602'),
          dict(ID='brine'),
          dict(ID='emissions'),
          dict(ID='ash_disposal')]

)
def create_sc_system(ins, outs):
    s = F.stream
    u = F.unit

    sugarcane, enzyme, H3PO4, lime, polymer, denaturant = ins
    ethanol, vent_R602, brine, emissions, ash_disposal = outs

    feedstock_handling_sys = sc.create_feedstock_handling_system(
        ins=[sugarcane],
        outs=[''],
        mockup=True,
    )

    juicing_sys = sc.create_juicing_system_with_fiber_screener(
        ins=[feedstock_handling_sys-0, enzyme, H3PO4, lime, polymer],
        mockup=True
    )

    ethanol_production_sys = sc.create_sucrose_to_ethanol_system(
        ins=(juicing_sys-0, denaturant), outs=(ethanol, 'vinasse'),
        mockup=True
    )

    M305 = bst.units.Mixer(
        'M305',
        ins=(juicing_sys-2, *ethanol_production_sys-[2, 3]),
        outs='wastewater'
    )

    ### Wastewater treatment ###
    create_wastewater_treatment_system(
        ins=[ethanol_production_sys-1, M305-0],
        outs=['biogas', vent_R602, 'S604_CHP', 'recycled_water', brine],
        mockup=True,
        IC_method='lumped',
        dry_flow_tpd=kph_to_tpd(s.sugarcane),
        need_ammonia=False
    )

    ### Facilities ###
    bst.units.BoilerTurbogenerator('BT',
        (juicing_sys-1, u.R601-0, 'boiler_makeup_water', 'natural_gas', '', ''),
        outs=(emissions, 'rejected_water_and_blowdown', ash_disposal),
        boiler_efficiency=0.80,
        turbogenerator_efficiency=0.85)
    bst.units.CoolingTower('CT')
    bst.units.ChilledWaterPackage('CWP')

    makeup_water_streams = (s.cooling_tower_makeup_water,
                            s.boiler_makeup_water)
    process_water_streams = (s.imbibition_water,
                             s.rvf_wash_water,
                             s.stripping_water,
                             *makeup_water_streams)
    makeup_water = bst.Stream('makeup_water', price=0.000254)
    bst.units.ProcessWaterCenter('PWC',
                                 (u.S605-0, # recycled wastewater from reverse osmosis
                                 makeup_water),
                                 (),
                                 None,
                                 makeup_water_streams,
                                 process_water_streams)

    plant_air = bst.Stream('plant_air', N2=83333, units='kg/hr')
    ADP = bst.facilities.AirDistributionPackage('ADP', plant_air)
    @ADP.add_specification(run=True)
    def adjust_plant_air():
        plant_air.imass['N2'] = 0.8 * feedstock_handling_sys.ins[0].F_mass

    F301 = u.F301
    D303 = u.D303
    bst.HeatExchangerNetwork('HXN', units=[F301, D303])


# %%

# =============================================================================
# Create the system
# =============================================================================

flowsheet = bst.Flowsheet('wwt_sugarcane')
F.set_flowsheet(flowsheet)
sugarcane_sys = create_sc_system()
u = F.unit

sugarcane_sys.simulate()
sugarcane_tea = sc.create_tea(sugarcane_sys)
ethanol = F.stream.ethanol

# Compare MESP
original_IRR = 0.1267
sugarcane_tea.IRR = sc.sugarcane_tea.IRR = original_IRR
print(f'\n\nIRR = {original_IRR:.0%}')
sc.wastewater.price = 0.
MESP_old = get_MESP(sc.ethanol, sc.sugarcane_tea, 'old sc sys w/o ww cost')
sc.wastewater.price = new_price['Wastewater']
MESP_old = get_MESP(sc.ethanol, sc.sugarcane_tea, 'old sc sys w ww cost')
MESP_new = get_MESP(ethanol, sugarcane_tea, 'new sc sys')

sugarcane_tea.IRR = sc.sugarcane_tea.IRR = 0.1
assert(sugarcane_tea.IRR==sc.sugarcane_tea.IRR)
print(f'\n\nIRR = {sugarcane_tea.IRR:.0%}')
sc.wastewater.price = 0.
MESP_old = get_MESP(sc.ethanol, sc.sugarcane_tea, 'old sc sys w/o ww cost')
sc.wastewater.price = new_price['Wastewater']
MESP_old = get_MESP(sc.ethanol, sc.sugarcane_tea, 'old sc sys w ww cost')
MESP_new = get_MESP(ethanol, sugarcane_tea, 'new sc sys')

# Ratios for IC design
# from _utils import get_CN_ratio
# R601_CN = get_CN_ratio(u.R601._inf)