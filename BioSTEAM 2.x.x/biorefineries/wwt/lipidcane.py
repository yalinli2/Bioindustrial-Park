#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021, Yalin Li <yalinli2@illinois.edu>
#
# Part of this module is based on the lipidcane biorefinery:
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/lipidcane
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import biosteam  as bst
from biosteam import main_flowsheet as F
from biorefineries import (
    lipidcane as lc,
    sugarcane as sc
    )

#!!! Need to enable relative importing
from _chemicals import create_lc_chemicals
from _utils import kph_to_tpd, get_MESP
from _settings import lc_price, new_price, load_lc_settings
from _wwt_sys import create_wastewater_treatment_system


# %%

# =============================================================================
# Function to make the system
# =============================================================================

load_lc_settings()
chems = create_lc_chemicals()
bst.settings.set_thermo(chems)

@bst.SystemFactory(
    ID='lipidcane_sys',
    ins=[*lc.create_juicing_and_lipid_extraction_system.ins,
         sc.create_sucrose_to_ethanol_system.ins[1]],
    outs=[dict(ID='ethanol', price=lc_price['Ethanol']),
          dict(ID='biodiesel', price=lc_price['Biodiesel']),
          dict(ID='crude_glycerol', price=lc_price['Crude glycerol']),
          dict(ID='vinasse'),
          dict(ID='wastewater'),
          dict(ID='emissions'),
          dict(ID='ash_disposal')]
)
def create_lc_system(ins, outs):
    s = F.stream
    u = F.unit

    lipidcane, enzyme, H3PO4, lime, polymer, denaturant = ins
    ethanol, biodiesel, crude_glycerol, vinasse, wastewater, emissions, ash_disposal = outs

    feedstock_handling_sys = sc.create_feedstock_handling_system(
        ins=lipidcane,
        mockup=True,
    )

    ### Oil and juice separation ###
    juicing_and_lipid_extraction_sys = lc.create_juicing_and_lipid_extraction_system(
        ins=[feedstock_handling_sys-0, enzyme, H3PO4, lime, polymer],
        mockup=True,
    )

    ### Ethanol section ###
    ethanol_production_sys = sc.create_sucrose_to_ethanol_system(
        ins=[juicing_and_lipid_extraction_sys-0, denaturant],
        outs=[ethanol, vinasse],
        mockup=True,
    )

    ### Biodiesel section ###
    # Fresh degummed oil
    oil = juicing_and_lipid_extraction_sys-1
    lc.create_transesterification_and_biodiesel_separation_system(
        ins=oil,
        outs=[biodiesel, crude_glycerol],
        mockup=True,
    )

    ### Wastewater treatment ###
    create_wastewater_treatment_system(
        ins=[vinasse, wastewater],
        mockup=True,
        IC_method='lumped',
        dry_flow_tpd=kph_to_tpd(s.lipidcane),
        need_ammonia=False
    )

    ### Facilities ###
    bst.Mixer('M305',
        ins=(juicing_and_lipid_extraction_sys-4,
             juicing_and_lipid_extraction_sys-3,
             *ethanol_production_sys-[1, 2, 3]),
        outs=wastewater,
    )

    # Burn bagasse from conveyor belt
    bst.units.BoilerTurbogenerator('BT',
                                   (juicing_and_lipid_extraction_sys-2, u.R601-0,
                                    'boiler_makeup_water', 'natural_gas', '', ''),
                                   (emissions, 'rejected_water_and_blowdown', ash_disposal),
                                   boiler_efficiency=0.80,
                                   turbogenerator_efficiency=0.85)

    bst.units.CoolingTower('CT')
    makeup_water_streams = (s.cooling_tower_makeup_water,
                            s.boiler_makeup_water)

    process_water_streams = (s.imbibition_water,
                             s.biodiesel_wash_water,
                             s.oil_wash_water,
                             s.rvf_wash_water,
                             s.stripping_water,
                             *makeup_water_streams)

    makeup_water = bst.Stream('makeup_water', price=0.000254)

    bst.units.ChilledWaterPackage('CWP')
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

flowsheet = bst.Flowsheet('wwt_lipidcane')
F.set_flowsheet(flowsheet)
lipidcane_sys = create_lc_system()
u = F.unit

lipidcane_sys.simulate()
lipidcane_tea = lc.create_tea(lipidcane_sys)
ethanol = F.stream.ethanol

# Compare MESP
lipidcane_tea.IRR = lc.lipidcane_tea.IRR
assert(lipidcane_tea.IRR==lc.lipidcane_tea.IRR)
print(f'\n\nIRR = {lipidcane_tea.IRR:.0%}')
MESP_old = get_MESP(lc.ethanol, lc.lipidcane_tea, 'old lc sys w/o ww cost')
lc.wastewater.price = new_price['Wastewater']
MESP_old = get_MESP(lc.ethanol, lc.lipidcane_tea, 'old lc sys w ww cost')
MESP_new = get_MESP(ethanol, lipidcane_tea, 'new lc sys')

lipidcane_tea.IRR = lc.lipidcane_tea.IRR = 0.1
assert(lipidcane_tea.IRR==lc.lipidcane_tea.IRR)
print(f'\n\nIRR = {lipidcane_tea.IRR:.0%}')
lc.wastewater.price = 0.
MESP_old = get_MESP(lc.ethanol, lc.lipidcane_tea, 'old lc sys w/o ww cost')
lc.wastewater.price = new_price['Wastewater']
MESP_old = get_MESP(lc.ethanol, lc.lipidcane_tea, 'old lc sys w ww cost')
MESP_new = get_MESP(ethanol, lipidcane_tea, 'new lc sys')

# Ratios for IC design
# from _utils import get_CN_ratio
# R601_CN = get_CN_ratio(u.R601._inf)