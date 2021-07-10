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


import biosteam  as bst
from biorefineries import sugarcane as sc

#!!! Need to enable relative importing
from _chemicals import create_sc_chemicals
from _settings import load_sc_settings
from _wwt_sys import create_wastewater_treatment_system


# %%

# =============================================================================
# Function for making the system
# =============================================================================

sc.sugarcane_sys.

@SystemFactory(
    ID='sugarcane_sys',
    ins=[*sc.create_juicing_system_with_fiber_screener.ins,
         sc.create_ethanol_purification_system.ins[1]], # denaturant
    outs=[create_ethanol_purification_system.outs[0], # ethanol
          dict(ID='vinasse'),
          dict(ID='wastewater'),
          dict(ID='emissions'),
          dict(ID='ash_disposal')]
)
def create_sugarcane_to_ethanol_system(ins, outs):
    s = f.stream
    u = f.unit

    sugarcane, enzyme, H3PO4, lime, polymer, denaturant = ins
    ethanol, vinasse, wastewater, emissions, ash_disposal = outs

    feedstock_handling_sys = create_feedstock_handling_system(
        ins=[sugarcane],
        outs=[''],
        mockup=True,
    )
    juicing_sys = create_juicing_system_with_fiber_screener(
        ins=[feedstock_handling_sys-0, enzyme, H3PO4, lime, polymer],
        mockup=True
    )
    ethanol_production_sys = create_sucrose_to_ethanol_system(
        ins=(juicing_sys-0, denaturant), outs=(ethanol, vinasse),
        mockup=True
    )
    M305 = units.Mixer('M305',
        ins=(juicing_sys-2, *ethanol_production_sys-[2, 3]),
        outs=wastewater
    )

    ### Facilities ###

    BT = units.BoilerTurbogenerator('BT',
        (juicing_sys-1, '', 'boiler_makeup_water', 'natural_gas', '', ''),
        outs=(emissions, 'rejected_water_and_blowdown', ash_disposal),
        boiler_efficiency=0.80,
        turbogenerator_efficiency=0.85
    )
    CT = units.CoolingTower('CT')
    makeup_water_streams = (s.cooling_tower_makeup_water,
                            s.boiler_makeup_water)
    process_water_streams = (s.imbibition_water,
                             s.rvf_wash_water,
                             s.stripping_water,
                             *makeup_water_streams)
    makeup_water = bst.Stream('makeup_water', price=0.000254)
    CWP = units.ChilledWaterPackage('CWP')
    PWC = units.ProcessWaterCenter('PWC',
                                   (bst.Stream(), makeup_water),
                                   (),
                                   None,
                                   makeup_water_streams,
                                   process_water_streams)

    F301 = u.F301
    D303 = u.D303
    HXN = bst.HeatExchangerNetwork('HXN', units=[F301, D303])