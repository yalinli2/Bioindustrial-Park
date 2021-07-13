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

import math
import biosteam as bst
from biosteam.design_tools.mechanical import brake_efficiency, motor_efficiency

#!!! Need to enable relative importing
from utils import auom, select_pipe

__all__ = ('AnMBRpump',)


class AnMBRpump(bst.Unit):
    '''
    Generic class for pumps used in water/wastewater treatment.

    References
    ----------
    .. [1] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
    Valorization of Dilute Organic Carbon Waste Streams.
    Energy Environ. Sci. 2016, 9 (3), 1102–1112.
    https://doi.org/10.1039/C5EE03715H.
    '''
    v = 3 # fluid velocity
    C = 110 # Hazen- Williams coefficient for stainless steel (SS)

    _default_equipment_lifetime = {'Pump': 15}
    _F_BM_default = bst.Pump._F_BM_default

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 N_pump=1, L_s=0., L_d=0., H_ts=0., H_p=0.,
                 **kwargs):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.N_pump = N_pump

        for k, v in kwargs.items():
            setattr(self, k, v)


    def _run(self):
        N_pump, v, L_s, L_d, C = \
            self.N_pump, self.v, self.L_s, self.L_d, self.C

        Q_cms = self.ins[0].F_vol # volumetric flow, [m3/s]
        Q_cfs = Q_cms * auom('m3').conversion_factor('ft3') # [ft3/s]

        ### Suction side ###
        # Suction pipe (permeate header) dimensions
        OD_s, t_s, ID_s = select_pipe(Q_cfs/N_pump, v) # [in]

        # Suction friction head, [ft]
        self._H_sf = 3.02 * L_s * (v**1.85) * (C**(-1.85)) * ((ID_s/12)**(-1.17))

        ### Discharge side ###
        # Discharge pipe (permeate collector) dimensions
        OD_d, t_d, ID_d = select_pipe(Q_cfs, v)

        # Discharge friction head, [ft]
        self._H_df = 3.02 * L_d * (v**1.85) * (C**(-1.85)) * ((ID_d/12)**(-1.17))

        ### Material usage ###
        # Pipe SS, assume stainless steel, density = 0.29 lbs/in3
        # SS volume for suction, [in3]
        V_s = N_pump * math.pi/4*((OD_s)**2-(ID_s)**2) * (L_s*12)
        # SS volume for discharge, [in3]
        V_d =math.pi/4*((OD_d)**2-(ID_d)**2) * (L_d*12)
        # Total SS mass, [kg]
        M_SS_pipe = 0.29 * (V_s+V_d) * auom('lb').conversion_factor('kg')

        # Pump SS (for pumps within 300-1000 gpm)
        # http://www.godwinpumps.com/images/uploads/ProductCatalog_Nov_2011_spread2.pdf
        # assume 50% of the product weight is SS
        M_SS_pump = N_pump * (725*0.5)

        D = self.design_results
        D['Total pipe stainless steel'] = M_SS_pipe
        D['Total pump stainless steel'] = M_SS_pump


    def _cost(self):
        self.power_utility.rate = self.BHP/self.motor_efficiency * \
            auom('hp').conversion_factor('kW')


    @property
    def H_sf(self):
        '''[float] Suction friction head, [ft].'''
        return self._H_sf

    @property
    def H_df(self):
        '''[float] Discharge friction head, [ft].'''
        return self._H_df

    @property
    def TDH(self):
        '''[float] Total dynamic head, [ft].'''
        return self.H_ts+self.H_sf+self.H_df+self.H_p

    @property
    def BHP(self):
        '''[float] Brake horsepower, [hp].'''
        return (self.TDH*self.Q_gpm) / 3960 / self.brake_efficiency

    @property
    def Q_gpm(self):
        '''[float] Volumetric flow rate in gallon per minute, [gpm].'''
        return self.ins[0].F_vol * auom('m3').conversion_factor('gallon') * 60

    @property
    def brake_efficiency(self):
        '''[float] Brake efficiency.'''
        return brake_efficiency(self.Q_gpm)

    @property
    def motor_efficiency(self):
        '''[float] Motor efficiency.'''
        return motor_efficiency(self.Q_gpm, self.BHP)