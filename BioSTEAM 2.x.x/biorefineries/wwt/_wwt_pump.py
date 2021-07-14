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
from utils import auom, select_pipe, format_str

__all__ = ('WWTpump',)

_hp_to_kW = auom('hp').conversion_factor('kW')
_lb_to_kg = auom('lb').conversion_factor('kg')
_m_to_ft = auom('m').conversion_factor('ft')
_ft3_to_gal = auom('ft3').conversion_factor('gallon')
_m3_to_gal = auom('m3').conversion_factor('gallon')

class WWTpump(bst.Unit):
    '''
    Generic class for pumps used in wastewater treatment.

    Parameters
    ----------
    reactor_type : str
        Can either be "CSTR" for continuous stirred tank reactor
        or "AF" for anaerobic filter.
    membrane_configuration : str
        Can either be "cross-flow" or "submerged".
    flow_type : str
        Can be "recirculation", "permeate", "retentate", or "lift".

    References
    ----------
    .. [1] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
    Valorization of Dilute Organic Carbon Waste Streams.
    Energy Environ. Sci. 2016, 9 (3), 1102–1112.
    https://doi.org/10.1039/C5EE03715H.
    '''
    _N_ins = 1
    _N_outs = 1

    v = 3 # fluid velocity, [ft/s]
    C = 110 # Hazen- Williams coefficient for stainless steel (SS)

    _default_equipment_lifetime = {'Pump': 15}
    _F_BM_default = bst.Pump._F_BM_default

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 reactor_type='CSTR',
                 membrane_configuration='cross-flow',
                 flow_type='recirculation',
                 Q_mgd=None,
                 **kwargs):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.reactor_type = reactor_type
        self.membrane_configuration = membrane_configuration
        self.flow_type = flow_type
        self.Q_mgd = Q_mgd

        for k, v in kwargs.items():
            setattr(self, k, v)


    def _run(self):
        self.outs[0].copy_like(self.ins[0])


    def _design(self):
        rx_type = format_str(self.reactor_type)
        design_func = getattr(self, f'_design_{rx_type}_{self.flow_type}_pump')
        M_SS_pipe, M_SS_pump = design_func() #!!! need to add in arguments

        D = self.design_results
        D['Total pipe stainless steel'] = M_SS_pipe
        D['Total pump stainless steel'] = M_SS_pump #!!! need to consider pump's lifetime in LCA


    def _cost(self):
        self.power_utility.rate = self.BHP/self.motor_efficiency * _hp_to_kW


    def _design_generic(self, Q_mgd, N_pump, L_s, L_d, H_ts, H_p):
        self.Q_mgd = Q_mgd
        v, C, Q_cfs = self.v, self.C, self.Q_cfs # [ft/s], -, [ft3/s]

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
        M_SS_pipe = 0.29 * (V_s+V_d) * _lb_to_kg

        # Pump SS (for pumps within 300-1000 gpm)
        # http://www.godwinpumps.com/images/uploads/ProductCatalog_Nov_2011_spread2.pdf
        # assume 50% of the product weight is SS
        M_SS_pump = N_pump * (725*0.5)

        return M_SS_pipe, M_SS_pump


    def _design_CSTR_recirculation_pump(self,
                                        Q_R_mgd,
                                        IRR, # internal recirculation ratio
                                        v_GAC, # additional upflow velocity for GAC, [m/hr]
                                        L_train, L_membrane_tank, W_membrane_tank,
                                        N_train, add_GAC=False):
        # Total IR pumping flow rate, [mgd]
        Q_IR_initial_mgd = max(0, (self.Q_mgd*IRR)-Q_R_mgd)

        if add_GAC: # to achieve adequate upflow velocity for GAC
            v_GAC *= _m_to_ft # [ft/hr]
            A_tank = L_membrane_tank * W_membrane_tank # [ft2]
            Q_upflow_req = v_GAC * A_tank * N_train # [ft3/hr]
            Q_upflow_req *= 24 * _ft3_to_gal/1e6 # [mgd]
            Q_additional_IR = max(0, Q_upflow_req-(self.Q_mgd+Q_IR_initial_mgd))
            Q_IR_mgd = Q_IR_initial_mgd+Q_additional_IR;
        else:
            Q_IR_mgd = Q_IR_initial_mgd

        M_SS_IR_pipe, M_SS_IR_pump = self._design_generic(
            Q_mgd=Q_IR_mgd,
            N_pump=1,
            L_s=0., # ignore suction side
            L_d=L_train, # pipe length per train
            H_ts=5., # H_ds_IR (5) - H_ss_IR (0)
            H_p=0. # no pressure
            )

        return M_SS_IR_pipe, M_SS_IR_pump


    def _design_CSTR_retentate_pump(self, Q_R_mgd,
                                    N_unit, # number of membrane units
                                    D_train):
        M_SS_IR_pipe, M_SS_IR_pump = self._design_generic(
            Q_mgd=Q_R_mgd,
            N_pump=N_unit,
            L_s=100, # pipe length per module
            L_d=30, # pipe length per filter (same as discharge side of lifting pump)
            H_ts=0., # H_ds_IR (D_train) - H_ss_IR (D_train)
            H_p=0. # no pressure
            )

        return M_SS_IR_pipe, M_SS_IR_pump


    def _design_cross_flow_permeate_pump(self,
                                         N_LU, # number of large membrane units
                                         TMP, # transmembrane pressure, [psi]
                                         D_train,
                                         include_aerobic_filteration=False):
        H_ts_PERM = D_train if include_aerobic_filteration else 0

        M_SS_IR_pipe, M_SS_IR_pump = self._design_generic(
            Q_mgd=self.Q_mgd,
            N_pump=N_LU,
            L_s=20, # based on a 30-module unit with a total length of 6 m, [ft]
            L_d_R=10*N_LU, # based on a 30-module unit with a total width of 1.6 m and extra space, [ft]
            H_ts=H_ts_PERM, #  H_ds_PERM (D_train) - H_ss_PERM (0 or D_train)
            H_p_R=TMP*2.31 # TMP in water head, [ft]
            )

        # # factor = 2.31 calculated by
        # factor = auom('psi').conversion_factor('Pa') # Pa is kg/m/s2, now in [Pa]
        # factor /= 9.81 # divided by the standard gravity in m/s2, now in [kg/m2]
        # factor /= 1e3 # divided by water's density in kg/m3, now in [m]
        # factor *= auom('m').conversion_factor('ft') # m to ft

        return M_SS_IR_pipe, M_SS_IR_pump

    def _design_cross_flow_lift_pump(self, N_train):
        M_SS_IR_pipe, M_SS_IR_pump = self._design_generic(
            Q_mgd=self.Q_mgd,
            N_pump=2, # one for waste sludge and one for return sludge/clarifier
            L_s=150, # length of suction pipe per filter, [ft]
            L_d_R=30, # pipe length per filter
            H_ts=12, # H_ds_LIFT (12) - H_ss_LIFT (0)
            H_p_R=0. # TMP in water head, [ft]
            )




    @property
    def reactor_type(self):
        '''
        [str] Can either be "CSTR" for continuous stirred tank reactor
        or "AF" for anaerobic filter.
        '''
        return self._reactor_type
    @reactor_type.setter
    def reactor_type(self, i):
        if not i.upper() in ('CSTR', 'AF'):
            raise ValueError('`reactor_type` can only be "CSTR", or "AF", '
                             f'not "{i}".')
        self._reactor_type = i.upper()

    @property
    def membrane_configuration(self):
        '''[str] Can either be "cross-flow" or "submerged".'''
        return self._membrane_configuration
    @membrane_configuration.setter
    def membrane_configuration(self, i):
        i = 'cross-flow' if i.lower() in ('cross flow', 'crossflow') else i
        if not i.lower() in ('cross-flow', 'submerged'):
            raise ValueError('`membrane_configuration` can only be "cross-flow", '
                             f'or "submerged", not "{i}".')
        self._membrane_configuration = i.lower()

    @property
    def flow_type(self):
        '''
        [str] Can either be "recirculation", "permeate", or "retentate".
        '''
        return self._flow_type
    @flow_type.setter
    def flow_type(self, i):
        if not i.lower() in ('recirculation', 'permeate', 'retentate', 'lift'):
            raise ValueError('`flow_type` can only be "recirculation", '
                             f'"permeate", "retentate", or "lift", not "{i}".')
        self._flow_type = i.lower()

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
        return (self.TDH*self.Q_gpm)/3960/self.brake_efficiency

    @property
    def Q_mgd(self):
        '''
        [float] Volumetric flow rate in million gallon per day, [mgd].
        Will return total volumetric flow through the unit if not provided.
        '''
        if self._Q_mgd:
            return self._Q_mgd
        return self.F_vol_in*_m3_to_gal*24/1e6
    @Q_mgd.setter
    def Q_mgd(self, i):
        self._Q_mgd = i

    @property
    def Q_cfs(self):
        '''[float] Volumetric flow rate in cubic feet per second, [cfs].'''
        return self.Q_mgd*1e6/24/60/60/_ft3_to_gal

    @property
    def Q_gpm(self):
        '''[float] Volumetric flow rate in gallon per minute, [gpm].'''
        return self.Q_mgd*1e6/24/60

    @property
    def brake_efficiency(self):
        '''[float] Brake efficiency.'''
        return brake_efficiency(self.Q_gpm)

    @property
    def motor_efficiency(self):
        '''[float] Motor efficiency.'''
        return motor_efficiency(self.Q_gpm, self.BHP)