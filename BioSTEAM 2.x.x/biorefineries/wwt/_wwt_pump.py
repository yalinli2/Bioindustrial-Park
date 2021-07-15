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
_ft_to_m = auom('ft').conversion_factor('m')
_ft3_to_gal = auom('ft3').conversion_factor('gallon')
_m3_to_gal = auom('m3').conversion_factor('gallon')

class WWTpump(bst.Unit):
    '''
    Generic class for pumps used in wastewater treatment.

    Parameters
    ----------
    pump_type : str
        The type of the pump that determines the design algorithms to use.
        The following combination is valid:
            - "cross-flow_permeate"
    Q_mgd : float
        Volumetric flow rate in million gallon per day.
        Will return total volumetric flow through the unit if not provided.
    inputs : dct
        Additional inputs that will be passed to the corresponding design alogrithm.
        Check the document for the design alogrithm for the specific input requirements.

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

    _valid_pump_types = (
        'cross-flow_permeate',
        )


    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 pump_type, Q_mgd=None, **inputs):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.pump_type = pump_type
        self.Q_mgd = Q_mgd
        self.inputs = inputs


    def _run(self):
        self.outs[0].copy_like(self.ins[0])


    def _design(self):
        pump_type = format_str(self.pump_type)
        design_func = getattr(self, f'design_{pump_type}')
        M_SS_pipe, M_SS_pump = design_func() #!!! need to add in arguments

        D = self.design_results
        D['Pipe stainless steel'] = M_SS_pipe
        D['Pump stainless steel'] = M_SS_pump #!!! need to consider pump's lifetime in LCA


    def _cost(self):
        self.power_utility.rate = self.BHP/self.motor_efficiency * _hp_to_kW


    # Generic algorithms that will be called by all design functions
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


    def design_cross_flow_permeate(self):
        '''
        Design pump for the permeate stream of cross-flow membrane configuration.

        Parameters
        ----------
        N_LU : int
            Number of large membrane units.
        D_tank: float
            Depth of the membrane tank, [ft].
        TMP : float
            Transmembrane pressure, [psi].
        include_aerobic_filter : bool
            Whether aerobic filter is included in the reactor design.

        '''
        N_LU, D_tank, TMP, include_aerobic_filter = self.inputs

        H_ts_PERM = D_tank if include_aerobic_filter else 0

        M_SS_IR_pipe, M_SS_IR_pump = self.design_generic(
            Q_mgd=self.Q_mgd,
            N_pump=N_LU,
            L_s=20, # based on a 30-module unit with a total length of 6 m, [ft]
            L_d_R=10*N_LU, # based on a 30-module unit with a total width of 1.6 m and extra space, [ft]
            H_ts=H_ts_PERM, #  H_ds_PERM (D_tank) - H_ss_PERM (0 or D_tank)
            H_p_R=TMP*2.31 # TMP in water head, [ft]
            )

        # # factor = 2.31 calculated by
        # factor = auom('psi').conversion_factor('Pa') # Pa is kg/m/s2, now in [Pa]
        # factor /= 9.81 # divided by the standard gravity in m/s2, now in [kg/m2]
        # factor /= 1e3 # divided by water's density in kg/m3, now in [m]
        # factor *= auom('m').conversion_factor('ft') # m to ft

        return M_SS_IR_pipe, M_SS_IR_pump


    def design_CSTR_recirculation(self,
                                        Q_R_mgd,
                                        IRR, # internal recirculation ratio
                                        v_GAC, # additional upflow velocity for GAC, [m/hr]
                                        L_train, L_membrane_tank, W_membrane_tank,
                                        N_train, add_GAC=False):
        # Total IR pumping flow rate, [mgd]
        Q_IR_initial_mgd = max(0, (self.Q_mgd*IRR)-Q_R_mgd)

        if add_GAC: # to achieve adequate upflow velocity for GAC
            v_GAC /= _ft_to_m # [ft/hr]
            A_tank = L_membrane_tank * W_membrane_tank # [ft2]
            Q_upflow_req = v_GAC * A_tank * N_train # [ft3/hr]
            Q_upflow_req *= 24 * _ft3_to_gal/1e6 # [mgd]
            Q_additional_IR = max(0, Q_upflow_req-(self.Q_mgd+Q_IR_initial_mgd))
            Q_IR_mgd = Q_IR_initial_mgd+Q_additional_IR;
        else:
            Q_IR_mgd = Q_IR_initial_mgd

        M_SS_IR_pipe, M_SS_IR_pump = self.design_generic(
            Q_mgd=Q_IR_mgd,
            N_pump=1,
            L_s=0., # ignore suction side
            L_d=L_train, # pipe length per train
            H_ts=5., # H_ds_IR (5) - H_ss_IR (0)
            H_p=0. # no pressure
            )

        return M_SS_IR_pipe, M_SS_IR_pump


    def design_CSTR_retentate(self, Q_R_mgd,
                                    N_unit, # number of membrane units
                                    D_tank):
        M_SS_IR_pipe, M_SS_IR_pump = self.design_generic(
            Q_mgd=Q_R_mgd,
            N_pump=N_unit,
            L_s=100, # pipe length per module
            L_d=30, # pipe length per filter (same as discharge side of lifting pump)
            H_ts=0., # H_ds_IR (D_tank) - H_ss_IR (D_tank)
            H_p=0. # no pressure
            )

        return M_SS_IR_pipe, M_SS_IR_pump


    def design_cross_flow_lift(self, N_train):
        M_SS_IR_pipe, M_SS_IR_pump = self.design_generic(
            Q_mgd=self.Q_mgd,
            N_pump=2, # one for waste sludge and one for return sludge/clarifier
            L_s=150, # length of suction pipe per filter, [ft]
            L_d_R=30, # pipe length per filter
            H_ts=12, # H_ds_LIFT (12) - H_ss_LIFT (0)
            H_p_R=0. # TMP in water head, [ft]
            )


    def design_cheimcal(self, Q_hr):
        # Hold two weeks of chemicals, assuming cubic in shape
        # V_naocl = self.ins[2].F_vol * 24 * 7 # [m3]
        # V_citric = self.ins[3].F_vol * 24 * 7 # [m3]
        V_CHEM = Q_hr * 24 * 7 # [m3]
        # Volume of HDPE, [m3], 0.003 is the thickness of the container in [m]
        V_HDPE = 0.003 * (V_CHEM**(1/3))**2*6
        # Mass of HDPE, [m3], 950 is the density of the HDPE in [kg/m3]
        M_HDPE = 950 * V_HDPE
        Q_CHEM_mgd = Q_hr*_m3_to_gal/1e6 # Q_hr in m3/hr

        H_ss_CHEM = V_CHEM**(1/3) / _ft_to_m
        # 9'-7" is the water level in membrane trains
        # 18" is the distance from C/L of the pump to the ground
        H_ds_CHEM = 9 + 7/12 - 18/12
        H_ts_CHEM = H_ds_CHEM - H_ss_CHEM

        M_SS_CHEM_pipe, M_SS_CHEM_pump = WWTpump.design_generic(
            Q_mgd=Q_CHEM_mgd,
            N_pump=1,
            L_s=0., # no suction pipe
            L_d=30.,
            H_ts=H_ts_CHEM, # H_ds_IR (5) - H_ss_IR (0)
            H_p=0. # no pressure
            )

        return M_SS_CHEM_pipe, M_SS_CHEM_pump



    @property
    def pump_type(self):
        '''
        [str] The type of the pump that determines the design algorithms to use.
        Use `valid_pump_type` to see acceptable pump types.
        '''
        return self._pump_type
    @pump_type.setter
    def pump_type(self, i):
        if i.lower() not in self.valid_pump_types:
            raise ValueError(f'The given `pump_type` "{i}" is not valid, '
                             'check `valid_pump_types` for acceptable pump types.')
        self._pump_type = i.lower()

    @property
    def valid_pump_types(self):
        '''[tuple] Acceptable pump types.'''
        return self._valid_pump_types

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