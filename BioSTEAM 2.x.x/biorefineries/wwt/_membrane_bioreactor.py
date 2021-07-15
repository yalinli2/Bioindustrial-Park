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
TODO:
    - Generic MBR algorithms for AnMBR and AeMBR
    - Add algorithms for other configurations

References
----------
[1] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
Valorization of Dilute Organic Carbon Waste Streams.
Energy Environ. Sci. 2016, 9 (3), 1102–1112.
https://doi.org/10.1039/C5EE03715H.
'''

import math
import biosteam as bst
import thermosteam as tmo
from biosteam.exceptions import DesignError
#!!! Need to enable relative importing
from _wwt_pump import WWTpump
from utils import (
    auom,
    compute_stream_COD,
    format_str,
    get_BD_dct,
    )

__all__ = ('AnMBR',)

_ft_to_m = auom('ft').conversion_factor('m')
_ft2_to_m2 = auom('ft2').conversion_factor('m2')
_ft3_to_m3 = auom('ft3').conversion_factor('m3')
_m3_to_gal = auom('m3').conversion_factor('gal')
_cmh_to_mgd = _m3_to_gal * 24 / 1e6 # cubic meter per hour to million gallon per day

_d_to_A = lambda d: math.pi/4*(d**2)
_A_to_d = lambda A: ((4*A)/math.pi)**0.5

def _get_d_AF(Q_cmd, R_AF, N_AF, HL_AF, A_AF, V_m_AF):
    # X-sectional area of each filter, [m2]
    A_AF = (Q_cmd/24) * (1+R_AF) / N_AF / HL_AF
    # Diameter of each filter, [m]
    d_AF = (4*A_AF/math.pi) ** 0.5
    # Depth of each filter, [m]
    D_AF = V_m_AF / A_AF
    return d_AF, D_AF


# %%

class AnMBR(bst.Unit):
    '''
    Anaerobic membrane bioreactor (AnMBR) for wastewater treatment as in
    Shoener et al. [1]_

    In addition to the anaerobic treatment, an optional second stage can be added,
    which can be aerobic filter or granular activated carbon (GAC).

    Parameters
    ----------
    reactor_type : str
        Can either be "CSTR" for continuous stirred tank reactor
        or "AF" for anaerobic filter.
    membrane_configuration : str
        Can either be "cross-flow" or "submerged".
    membrane_type : str
        Can be "hollow fiber" ("submerged" configuration only),
        "flat sheet" (either "cross-flow" or "submerged" configuration),
        or "multi-tube" ("cross-flow" configuration only).
    membrane_material : str
        Can be any of the plastics ("PES", "PVDF", "PET", "PTFE")
        for any of the membrane types ("hollow fiber", "flat sheet", "multi-tube"),
        or "sintered steel" for "flat sheet",
        or "ceramic" for "multi-tube".
    include_aerobic_filter : bool
        Whether to include an aerobic filtration process in this AnMBR,
        can only be "True" in "AF" (not "CSTR") reactor.
    add_GAC : bool
        If to add granual activated carbon to enhance biomass retention,
        can only be "True" for the "submerged" configuration.
    include_degassing_membrane : bool
        If to include a degassing membrane to enhance methane
        (generated through the digestion reaction) recovery.

    References
    ----------
    .. [1] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
    Valorization of Dilute Organic Carbon Waste Streams.
    Energy Environ. Sci. 2016, 9 (3), 1102–1112.
    https://doi.org/10.1039/C5EE03715H.
    '''
    _N_ins = 5
    _N_outs = 3

    # Equipment-related parameters
    _N_train = 2
    _cas_per_tank_spare = 2

    _mod_surface_area = {
        'hollow fiber': 370,
        'flat sheet': 1.45/_ft2_to_m2,
        'multi-tube': 32/_ft2_to_m2
        }

    _mod_per_cas = None
    _mod_per_cas_range = {
        'hollow fiber': (30, 48), # min, max
        'flat sheet': (150, 200),
        'multi-tube': (44, 48)
        }

    _cas_per_tank = None
    _cas_per_tank_range = (16, 22)

    _N_blower = 0

    _W_tank = 21
    _D_tank = 12

    _W_dist = 4.5
    _W_eff = 4.5

    _L_well = 8
    _W_well = 8
    _D_well = 8

    _excav_slope = 1.5 # horizontal over vertical
    _excav_access = 3 # construction access, [ft]

    # Operation-related parameters
    _HRT = 10
    _J_max = 8.5
    _TMP_dct = {
        'cross-flow': 21,
        'submerged': 25,
        }
    _TMP_aerobic = None
    _SGD = 1.7
    _v_cross_flow = 3






    #!!! Related to the A_B_script, might not need this
    # # Kinetic and stoichiometric parameters
    # q_hat = 12 # [mg BOD/mg VSS-d]
    # K = 20 # [mg COD/L]
    # Y = 0.5 # biomass yield [mg TSS/mg BOD]
    # b = 0.396 # biomass decay, [1/d]
    # f_d = 0.2 # fraction of biomass that remains as cell debris
    # q_UAP = 1.8 # utilization-associated products, [mg COD/mg VSS-d]
    # q_BAP = 0.1 # biomass-associated products, [mg COD/mg VSS-d]
    # K_UAP = 100 # [mg COD/L]
    # K_BAP = 85 # [mg COD/L]
    # k1 = 0.12 # [mg COD/mg BOD]
    # k2 = 0.09 # [mg COD/mg VSS-d]


    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 reactor_type='CSTR',
                 membrane_configuration='cross-flow',
                 membrane_type='multi-tube',
                 membrane_material='ceramic',
                 include_aerobic_filter=False,
                 add_GAC=False,
                 include_degassing_membrane=False,
                 **kwargs):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.reactor_type = reactor_type
        self.include_aerobic_filter = include_aerobic_filter
        self.membrane_configuration = membrane_configuration
        self.membrane_type = membrane_type
        self.membrane_material = membrane_material
        self.add_GAC = add_GAC
        self.include_degassing_membrane = include_degassing_membrane

        self.Pumps = dict.fromkeys(
            ('Reactor', 'Membrane', 'Degassing'),
            (None, None, None)) # each of the value should be a `WWTpump` object

        for k, v in kwargs.items():
            setattr(self, k, v)

        self._check_design()


    def _check_design(self):
        reactor_type = self.reactor_type
        config = self.membrane_configuration
        m_type = self.membrane_type
        m_material = self.membrane_material

        if reactor_type == 'CSTR':
            if self.aerobic_filter:
                raise DesignError('Aerobic filteration cannot be used in CSTR.')

        if config == 'submerged':
            if not m_type in ('hollow fiber', 'flat sheet'):
                raise DesignError('Only "hollow fiber" or "flat sheet" is allowed '
                                  'for "submerged" membrane, not {m_type}.')
        else: # cross-flow
            if not m_type in ('flat sheet', 'multi-tube'):
                raise DesignError('Only "flat sheet" or "multi-tube" is allowed '
                                  'for "cross-flow" membrane, not {m_type}.')
            if self.add_GAC:
                raise DesignError('No GAC should be added '
                                  '(i.e., `add_GAC` can only be "False"'
                                  'for "cross-flow" membrane.')

        plastics = ('PES', 'PVDF', 'PET', 'PTFE')
        if m_type == 'hollow fiber':
            if not m_material in plastics:
                raise DesignError(f'Only plastic materials {plastics} '
                                  'allowed for "hollow fiber" membrane',
                                  f'not "{m_material}".')
        elif m_type == 'flat sheet':
            if not m_material in (*plastics, 'sintered steel'):
                raise DesignError(f'Only plastic materials {plastics} and "sintered steel"'
                                  'allowed for "flat sheet" membrane',
                                  f'not "{m_material}".')
        else: # multi-tube
            if not m_material in (*plastics, 'ceramic'):
                raise DesignError(f'Only plastic materials {plastics} and "ceramic"'
                                  'allowed for "multi-tube" membrane',
                                  f'not "{m_material}".')



    def _run(self):
        inf, recycled, naocl, citric, air_in = self.ins
        perm, retent, recir, biogas, air_out = self.outs


        ### Compute reactor numbers ###
        self._set_mod_case_tank_N()

        ### Set reactor configuration ###




        ### Clean in place (CIP) ###
        # 2200 gal/yr/mgd of 12.5 wt% solution, 15% by volume
        #!!! The original codes seem to have a bug
        dose_naocl = (2200*self.Q_mgd) / _m3_to_gal * 1e3 / 365 / 24 # kg/hr solution
        naocl.imass['NaOCl'] = dose_naocl * 0.125
        naocl.imass['H2O'] = dose_naocl - naocl.imass['NaOCl']

        # 600 gal/yr/mgd of 100% solution, 13.8 lb/gal
        citric.imass['CitricAcid'] = (600*self.Q_mgd) * \
            (13.8*auom('lb').conversion_factor('kg')) / 365 / 24 # kg/hr solution

        air_out.link(air_in)
        if not self.add_GAC and self.membrane_configuration=='submerged':
            gas_ft3_min = self._compute_gas_demand()
            self._design_blower(gas_ft3_min)

            gas_m3_hr = gas_ft3_min / _ft3_to_m3 * 60 # ft3/min to m3/hr
            air_in.ivol['N2'] = 0.79
            air_in.ivol['O2'] = 0.21
            air_in.F_vol = gas_m3_hr

        else:
            air_in.empty()



    def _set_mod_case_tank_N(self):
        mod_per_cas_range, cas_per_tank_range = \
            self.mod_per_cas_range, self.cas_per_tank_range

        mod_per_cas = self.mod_per_cas = mod_per_cas_range[0]
        cas_per_tank = self.cas_per_tank = cas_per_tank_range[0]

        J, J_max, N_train = self.J, self.J_max, self.N_train
        while J > J_max:
            mod_per_cas += 1
            if mod_per_cas == mod_per_cas_range[1] + 1:
                if cas_per_tank == cas_per_tank_range[1] + 1:
                    N_train += 1
                    mod_per_cas, cas_per_tank = \
                        mod_per_cas_range[0], cas_per_tank_range[0]
                else:
                    cas_per_tank += 1
                    mod_per_cas = mod_per_cas_range[0]


    #!!! No sparging if submerged or using GAC
    # Maybe make blower an another class
    def _compute_gas_demand(self): # for sparging
        gas = self.SGD * self.mod_surface_area*_ft2_to_m2 # [m3/h]
        gas /= (_ft3_to_m3 * 60) # [ft3/min]
        gas_train = gas * self.N_train*self.cas_per_tank*self.mod_per_cas
        return math.ceil(gas_train)


    def _compute_blower_N(self, gas_ft3_min):
        TCFM = gas_ft3_min
        N = 1
        if TCFM <= 30000:
            CFMB = TCFM / N
            while CFMB > 7500:
                N += 1
                CFMB = TCFM / N
        elif 30000 < TCFM <= 72000:
            CFMB = TCFM / N
            while CFMB > 18000:
                N += 1
                CFMB = TCFM / N
        else:
            CFMB = TCFM / N
            while CFMB > 100000:
                N += 1
                CFMB = TCFM / N

        self._N_blower = N + 1 # add a spare





    def _design(self):
        # Reactor and membrane tanks
        self._design_reactor()

        # Membrane
        self._design_membrane()

        # Pumps


    ### Reactor and membrane tanks ###
    def _design_reactor(self):
        # Call the corresponding design function
        # (_design_CSTR or _design_AF)
        func = getattr(self, f'_design_{self.reactor_type}')
        func()


    def _design_CSTR(self):
        N_train = self.N_train
        L_dist, L_CSTR , L_eff, L_membrane_tank = \
            self.L_dist, self.L_CSTR, self.L_eff, self.L_membrane_tank
        W_PB, W_BB, D_tank = self.W_PB, self.W_BB, self.D_tank
        t_wall, t_slab = self.t_wall, self.t_slab
        SL, CA = self.excav_slope, self.constr_access
        D = self.design_results

        ### Concrete calculation ###
        W_N_trains = (self.W_train+2*t_wall)*N_train - t_wall*(N_train-1)

        D = D_tank + 2 # 2add  ft of freeboard
        t = t_wall + t_slab

        get_VWC = lambda L1, N: N * t_wall * L1 * D
        get_VSC = lambda L2: t * L2 * W_N_trains

        # Concrete for distribution channel, [ft3]
        VWC_dist = get_VWC(L1=(W_N_trains+L_dist), N=2)
        VSC_dist = get_VSC(L2=(L_dist+2*t_wall))

        # Concrete for CSTR tanks, [ft3]
        VWC_CSTR = get_VWC(L1=L_CSTR, N=(N_train+1))
        VSC_CSTR = get_VSC(L2=L_CSTR)

        # Concrete for effluent channel, [ft3]
        VWC_eff = get_VWC(L1=(W_N_trains+L_eff), N=2)
        VSC_eff = get_VSC(L2=(L_eff+2*t_wall))

        # Concrete for the pump/blower building, [ft3]
        VWC_PBB = get_VWC(L1=(W_N_trains+W_PB+W_BB), N=2)
        VSC_PBB = get_VSC(L2=(W_PB+t_wall+W_BB))

        if self.membrane_configuration == 'submerged':
            VWC_membrane_tank, VSC_membrane_tank, VWC_well, VSC_well = \
                self._design_membrane_tank(self, D, N_train, W_N_trains,
                                           self.L_membrane_tank)
        else:
            VWC_membrane_tank, VSC_membrane_tank, VWC_well, VSC_well = 0., 0., 0., 0.

        # Total volume of wall concrete, [ft3]
        VWC = VWC_dist + VWC_CSTR + VWC_eff + VWC_PBB + VWC_membrane_tank + VWC_well

        # Total volume of slab concrete [ft3]
        VSC = VSC_dist + VSC_CSTR + VSC_eff + VSC_PBB + VSC_membrane_tank + VSC_well

        D['Tank concrete'] = VWC + VSC

        ### Excavation calculation ###
        get_VEX = lambda L_bttom, W_bottom, diff: \
            0.5 * D_tank * (L_bottom*W_bottom+(L_bottom+diff)*(W_bottom+diff)) # bottom+top

        L_bottom = L_dist + L_CSTR + L_eff + L_membrane_tank + 2*CA
        W_bottom = W_N_trains + 2*CA
        diff = D_tank * SL

        # Excavation volume for membrane tanks, [ft3]
        VEX_membrane_tank = get_VEX(L_bottom, W_bottom, diff)

        # Excavation volume for pump/blower building, [ft3]
        VEX_PBB = get_VEX((W_PB+W_BB+2*CA), W_bottom, diff)

        D['Excavation'] = VEX_membrane_tank + VEX_PBB


    def _design_AF(self):
        '''NOT READY YET.'''
    #!!! Update/recycle this for AF
    def _design_anaerobic_filter(
            self, Q_mgd,
            Ss, # readily biodegradable (soluble) substrate concentration, [kg COD/m3]
            Sp, # slowly biodegradable (particulate) substrate concentration, [kg COD/m3]
            OLR_AF, # organic loading rate, [kg-COD/m3/day]
            HL_AF, # hydraulic loading rate, [m3/m2/hr]
            R_AF # recirculation ratio
            ):

        ### Filter material ###
        N_AF = 2
        Q_cmd = Q_mgd *1e6/_m3_to_gal # [m3/day]
        # Volume of packing media in each filter, [m3]
        V_m_AF = (Q_cmd/N_AF) * (Ss+Sp) / OLR_AF
        # Diameter (d) / depth (D) of each filter, [m]
        d_AF, D_AF = _get_d_AF(Q_cmd, R_AF, N_AF, HL_AF, V_m_AF)

        while D_AF > 6: # assumed maximum depth assumption, [m]
            R_AF = R_AF + 0.1;
            d_AF, D_AF = _get_d_AF(Q_cmd, R_AF, N_AF, HL_AF, V_m_AF)

            while d_AF > 12: # assumed maximum diameter, [m]
                N_AF = N_AF + 1;
                d_AF, D_AF = _get_d_AF(Q_cmd, R_AF, N_AF, HL_AF)

        # Unit conversion
        d_AF /= _ft_to_m # [ft]
        D_AF /= _ft_to_m # [ft]
        V_m_AF /= _ft3_to_m3 # [ft3]

        ### Concrete material ###
        # External wall concrete, [ft3]
        # 6/12 is wall thickness and 3 is freeboard
        VWC_AF = N_AF * 6/12 * math.pi * d_AF * (D_AF+3)
        VWC_AF *= N_AF
        # Floor slab concrete, [ft3]
        # 8/12 is slab thickness
        VSC_AF = _d_to_A(d_AF)+ 8/12 * _d_to_A(d_AF)
        VSC_AF *= N_AF

        ### Excavation ###
        SL = 1.5 # slope = horizontal/vertical
        CA = 3 # construction Access, [ft]
        #  Excavation of pump building
        PBL, PBW, PBD = 50, 30, 10 # pump building length, width, depth, [ft]
        Area_B_P = (PBL+2*CA) * (PBW+2*CA) # bottom area of frustum, [ft2]
        Area_T_P = (PBL+2*CA+PBW*SL) * (PBW+2*CA+PBD*SL) # top area of frustum, [ft2]
        VEX_PB = 0.5 * (Area_B_P+Area_T_P) * PBD # total volume of excavaion, [ft3]

        return N_AF, d_AF, D_AF, V_m_AF, VWC_AF, VWC_AF, VEX_PB



    def _design_packing_media(self, V):
        # Assume 50%/50% wt/wt LDPE/HDPE
        # 0.9 is void fraction, usually 85% - 95% for plastic packing media
        # 925 is density of LDPE (910-940), [kg/m3]
        # 950 is density of LDPE (930-970), [kg/m3]
        # M_LDPE_kg = 0.5 * (1-0.9) * 925 * V_m
        # M_HDPE_kg = 0.5 * (1-0.9) * 950 * V_m
        return 46.25*V, 47.5*V



    # Called by _design_CSTR/_design_AF
    def _design_membrane_tank(self, D, N_train, W_N_trains, L_membrane_tank,
                              t_wall, t_slab):
        L_well, W_well, D_well = self.L_well, self.W_well, self.D_well

        # Concrete for membrane tanks, [ft3]
        t = t_wall + t_slab
        VWC_membrane_tank = (N_train+1) * t_wall * L_membrane_tank * D
        VSC_membrane_tank = t * L_membrane_tank * W_N_trains

        # Concrete for wet well (mixed liquor storage), [ft3]
        L = L_well + 2*t_wall
        W = W_well + 2*t_wall

        VWC_well = 2 * t_wall * (L_well+W) * D_well
        VSC_well = (t_slab+t_wall) * L * W

        return VWC_membrane_tank, VSC_membrane_tank, VWC_well, VSC_well


    ### Membrane ###
    def _design_membrane(self):
        # Call the corresponding design function
        # (_design_hollow_fiber, _design_flat_sheet, or _design_multi_tube)
        m_type = format_str(self.membrane_type)
        func = getattr(self, f'_design_{m_type}')
        func()


    def _design_hollow_fiber(self):
        '''NOT READY YET.'''


    def _design_flat_sheet(self):
        '''NOT READY YET.'''


    def _design_multi_tube(self,):
        D = self.design_results
        # Cross-flow flow rate per module,
        # based on manufacture specifications for compact 33, [m3/hr]
        Q_cross_flow = 53.5 * self.v_cross_flow

        N_tot = self.N_train * self.cas_per_tank * self.mod_per_cas # total number of modules

        Q_R_cmh = N_tot * Q_cross_flow # total retentate flow rate, [m3/hr]
        Q_R_mgd = Q_R_cmh * _cmh_to_mgd # [mgd]

        M_tot = N_tot * 263.05
        # # 263.05 is volume of material for each membrane tube, [m3]
        # #      L_tube               OD        ID
        # V_tube = 3 * math.pi/4 * ((6e-3)**2-(5.2e-3)**2)
        # V_SU = 700 * V_tube # V for each small unit [m3]
        # M_SU = 1.78*10e3 * V_SU # mass = density*volume, [kg/m3]

        D['Membrane'] = M_tot
        D['Retentate flow'] = Q_R_mgd


    ### Pumps ###
    def _design_pump(self):
        pumps = (None, None, None, None)

        # Permeate
        inputs = (self.cas_per_tank, self.D_tank, self.TMP_anaerobic,
                  self.include_aerobic_filter)
        pumps[0] = WWTpump(ID=f'{self.ID}_PERM', ins=self.outs[0].proxy(),
                           pump_type='cross-flow_permeate', inputs=inputs)

        PAUSED, ADD OTHER PUMPS

        self._pump_dct = dict.fromkeys(
            ('permeate', 'retentate', 'recirculation', 'lift'), pumps)

        pipe_ss, pump_ss = 0., 0.
        for p in pumps:
            if p is not None:
                p.simulate()
                pipe_ss += p.design_results['Pipe stainless steel']
                pump_ss += p.design_results['Pump stainless steel']

        D = self.design_results
        D['Pipe stainless steel'] = pipe_ss
        D['Pump stainless steel'] = pump_ss



    def _cost(self):
        pump_power = sum(v.power_utility.rate for v in self.pump_dct.values())

        self.power_utility.rate = pump_power #!!! also sparging, etc.







    ### Reactor configuration ###
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
    def membrane_type(self):
        '''
        [str] Can be "hollow fiber" ("submerged" configuration only),
        "flat sheet" (either "cross-flow" or "submerged" configuration),
        or "multi-tube" ("cross-flow" configuration only).
        '''
        return self._membrane_type
    @membrane_type.setter
    def membrane_type(self, i):
        i = 'multi-tube' if i.lower() in ('multi tube', 'multitube') else i
        if not i.lower() in ('hollow fiber', 'flat sheet', 'multi-tube'):
            raise ValueError('`membrane_type` can only be "hollow fiber", '
                             f'"flat sheet", or "multi-tube", not "{i}".')
        self._membrane_type = i.lower()

    @property
    def membrane_material(self):
        '''
        [str] Can be any of the plastics ("PES", "PVDF", "PET", "PTFE")
        for any of the membrane types ("hollow fiber", "flat sheet", "multi-tube"),
        or "sintered steel" for "flat sheet",
        or "ceramic" for "multi-tube".
        '''
        return self._membrane_material
    @membrane_material.setter
    def membrane_material(self, i):
        plastics = ('PES', 'PVDF', 'PET', 'PTFE')
        if i.upper() in plastics:
            self._membrane_material = i.upper()
        elif i.lower() in ('sintered steel', 'ceramic'):
            self._membrane_material = i.lower()
        else:
            raise ValueError(f'`membrane_material` can only be plastics materials '
                             f'{plastics}, "sintered steel", or "ceramic", not {i}.')


    ### Reactor/membrane tank ###
    @property
    def N_train(self):
        '''[int] Number of treatment train.'''
        return self._N_train
    @N_train.setter
    def N_train(self, i):
        self._N_train = math.ceil(i)

    @property
    def cas_per_tank_spare(self):
        '''[int] Number of spare cassettes per train.'''
        return self._cas_per_tank_spare
    @cas_per_tank_spare.setter
    def cas_per_tank_spare(self, i):
        self._cas_per_tank_spare = i

    @property
    def mod_per_cas_range(self):
        '''
        [tuple] Range (min, max) of the number of membrane modules per cassette
        for the current membrane type.
        '''
        return self._mod_per_cas_range
    @mod_per_cas_range.setter
    def mod_per_cas_range(self, i):
        self._mod_per_cas_range[self.membrane_type] = tuple(i)

    @property
    def mod_per_cas(self):
        '''
        [float] Number of membrane modules per cassette for the current membrane type.
        '''
        return self._mod_per_cas or self._mod_per_cas_range[self.membrane_type][0]

    @property
    def cas_per_tank_range(self):
        '''
        [tuple] Range (min, max) of the number of membrane cassette per tank
        (same for all membrane types).
        '''
        return self._cas_per_tank_range
    @cas_per_tank_range.setter
    def cas_per_tank_range(self, i):
        self._cas_per_tank_range = tuple(i)

    @property
    def cas_per_tank(self):
        '''
        [float] Number of membrane cassettes per tank for the current membrane type.
        '''
        return self._cas_per_tank or self._cas_per_tank_range [0]

    @property
    def mod_surface_area(self):
        '''
        [float] Surface area of the membrane for the current membrane type, [m2/module].
        Note that one module is one sheet for plat sheet and one tube for multi-tube.
        '''
        return self._mod_surface_area[self.membrane_type]
    @mod_surface_area.setter
    def mod_surface_area(self, i):
        self._mod_surface_area[self.membrane_type] = float(i)

    @property
    def L_CSTR(self):
        '''[float] Length of the CSTR tank, [ft].'''
        if self.reactor_type == 'AF':
            return 0
        return self.ins[0].F_vol/_ft3_to_m3*self.HRT/(self.N_train*self.W_tank*self.D_tank)

    @property
    def L_membrane_tank(self):
        '''[float] Length of the membrane tank, [ft].'''
        return math.ceil((self.cas_per_tank+self.cas_per_tank_spare)*3.4)

    @property
    def W_tank(self):
        '''[float] Width of the reactor/membrane tank (same value), [ft].'''
        return self._W_tank
    @W_tank.setter
    def W_tank(self, i):
        self._W_tank = float(i)

    @property
    def D_tank(self):
        '''[float] Depth of the reactor/membrane tank (same value), [ft].'''
        return self._D_tank
    @D_tank.setter
    def D_tank(self, i):
        self._D_tank = float(i)

    @property
    def W_dist(self):
        '''[float] Width of the distribution channel, [ft].'''
        return self._W_dist
    @W_dist.setter
    def W_dist(self, i):
        self._W_dist = float(i)

    @property
    def W_eff(self):
        '''[float] Width of the effluent channel, [ft].'''
        return self._W_eff
    @W_eff.setter
    def W_eff(self, i):
        self._W_eff = float(i)

    @property
    def t_wall(self):
        '''
        [float] Concrete wall thickness, [ft].
        Minimum of 1 ft with 1 in added for every ft of depth over 12 ft.
        '''
        return 1 + max(self.D_tank-12, 0)/12

    @property
    def t_slab(self):
        '''
        [float] Concrete slab thickness, [ft].
        2 in thicker than the wall thickness.
        '''
        return self.t_wall+2/12


    ### Pump/blower ###
    @property
    def pump_dct(self):
        '''
        [dict] All pumps included in this unit, will be automatically updated
        during simulation.
        Keys are "permeate", "retentate", "recirculation", and "lift",
        values are :class:`WWTpump` objects (or None if not applicable).
        '''
        return self._pump_dct


    @property
    def N_blower(self):
        '''
        [int] Number of blowers needed for gas sparging
        (not needed for some designs).
        '''
        if not self.add_GAC and self.membrane_configuration=='submerged':
            return self._N_blower
        return 0

    @property
    def W_PB(self):
        '''[float] Width of the pump building, [ft].'''
        if self.membrane_configuration == 'submerged':
            N = self.cas_per_tank
        else: # cross-flow
            N = math.ceil(self.L_CSTR/((1+8/12)+(3+4/12)))

        if 0 <= N <= 10:
            W_PB = 27 + 4/12
        elif 11 <= N <= 16:
            W_PB = 29 + 6/12
        elif 17 <= N <= 22:
            W_PB = 31 + 8/12
        elif 23 <= N <= 28:
            W_PB = 35
        elif N >= 29:
            W_PB = 38 + 4/12
        else:
            W_PB = 0

        return W_PB

    @property
    def L_BB(self):
        '''[float] Length of the blower building, [ft].'''
        if self.membrane_configuration == 'submerged':
            return (69+6/12) if self.cas_per_tank<=18 else (76 + 8/12)
        return 0

    @property
    def W_BB(self):
        '''[float] Width of the blower building, [ft].'''
        if self.membrane_configuration == 'submerged':
            return (18+8/12) if self.cas_per_tank<=18 else 22
        return 0


    ### Wet well (submerged only) ###
    @property
    def L_well(self):
        '''
        [float] Length of the wet well, [ft].
        Only needed for submerged configuration.
        '''
        return self._L_well if self.membrane_configuration == 'submerged' else 0
    @L_well.setter
    def L_well(self, i):
        self._L_well = float(i)

    @property
    def W_well(self):
        '''
        [float] Width of the wet well, [ft].
        Only needed for submerged configuration.
        '''
        return self._W_well if self.membrane_configuration == 'submerged' else 0
    @W_well.setter
    def W_well(self, i):
        self._W_well = float(i)

    @property
    def D_well(self):
        '''
        [float] Depth of the wet well, [ft].
        Only needed for submerged configuration.
        '''
        return self._D_well if self.membrane_configuration == 'submerged' else 0
    @D_well.setter
    def D_well(self, i):
        self._D_well = float(i)


    ### Excavation ###
    @property
    def excav_slope(self):
        '''[float] Slope for excavation (horizontal/vertical).'''
        return self._excav_slope
    @excav_slope.setter
    def excav_slope(self, i):
        self._excav_slope = float(i)

    @property
    def constr_access(self):
        '''[float] Extra room for construction access, [ft].'''
        return self._constr_access
    @constr_access.setter
    def constr_access(self, i):
        self._constr_access = float(i)



    ### Operation-related parameters ###
    @property
    def Q_mgd(self):
        '''
        [float] Volumetric flow rate in million gallon per day, [mgd].
        '''
        return self.ins[0].F_vol*_m3_to_gal*24/1e6

    @property
    def HRT(self):
        '''
        [float] Hydraulic retention time, [hr].
        '''
        return self._HRT
    @HRT.setter
    def HRT(self, i):
        self._HRT = float(i)

    @property
    def J_max(self):
        '''[float] Maximum membrane flux, [L/m2/hr].'''
        return self._J_max
    @J_max.setter
    def J_max(self, i):
        self._J_max = float(i)

    @property
    def J(self):
        '''[float] Membrane flux, [L/m2/hr].'''
        # Based on the flux of one train being offline
        SA = (self.N_train-1) * self.cas_per_tank * self.mod_per_cas * self.mod_surface_area
        return self.ins[0].F_vol*1e3/SA # 1e3 is conversion from m3 to L

    @property
    def TMP_anaerobic(self):
        '''[float] Transmembrane pressure in the anaerobic reactor, [psi].'''
        return self._TMP_dct[self.membrane_configuration]
    @TMP_anaerobic.setter
    def TMP_anaerobic(self, i):
        self._TMP_dct[self.membrane_configuration] = float(i)

    @property
    def TMP_aerobic(self):
        '''
        [float] Transmembrane pressure in the aerobic filter, [psi].
        Defaulted to half of the reactor TMP.
        '''
        if not self._include_aerobic_filter:
            return 0.
        else:
            return self._TMP_aerobic or self._TMP_dct[self.membrane_configuration]/2
    @TMP_aerobic.setter
    def TMP_aerobic(self, i):
        self._TMP_aerobic = float(i)

    @property
    def SGD(self):
        '''[float] Specific gas demand, [m3 gas/m2 membrane area/h].'''
        return self._SGD
    @SGD.setter
    def SGD(self, i):
        self._SGD = float(i)

    @property
    def v_cross_flow(self):
        '''
        [float] Cross-flow velocity, [m/s].
        '''
        return self._v_cross_flow if self.membrane_configuration=='cross-flow' else 0
    @v_cross_flow.setter
    def v_cross_flow(self, i):
        self._v_cross_flow = float(i)