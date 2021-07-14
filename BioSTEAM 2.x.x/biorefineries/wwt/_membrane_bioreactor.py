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
import sympy as sp
import biosteam as bst
import thermosteam as tmo
from biosteam.exceptions import DesignError
#!!! Need to enable relative importing
from utils import (
    auom,
    compute_stream_COD,
    format_str,
    get_BD_dct,
    )

__all__ = ('AnMBR',)

_ft_to_m = auom('ft').conversion_factor('m')
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

    Parameters
    ----------
    reactor_type : str
        Can either be "CSTR" for continuous stirred tank reactor
        or "AF" for anaerobic filter.
    aerobic_filteration : bool
        If to include an aerobic filtration process in this AnMBR,
        can only be "True" in "AF" (not "CSTR") reactor.
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
    _N_ins = 4
    _N_outs = 3

    HRT = 10 # hydraulic retention time, [hr]
    N_train = 2 #!!! number of trains, need to be updated in the algorithm

    W_train = 21 # width of one train, [ft]
    D_train = 12 # depth of one train, [ft]
    L_dist = 4.5 # width of distribution channel, [ft]
    L_eff = 4.5 # width of effluent channel, [ft]

    excav_slope = 1.5 # horizontal over vertical
    excav_access = 3 # construction access, [ft]


    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 reactor_type='CSTR',
                 membrane_configuration='cross-flow',
                 membrane_type='multi-tube',
                 membrane_material='ceramic',
                 include_aerobic_filteration=False,
                 add_GAC=False,
                 include_degassing_membrane=False,
                 **kwargs):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.reactor_type = reactor_type
        self.include_aerobic_filteration = include_aerobic_filteration
        self.membrane_configuration = membrane_configuration
        self.membrane_type = membrane_type
        self.membrane_material = membrane_material
        self.add_GAC = add_GAC
        self.include_degassing_membrane = include_degassing_membrane

        for k, v in kwargs.items():
            setattr(self, k, v)

        self._check_design()


    def _check_design(self):
        rx_type = self.reactor_type
        config = self.membrane_configuration
        m_type = self.membrane_type
        m_material = self.membrane_material

        if rx_type == 'CSTR':
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
        inf, recycled, naocl, citric = self.ins


        ### Clean in place (CIP) ###
        # 2200 gal/yr/mgd of 12.5 wt% solution, 15% by volume
        #!!! The original codes seem to have a bug
        dose_naocl = (2200*self.Q_mgd) / _m3_to_gal * 1e3 / 365 / 24 # kg/hr solution
        naocl.imass['NaOCl'] = dose_naocl * 0.125
        naocl.imass['H2O'] = dose_naocl - naocl.imass['NaOCl']

        # 600 gal/yr/mgd of 100% solution, 13.8 lb/gal
        citric.imass['CitricAcid'] = (600*self.Q_mgd) * \
            (13.8*auom('lb').conversion_factor('kg')) / 365 / 24 # kg/hr solution



    #!!! No sparging if submerged and using GAC

    def _design(self):
        N_train, D_train, L_dist, L_CSTR , L_eff, W_PB, W_BB, SL, CA = \
            self.N_train, self.D_train, self.L_dist, self.L_CSTR, self.L_eff, \
            self.W_PB, self.L_membrane_tank, self.excav_slope, self.excav_access

        D = self.design_results

        ### Concrete calculation ###
        t_wall = 1 + max(self.D_train-12, 0)/12 # concrete wall thickness, [ft]
        t_slab = t_wall + 2/12 # concrete slab thickness, [ft]
        W_N_trains = (self.W_train+2*t_wall)*N_train - t_wall*(N_train-1)

        # Concrete for distribution channel, [ft3]
        VWC_dist = (D_train+2) * t_wall * (2*W_N_trains+2*L_dist)
        VSC_dist = W_N_trains * (L_dist+2*t_wall)*t_wall + \
            W_N_trains*(L_dist+2*t_wall)*t_slab

        # Concrete for CSTR trains, [ft3]
        VWC_CSTR = (D_train+2) * t_wall * (N_train+1) * L_CSTR;
        VSC_CSTR = W_N_trains*L_CSTR*t_wall + W_N_trains*t_slab*L_CSTR

        # Concrete for effluent channel, [ft3]
        VWC_eff = (D_train+2) * t_wall * (2*W_N_trains+2*L_eff)
        VSC_eff = W_N_trains * (L_eff+2*t_wall)*t_wall+ \
            W_N_trains*(L_eff+2*t_wall)*t_slab;

        # Concrete for the pump/blower building, [ft3]
        VWC_PBB = (D_train+2) * t_wall * (2*W_N_trains+2* W_PB+2*W_BB)
        VSC_PBB = W_N_trains * (W_PB+t_wall+W_BB)*t_wall + \
            W_N_trains*(W_PB+t_wall+W_BB)*t_slab

        # Total volume of wall concrete, [ft3]
        VWC = VWC_dist + VWC_CSTR + VWC_eff + VWC_PBB

        # Total volume of slab concrete [ft3]
        VSC = VSC_dist + VSC_CSTR + VSC_eff + VSC_PBB

        D['Total concrete'] = VWC + VSC

        ### Excavation calculation ###
        # Excavation volume for membrane trains, [ft3]
        Area_B_train = (L_dist+L_CSTR+L_eff+2*CA) * (W_N_trains+2*CA) # bottom
        Area_T_train = (L_dist+L_CSTR+L_eff+2*CA+D_train*SL) * \
            (W_N_trains+2*CA+D_train*SL) # top
        VEX_train = 0.5 * (Area_B_train+Area_T_train) * D_train

        # Excavation volume for membrane trains, [ft3]
        Area_B_PBB = (W_PB+W_BB+2*CA) * (W_N_trains+2*CA);
        Area_T_PBB = (W_PB+W_BB+2*CA+D_train*SL) * \
            (W_N_trains+2*CA+D_train*SL)
        VEX_PB = 0.5 * (Area_B_PBB+Area_T_PBB) * D_train

        D['Total excavation'] = VEX_train + VEX_PB



    def _design_membrane(self, Q_mgd):
        m_type = format_str(self.membrane_type)
        func = getattr(self, f'_design_{m_type}_membrane')
        m_tot, Q_R_mgd = func()



    def _design_flat_sheet_membrane(self):
        pass


    def _design_hollow_fiber_membrane(self):
        pass


    def _design_multi_tube_membrane(
            self,Q_mgd,
            T_bw, # backwashing time per day, [hr]
            v_xflow, # cross-flow velocity, [m/s]
            N_train,
            cas_per_tank, # number of membrane cassette per tank
            mod_per_cas, # number of membrane module per cassette
            module_SA
            ):
        # Cross-flow flow rate per module,
        # based on manuf. specs. for compact 33, [m3/hr]
        Q_xflow = 53.5 * v_xflow

        N_tot = mod_per_cas * cas_per_tank * N_train # total number of modules

        # Total retentate flow rate, [m3/hr]
        Q_R_cmh = N_tot * Q_xflow
        Q_R_mgd = Q_R_cmh * _cmh_to_mgd # [mgd]

        M_tot = N_tot * 263.05
        # # 263.05 is volume of material for each membrane tube, [m3]
        # #      L_tube               OD        ID
        # V_tube = 3 * math.pi/4 * ((6e-3)**2-(5.2e-3)**2)
        # V_SU = 700 * V_tube # V for each small unit [m3]
        # M_SU = 1.78*10e3 * V_SU # mass = density*volume, [kg/m3]

        return Q_R_mgd, M_tot


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
        Area_B_P = (PBL+2*CA) * (PBW+2*CA); # bottom area of frustum, [ft2]
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

    @property
    def L_CSTR(self):
        '''[float] Length of the CSTR tank, [ft].'''
        return self.ins[0].F_vol/_ft3_to_m3*self.HRT/(self.N_train*self.W_train*self.D_train)

    @property
    def N_pump(self):
        '''[int] Number of required pump.'''
        return math.ceil(self.L_CSTR/((1+8/12)+(3+4/12)))

    @property
    def W_PB(self):
        '''[float] Width of pump building, [ft].'''
        N_pump = self.N_pump
        if 0<=N_pump<=10:
            W_PB = 27 + 4/12
        elif 11<=N_pump<=16:
            W_PB = 29 + 6/12
        elif 17<=N_pump<=22:
            W_PB = 31 + 8/12
        elif 23<=N_pump<=28:
            W_PB = 35
        elif N_pump>=29:
            W_PB = 38 + 4/12
        else:
            W_PB = 0.

        return W_PB

    # @property
    # def Q_mgd(self):
    #     '''
    #     [float] Volumetric flow rate in million gallon per day, [mgd].
    #     Will return total volumetric flow through the unit if not provided.
    #     '''
    #     if self._Q_mgd:
    #         return self._Q_mgd
    #     return self.ins[0].F_vol*_m3_to_gal*24/1e6
    # @Q_mgd.setter
    # def Q_mgd(self, i):
    #     self._Q_mgd = i