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
import thermosteam as tmo

# from biorefineries.wwt import WWTpump
# from biorefineries.wwt.utils import (
#     auom,
#     compute_stream_COD,
#     get_digestion_rxns,
#     get_split_dct
#     )
from _wwt_pump import WWTpump
from utils import (
    auom,
    compute_stream_COD,
    get_digestion_rxns,
    get_split_dct
    )

__all__ = ('FilterTank',)

_ft_to_m = auom('ft').conversion_factor('m')
_ft2_to_m2 = auom('ft2').conversion_factor('m2')
_ft3_to_m3 = auom('ft3').conversion_factor('m3')
_ft3_to_gal = auom('ft3').conversion_factor('gallon')
_m3_to_gal = auom('m3').conversion_factor('gal')
_cmh_to_mgd = _m3_to_gal * 24 / 1e6 # cubic meter per hour to million gallon per day
_lb_to_kg = auom('lb').conversion_factor('kg')


_d_to_A = lambda d: math.pi/4*(d**2)
_A_to_d = lambda A: ((4*A)/math.pi)**0.5


# %%

class FilterTank(bst.Unit):
    '''
    A superclass for anaerobic and aerobic filter treatment as in
    Shoener et al. [1]_ Some assumptions adopted from Humbird et al. [2]_

    Parameters
    ----------
    filter_type : str
        Can either be "anaerobic" or "aerobic".
    OLR : float
        Organic loading rate of influent, [kg COD/m3/hr].
    HLR : float
        Hydraulic loading rate of influent, [m3/m2/hr].
    X_decomp : float
        Fraction of the influent COD converted to biogas (`filter_type`=="anaerobic")
        or CO2 (`filter_type`=="aerobic").
    X_growth : float
        Fraction of the influent COD converted to biomass growth.
    split : dict
        Component-wise split to the treated water.
        E.g., {'Water':1, 'WWTsludge':0} indicates all of the water goes to
        the treated water and all of the WWTsludge goes to the wasted sludge.
        Default splits (based on the membrane bioreactor in [2]_) will be used
        if not provided.
    T : float
        Temperature of the filter tank.

    References
    ----------
    .. [1] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
    Valorization of Dilute Organic Carbon Waste Streams.
    Energy Environ. Sci. 2016, 9 (3), 1102–1112.
    https://doi.org/10.1039/C5EE03715H.

    .. [2] Humbird et al., Process Design and Economics for Biochemical Conversion of
    Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic
    Hydrolysis of Corn Stover; Technical Report NREL/TP-5100-47764;
    National Renewable Energy Lab (NREL), 2011.
    https://www.nrel.gov/docs/fy11osti/47764.pdf
    '''
    _N_ins = 3 # influent, recycle, air (optional)
    _N_outs = 4 # biogas (optional), effluent, waste sludge, air (optional)

    _N_filter_min = 2
    _d_max = 12
    _excav_slope = 1.5
    _constr_access = 3


    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 filter_type='aerobic',
                 OLR=1.67, # 2/24/(1-0.95) based on the 2 kg COD/m3/day of effluent in ref [1]
                 HLR=2, X_decomp=0.74, X_growth=0.22, # X_decomp & X_growth from ref[2]
                 split=None, T=30+273.15):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.filter_type = filter_type
        self.OLR = OLR
        self.HLR = HLR
        self.X_decomp = X_decomp
        self.X_growth = X_growth
        self.split = split if split else get_split_dct(self.chemicals)
        self.T = T
        self._refresh_rxns()


    def _refresh_rxns(self, X_decomp=None, X_growth=None):
        X_decomp = X_decomp if X_decomp else self.X_decomp
        X_growth = X_growth if X_growth else self.X_growth

        self._decomp_rxns = get_digestion_rxns(self.ins[0], 1.,
                                               X_decomp, 0., 'WWTsludge')
        self._growth_rxns = get_digestion_rxns(self.ins[0], 1.,
                                               0., X_growth, 'WWTsludge')
        self._i_rm = self._decomp_rxns.X_net.data + self._growth_rxns.X_net.data


    def _run(self):
        raw, recycled, air_in = self.ins
        biogas, eff, waste, air_out = self.outs

        inf = raw.copy()
        inf.mix_from((raw, recycled))
        self._inf = inf.copy() # this stream will be preserved (i.e., no reaction)

        self.growth_rxns(inf.mol)
        self.decomp_rxns.force_reaction(inf.mol)
        tmo.separations.split(inf, eff, waste, self._isplit.data)

        biogas.phase = air_in.phase = air_out.phase = 'g'

        if inf.imol['O2'] < 0:
            air_in.imol['O2'] = - inf.imol['O2']
            air_in.imol['N2'] = - 0.79/0.21 * inf.imol['O2']
            inf.imol['O2'] = 0
        else:
            air_in.empty()

        if self.filter_type == 'anaerobic':
            biogas.receive_vent(eff)
            biogas.receive_vent(waste)
            air_in.empty()
            air_out.empty()
        else:
            biogas.empty()
            air_out.receive_vent(eff)
            air_out.receive_vent(waste)
            air_out.imol['N2'] += air_in.imol['N2']
            self._recir_ratio = None


    def _design(self):
        D = self.design_results
        func = self._design_anaerobic if self.filter_type=='anaerobic' \
            else self._design_aerobic

        ### Concrete and excavation ###
        V, VWC, VSC, VEX = func()
        D['Volume [ft3]'] = V
        D['Wall concrete [ft3]'] = VWC
        D['Slab concrete [ft3]'] = VSC
        D['Excavation [ft3]'] = VEX

        ### Pump ###
        add_inputs = (self.N_filter, self.D)
        Pump =  WWTpump(ID=f'{self.ID}_lift',
                        ins=self.ins[0].proxy(),
                        pump_type='lift',
                        add_inputs=add_inputs)

        Pump.simulate()
        D['Pipe stainless steel [kg]'] = Pump.design_results['Pipe stainless steel [kg]']
        D['Pump stainless steel [kg]'] = Pump.design_results['Pump stainless steel [kg]']
        self._Pump = Pump

        ### Packing ###
        # Assume 50%/50% vol/vol LDPE/HDPE
        # 0.9 is void fraction, usually 85% - 95% for plastic packing media
        # 925 is density of LDPE (910-940), [kg/m3] (not used)
        # 950 is density of LDPE (930-970), [kg/m3] (not used)
        # M_LDPE_kg = 0.5 * (1-0.9) * 925 * V_m
        # M_HDPE_kg = 0.5 * (1-0.9) * 950 * V_m
        D['Packing LDPE [m3]'] = D['Packing HDPE [m3]'] = 0.05 * V


    def _design_aerobic(self):
        inf, N, SL, CA \
            = self._inf, self._N_filter_min, self.excav_slope, self.constr_access
        Q = inf.F_vol

        ### Concrete ###
        get_V = lambda N: ((Q/N)*self.compute_COD(inf)) / self.OLR # m3
        get_A = lambda N: Q/N/self.HLR

        V, A = get_V(N), get_A(N)
        d = _A_to_d(A)
        D = V / d # D is depth

        # Check if more than one filter is needed
        while d > self.d_max:
            N += 1
            V, A = get_V(N), get_A(N)
            d = _A_to_d(A)
            D = V / d

        V_ft3 = V / _ft3_to_m3
        d_ft = d / _ft_to_m
        D_ft = D / _ft_to_m
        self._N_filter, self._d, self._D = N, d, D_ft

        # Volume of wall/slab concrete, [ft3]
        # 8/12 is wall thickness, 3 is freeboard
        VWC = 8/12 * math.pi * d_ft * (D_ft+3)
        VWC *= N

        # 1 is slab thickness
        VSC = 2 * 1 * _d_to_A(d_ft)
        VSC *= N

        ### Excavation ###
        # 50/30/10 are building length/width/depth, [ft]
        L_B, W_B, diff = (50+2*CA), (30+2*CA), (10*SL)
        Area_B = L_B * W_B
        Area_T = (L_B+diff) * (W_B+diff)
        VEX = 0.5 * (Area_B+Area_T) * 10 # [ft3]

        return V_ft3, VWC, VSC, VEX


    # def _design_anaerobic_filter(
    #         self, Q_mgd,
    #         Ss, # readily biodegradable (soluble) substrate concentration, [kg COD/m3]
    #         Sp, # slowly biodegradable (particulate) substrate concentration, [kg COD/m3]
    #         OLR_AF, # organic loading rate, [kg-COD/m3/day]
    #         HL_AF, # hydraulic loading rate, [m3/m2/hr]
    #         R_AF # recirculation ratio
    #         ):

    #     ### Filter material ###
    #     N_AF = 2
    #     Q_cmd = self.Q_cmd
    #     # Volume of the filter packing media in each filter, [m3]
    #     V_m_AF = (Q_cmd/N_AF) * (Ss+Sp) / OLR_AF
    #     # Diameter (d) / depth (D) of each filter, [m]
    #     d_AF, D_AF = _get_d_AF(Q_cmd, R_AF, N_AF, HL_AF, V_m_AF)

    #     while D_AF > 6: # assumed maximum depth assumption, [m]
    #         R_AF = R_AF + 0.1;
    #         d_AF, D_AF = _get_d_AF(Q_cmd, R_AF, N_AF, HL_AF, V_m_AF)

    #         while d_AF > 12: # assumed maximum diameter, [m]
    #             N_AF = N_AF + 1;
    #             d_AF, D_AF = _get_d_AF(Q_cmd, R_AF, N_AF, HL_AF)

    #     # Unit conversion
    #     d_AF /= _ft_to_m # [ft]
    #     D_AF /= _ft_to_m # [ft]
    #     V_m_AF /= _ft3_to_m3 # [ft3]

    #     ### Concrete material ###
    #     # External wall concrete, [ft3]
    #     # 6/12 is wall thickness and 3 is freeboard
    #     VWC_AF = N_AF * 6/12 * math.pi * d_AF * (D_AF+3)
    #     VWC_AF *= N_AF
    #     # Floor slab concrete, [ft3]
    #     # 8/12 is slab thickness
    #     VSC_AF = _d_to_A(d_AF)+ 8/12 * _d_to_A(d_AF)
    #     VSC_AF *= N_AF

    #     ### Excavation ###
    #     SL = 1.5 # slope = horizontal/vertical
    #     CA = 3 # construction Access, [ft]
    #     #  Excavation of pump building
    #     PBL, PBW, PBD = 50, 30, 10 # pump building length, width, depth, [ft]
    #     Area_B_P = (PBL+2*CA) * (PBW+2*CA) # bottom area of frustum, [ft2]
    #     Area_T_P = (PBL+2*CA+PBW*SL) * (PBW+2*CA+PBD*SL) # top area of frustum, [ft2]
    #     VEX_PB = 0.5 * (Area_B_P+Area_T_P) * PBD # total volume of excavaion, [ft3]

    #     return N_AF, d_AF, D_AF, V_m_AF, VWC_AF, VWC_AF, VEX_PB

    def _cost(self):
        D, C, F_BM, lifetime = self.design_results, self.baseline_purchase_costs, \
            self.F_BM, self._default_equipment_lifetime

        ### Capital ###
        # Concrete and excavaction
        VEX, VWC, VSC = \
            D['Excavation [ft3]'], D['Wall concrete [ft3]'], D['Slab concrete [ft3]']
        C['Excavation'] = VEX / 27 * 8 # 27 is to convert the VEX from ft3 to yard3
        C['Wall concrete'] = VWC / 27 * 650
        C['Slab concrete'] = VSC / 27 * 350

        # Packing material
        # 195 is the cost of both LDPE and HDPE in $/m3
        C['Packing LDPE'] = 195 * D['Packing LDPE [m3]']
        C['Packing HDPE'] = 195 * D['Packing HDPE [m3]']

        # Pump
        # Note that maintenance and operating costs are included as a lumped
        # number in the biorefinery thus not included here
        # TODO: considering adding the O&M and letting user choose if to include
        pumps, building = self._cost_pump()
        C['Pumps'] = pumps
        C['Pump building'] = building
        C['Pump excavation'] = VEX / 27 * 0.3

        F_BM['Pumps'] = F_BM['Pump building'] = F_BM['Pump excavation'] = \
            1.18 * (1+0.007) # 0.007 is for  miscellaneous costs
        lifetime['Pumps'] = 15

        # Set bare module factor to 1 if not otherwise provided
        for k in C.keys():
            F_BM[k] = 1 if not F_BM.get(k) else F_BM.get(k)

        self.power_utility.rate = self.Pump.power_utility.rate


    def _cost_pump(self):
        Q_mgd, recir_ratio = self.Q_mgd, self.recir_ratio

        # Installed pump cost, this is a fitted curve
        pumps = 2.065e5 + 7.721*1e4*Q_mgd

        # Design capacity of intermediate pumps, gpm,
        # 2 is the excess capacity factor to handle peak flows
        GPMI = 2 * Q_mgd * 1e6 / 24 / 60

        # Design capacity of recirculation pumps, gpm
        GPMR = recir_ratio * Q_mgd * 1e6 / 24 / 60

        building = 0.
        for GPM in (GPMI, GPMR):
            if GPM == 0:
                N = 0
            else:
                N = 1 # number of buildings
                GPMi = GPM
                while GPMi > 80000:
                    N += 1
                    GPMi = GPM / N

            PBA = N * (0.0284*GPM+640) # pump building area, [ft]
            building += 90 * PBA

        return pumps, building


    # Util function
    @staticmethod
    def compute_COD(stream):
        return compute_stream_COD(stream)


    # TODO: now need to add pumps externally,
    # consider using algorithms in ref [2] to add influent, recirculation,
    # and effluent pumps
    @property
    def Pump(self):
        '''[:class:`~.WWTpump`] Lift pump of the filter tank.'''
        return self._Pump

    @property
    def filter_type(self):
        '''[str] Can either be "anaerobic" or "aerobic".'''
        return self._filter_type
    @filter_type.setter
    def filter_type(self, i):
        if i.lower() in ('anaerobic', 'aerobic'):
            self._filter_type = i.lower()
        else:
            raise ValueError('`filter_type` can only be "anaerobic" or "aerobic", '
                             f'not "{i}".')

    @property
    def OLR(self):
        '''[float] Organic loading rate, [kg COD/m3/hr].'''
        return self._OLR
    @OLR.setter
    def OLR(self, i):
        if i < 0:
            raise ValueError('`OLR` should be >=0, '
                             f'the input value {i} is outside the range.')
        self._OLR = i

    @property
    def HLR(self):
        '''[float] Hydraulic loading rate, [m3/m2/hr].'''
        return self._HLR
    @HLR.setter
    def HLR(self, i):
        self._HLR = i

    @property
    def d_max(self):
        '''[float] Maximum filter diameter, [m].'''
        return self._d_max
    @d_max.setter
    def d_max(self, i):
        self._d_max = i

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

    @property
    def N_filter(self):
        '''[int] Number of filter tanks.'''
        return self._N_filter

    @property
    def d(self):
        '''[float] Diameter of the filter tank, [ft].'''
        return self._d

    @property
    def D(self):
        '''[float] Depth of the filter tank, [ft].'''
        return self._D

    @property
    def Q_mgd(self):
        '''
        [float] Influent volumetric flow rate in million gallon per day, [mgd].
        '''
        return self.ins[0].F_vol*_m3_to_gal*24/1e6

    # @property
    # def Q_gpm(self):
    #     '''[float] Influent volumetric flow rate in gallon per minute, [gpm].'''
    #     return self.Q_mgd*1e6/24/60

    # @property
    # def Q_cmd(self):
    #     '''
    #     [float] Influent volumetric flow rate in cubic meter per day, [cmd].
    #     '''
    #     return self.Q_mgd *1e6/_m3_to_gal # [m3/day]

    @property
    def recir_ratio(self):
        '''
        [float] Internal recirculation ratio, will be updated in simulation
        if the originally set ratio is not adequate for the desired flow
        considering retentate and GAC (if applicable).
        '''
        return self._recir_ratio or self.ins[1].F_vol/self.ins[0].F_vol
    @recir_ratio.setter
    def recir_ratio(self, i):
        self._recir_ratio = float(i)

    @property
    def i_rm(self):
        '''[:class:`np.array`] Removal of each chemical in this filter tank.'''
        return self._i_rm

    @property
    def split(self):
        '''Component-wise split to the treated water.'''
        return self._split
    @split.setter
    def split(self, i):
        self._split = i
        self._isplit = self.chemicals.isplit(i, order=None)

    @property
    def X_decomp(self):
        '''
        [float] Fraction of the influent COD converted to biogas
        (`filter_type`=="anaerobic") or CO2 (`filter_type`=="aerobic").
        '''
        return self._X_decomp
    @X_decomp.setter
    def X_decomp(self, i):
        if not 0 <= i <= 1:
            raise ValueError('`X_decomp` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._X_decomp = i

    @property
    def X_growth(self):
        '''
        [float] Fraction of the influent COD converted to biomass growth.
        '''
        return self._X_growth
    @X_growth.setter
    def X_growth(self, i):
        if not 0 <= i <= 1:
            raise ValueError('`X_growth` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._X_growth = i

    @property
    def decomp_rxns(self):
        '''
        [:class:`tmo.ParallelReaction`] Organics to biogas (`filter_type`=="anaerobic")
        or CO2 (`filter_type`=="aerobic") reactions.
        '''
        return self._decomp_rxns

    @property
    def growth_rxns(self):
        '''
        [:class:`tmo.ParallelReaction`] Biomass (WWTsludge) growth reactions.
        '''
        return self._growth_rxns

    @property
    def organic_rm(self):
        '''[float] Overall organic (COD) removal rate.'''
        Qi, Qe = self._inf.F_vol, self.outs[1].F_vol
        Si, Se = self.compute_COD(self._inf), self.compute_COD(self.outs[1].F_vol)
        return 1 - Qe*Se/(Qi*Si)