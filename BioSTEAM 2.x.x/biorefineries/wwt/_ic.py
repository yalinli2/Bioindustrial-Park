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

import sympy as sp
import flexsolve as flx
import biosteam as bst
import thermosteam as tmo
from biosteam.exceptions import DesignError
#!!! Need to enable relative importing
from _utils import (
    get_BD_dct,
    compute_stream_COD,
    get_AD_rxns,
    IC_purchase_cost_algorithms
    )

__all__ = ('IC',)

def _check_if_relevant(attr, method, relevant_method):
    if not relevant_method in method:
        raise AttributeError(f'`{attr[1:]}` is not relevant when `method` is '
                             f'"{method}".')
    return attr


# %%

class IC(bst.MixTank):
    '''
    Internal circulation (IC) reactor for anaerobic digestion (AD),
    including a high-rate bottom reactor for rapid organic removal and
    a low-rate top reactor for polishing.
    Both reactors are similar to upflow anaerobic blanket reactor (UASB).

    Design of the reactor follows steps described in [1]_
    (assuming steady state and pseudo-zeroth-order kinetics),
    where two methods are used based on Irizar et al.[2]_ and
    Tchobanoglous et al.[3]_.

    Parameters
    ----------
    method : str
        Either "separate" to design the bottom and top reactors separately as in [2]_,
        or "lumped" to design the entire IC reactor as a black box following [3]_.

        In "separate" method, design parameters include:
            - OLRall, biodegradability, Y, q_Qw, mu_max, b, Fxt, Fxb

        In "lumped" method, design parameters include:
            - OLRall, biodegradability, Y, q_Qw ("lumped-q_Qw") or q_Xw ("lumped-q_Xw")
    OLRall : float
        Overall organic loading rate, [kg COD/m3/hr].
    biodegradability : dict
        Chemical removal
        (chemicals not in the dict are conserved in this reactor).
    Y : float
        Biomass yield, [kg biomass/kg substrate COD].
    q_Qw : float
        Ratio between the bottom reactor waste flow and the influent.
    q_Xw : float
        Ratio between the biomass concentration in the reactor and the waste flow.
    mu_max : float
        Maximum specific growth rate, [/hr].
    b : float
        Specific endogenous decay coefficient, [/hr].
    V_wf : float
        Fraction of working volume over total volume.
    vessel_type : str
        Can be "IC" to use the reactor size constraints according to [1]_,
        or "Conventional" based on :class:`biosteam.MixTank`
        (much smaller tank size, not recommended).
    vessel_material : str
        Vessel material.
    kW_per_m3 : float
        Electricity requirement per unit volume, [kW/m^3].
        Default to 0 as IC reactors realizes mixing through internal circulation
        caused by the rising force of the generated biogas.
    T : float
        Temperature of the reactor.
    kwargs : dict
        Other keyword arguments (e.g., Fxb, Fxt).

    References
    ----------
    .. [1] Kontos, G. A. Advanced Anaerobic Treatment for Energy Recovery and
    Improved Process Economics in the Management of Biorefinery Wastewaters,
    University of Illinois at Urbana-Champaign, Champaign, IL, 2021.

    .. [2] Irizar et al., Model-Based Design of a Software Sensor for Real-Time
    Diagnosis of the Stability Conditions in High-Rate Anaerobic Reactors –
    Full-Scale Application to Internal Circulation Technology.
    Water Research 2018, 143, 479–491.
    `<https://doi.org/10.1016/j.watres.2018.06.055>`_.

    .. [3] Tchobanoglous et al., Wastewater Engineering: Treatment and Resource Recovery,
    5th ed.; McGraw-Hill Education: New York, 2013.
    '''
    _N_ins = 2
    _N_outs = 3
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True

    # Assumptions
    _q_Qw = 0.05
    _q_Xw = 1.5
    _mu_max = 0.01
    _b = 0.00083
    _Fxb = 0.0032
    _Fxt = 0.0281

    # Related to cost algorithm
    _default_vessel_type = 'IC'
    _default_vessel_material = 'Stainless steel'
    purchase_cost_algorithms = IC_purchase_cost_algorithms

    # Heating utilities
    auxiliary_unit_names = ('heat_exchanger',)

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 method='lumped-Xw', OLRall=1.25, biodegradability={}, Y=0.07,
                 vessel_type='IC', vessel_material='Stainless steel',
                 V_wf=0.8, kW_per_m3=0., T=35+273.15, **kwargs):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.method = method
        self.OLRall = OLRall
        self.biodegradability = \
            biodegradability if biodegradability else get_BD_dct(self.chemicals)
        self.Y = Y
        self.V_wf = V_wf or self._default_V_wf
        self.vessel_type = 'IC'
        self.vessel_material = vessel_material
        self.kW_per_m3 = kW_per_m3
        self.T = T

        # Initiate the attributes
        self.heat_exchanger = hx = bst.HXutility(None, None, None, T=T)
        self.heat_utilities = hx.heat_utilities
        self._AD_rxns = None

        # Note that the reaction conversion need to be adjusted when using
        self._decay_rxn = self.chemicals.WWTsludge.get_combustion_reaction()

        for k, v in kwargs.items():
            setattr(self, k, v)


    def get_AD_rxns(self):
        '''Refresh the auto-generated biogas and growth reactions.'''
        rxns = get_AD_rxns(self.ins[0], self.biodegradability,
                           1-self.Y, self.Y, 'WWTsludge')
        self._i_rm = rxns.X_net.data
        return rxns


    def _run(self):
        inf = tmo.Stream()
        inf.mix_from(self.ins)
        self._inf = inf.copy()
        biogas, eff, waste  = self.outs

        # Initiate the streams
        biogas.T = eff.T = waste.T = self.T
        biogas.phase = 'g'
        biogas.empty()

        method = self.method.lower()
        if method == 'lumped-q_xw':
            self.q_Qw = flx.IQ_interpolation(
                f=self._q_Qw_obj_func,
                x0=0., x1=1., xtol=1e-5, ytol=1e-5,
                checkbounds=False)

        breakpoint()
        Qi, Si, Xi, Y, Qw = self.Qi, self.Si, self.Xi, self.Y, self.Qw

        inf.split_to(waste, eff, self.q_Qw)
        AD_rxns = self.AD_rxns
        AD_rxns(inf.mol)



        if method == 'separate':
            run_inputs = (Qi, Si, Xi, self.Vliq, Y, Qw,
                          self.mu_max, self.b, self.Fxb, self.Fxt)
            Xb, Xe = self._run_separate(run_inputs)

            rxn_dct = {
                'waste': (self.Sw, waste),
                'eff': (self.Se, eff)
                }
            for k, v in rxn_dct.items():
                # to calculate decay conversion
                stream_copy = v[1].copy()
                AD_rxns(v[1].mol)
                # This is biomass change for mu_max-b
                net_sludge = v[1].imol['WWTsludge'] - stream_copy.imol['WWTsludge']

                stream_copy.empty()
                # Adjust to the amount undergoes decay
                stream_copy.imol['WWTsludge'] = \
                    net_sludge/(self.mu_max-self.b)*self.b

                self.decay_rxn.force_reaction(stream_copy.mol)
                stream_copy.imol['O2'] = 0
                v[1].mix_from((v[1], stream_copy))

        else:
            Xw = (Qi*Xi+Qi*(Si-self.Se)*Y-(Qi-Qw)*Xe)/Qw
            Xe = inf.imass['WWTsludge'] - Xw
            AD_rxns(waste.mol)
            AD_rxns(eff.mol)

        # Changes due to biomass settling
        eff.imass['WWTsludge'] = Xe
        waste.imass['WWTsludge'] = Xw

        biogas.receive_vent(eff, accumulate=True)
        biogas.receive_vent(waste, accumulate=True)


    def _q_Qw_obj_func(self, q_Qw):
        eff_temp = self._inf.copy()
        self.AD_rxns(eff_temp.mol)
        Se = compute_stream_COD(eff_temp)

        Qi, Si, Xi, Xe = self.Qi, self.Si, self.Xi, self.Xe
        Qw = Qi * q_Qw
        Xw = (Qi*Xi+Qi*(Si-Se)*self.Y-(Qi-Qw)*Xe)/Qw

        return Xw/Xe - self.q_Xw


    def _run_separate(self, run_inputs):
        Qi, Si, Xi, Qe, Se, Vliq, Y, Qw, mu_max, b, Fxb, Fxt = run_inputs

        Xb, Xe, Sb, Vb = sp.symbols('Xb, Xe, Sb, Vb', real=True)

        # Mass balances based on biomass/substrate changes in the bottom/top rx,
        # (0 at steady state)
        biomass_b = Qi*Xi - (Qe*Xb*Fxb+Qw*Xb) + Xb*Vb*(mu_max-b)
        biomass_t = Qe*(Fxb*Xb-Fxt*Xe) + Xe*(Vliq-Vb)*(mu_max-b)
        substrate_b = Qi*(Si-Sb) - mu_max*(Xb*Vb/Y)
        substrate_t = Qe*(Sb-Se) - mu_max*((Vliq-Vb)*Xe/Y)

        parameters = (Qi, Qe, Si, Se, Vliq)

        results = sp.solve(
            (sp.Eq(biomass_b, 0),
             sp.Eq(biomass_t, 0),
             sp.Eq(substrate_b, 0),
             sp.Eq(substrate_t, 0)), (Xb, Xe, Sb, Vb))

        Xb, Xe, Sb, Vb = self._filter_results('separate', parameters, results)

        Vt = Vliq - Vb # volume of the top rx, m3
        self._Vb, self._Vt = Vb, Vt
        return Xb, Xe


    @staticmethod
    def _filter_results(method, parameters, results):
        '''Check if the solution satisfies the design constraints.'''
        Qi, Qe, Si, Se, Vliq = parameters
        solutions = []
        for result in results:
            Xb, Xe, Sb, Vb = result
            Vt = Vliq - Vb
            OLRt = Qe*Sb / Vt
            OLRb = Qi*Si / Vb

            if (
                    0 <= OLRt <= OLRb and
                    0 <= Se <= Sb <= Si and
                    0 <= Xe <= Xb and
                    0 <= Vb <= Vliq
                ):
                solutions.append(result)

        if len(solutions) == 0 :
            raise DesignError('No feasible design found for the given parameters.')

        elif len(solutions) >1: # find more than one solution
            Xbs = [i[1] for i in solutions]
            index = Xbs.index(min(Xbs)) # choose the one with lowest effluent biomass
            return solutions[index]


    _units = {
        'HRT': 'hr',
        'SRT': 'hr',
        'Single reactor liquid volume': 'm3',
        'Bottom reactor volume': 'm3',
        'Top reactor volume': 'm3',
        'Gas chamber volume': 'm3'
        }
    def _design(self):
        D = self.design_results
        D['HRT'] = D['Residence time'] = self.HRT
        D['SRT'] = self.SRT
        D['Total volume'] = self.Vliq / self.V_wf #!!! if no gas headspace, then change it to Vtot
        D['Total liquid volume'] = self.Vliq
        if self.method == 'separate':
            D['Bottom reactor volume'] = self.Vb
            D['Top reactor volume'] = self.Vt

    def _cost(self):
        bst.MixTank._cost(self)

        inf = self._inf
        H_at_T = inf.thermo.mixture.H(mol=inf.mol, phase='l', T=self.T, P=101325)
        duty = -(inf.H - H_at_T)
        self.heat_exchanger.simulate_as_auxiliary_exchanger(duty, inf)


    @property
    def method(self):
        '''[str] Design method, can be "separate", "lumped-q_Qw", or "lumped-q_Xw".'''
        return self._method
    @method.setter
    def method(self, i):
        if not i.lower() in ('separate', 'lumped-q_qw', 'lumped-q_xw'):
            raise ValueError('`method` can only be "separated", "lumped-q_Qw", '
                             f'or "lumped-q_Xw" not "{i}".')
        self._method = i

    @property
    def Qi(self):
        '''[float] Influent volumetric flow rate, [m3/hr].'''
        return self._inf.F_vol

    @property
    def Qe(self):
        '''[float] Effluent volumetric flow rate, [m3/hr].'''
        return self.outs[1].F_vol

    @property
    def Qw(self):
        '''[float] Waste flow volumetric flow rate, [m3/hr].'''
        return self.outs[2].F_vol

    @property
    def Si(self):
        '''
        [float] Influent substrate (i.e., biodegradable chemicals)
        concentration, [kg/m3].
        '''
        return compute_stream_COD(self._inf)

    @property
    def Se(self):
        '''
        [float] Effluent substrate (i.e., biodegradable chemicals)
        concentration, [kg/m3].
        '''
        return compute_stream_COD(self.outs[1])

    @property
    def Sw(self):
        '''
        [float] Waste flow substrate (i.e., biodegradable chemicals)
        concentration, [kg/m3].
        '''
        return compute_stream_COD(self.outs[2])

    @property
    def organic_rm(self):
        '''[float] Overall organic removal rate.'''
        return 1 - self.Qe*self.Se/(self.Qi*self.Si)

    @property
    def Xi(self):
        '''[float] Effluent biomass (i.e., `WWTsludge`) concentration, [kg/m3].'''
        return self._inf.imass['WWTsludge']/self._inf.F_vol

    @property
    def Xe(self):
        '''[float] Influent biomass (i.e., `WWTsludge`) concentration, [kg/m3].'''
        return self.outs[1].imass['WWTsludge']/self.outs[1].F_vol

    @property
    def Xw(self):
        '''[float] Waste flow biomass (i.e., `WWTsludge`) concentration, [kg/m3].'''
        return self.outs[2].imass['WWTsludge']/self.outs[2].F_vol

    @property
    def OLRall(self):
        '''
        [float] Overall organic loading rate, [kg COD/m3/hr].
        '''
        return self._OLRall
    @OLRall.setter
    def OLRall(self, i):
        if i < 0:
            raise ValueError('`OLRall` should be >=0, '
                             f'the input value {i} is outside the range.')
        self._OLRall = i

    @property
    def biodegradability(self):
        '''
        [float of dict] Biodegradability of chemicals,
        when shown as a float, all biodegradable chemicals are assumped to have
        the same degradability.
        '''
        return self._biodegradability
    @biodegradability.setter
    def biodegradability(self, i):
        if isinstance(i, float):
            if not 0<=i<=1:
                raise ValueError('`biodegradability` should be within [0, 1], '
                                 f'the input value {i} is outside the range.')
            self._biodegradability = i
            return

        for k, v in i.items():
            if not 0<=v<=1:
                raise ValueError('`biodegradability` should be within [0, 1], '
                                 f'the input value for chemical "{k}" is '
                                 'outside the range.')
        self._biodegradability = i

    @property
    def i_rm(self):
        '''[:class:`np.array`] Removal of each chemical in this reactor.'''
        return self._i_rm

    @property
    def Y(self):
        '''
        [float] Biomass yield, [kg biomass/kg substrate COD].
        .. note::
            This yield is considered as the "synthesis"

        '''
        return self._Y
    @Y.setter
    def Y(self, i):
        if i < 0:
            raise ValueError('`Y` should be >= 0, '
                             f'the input value {i} is outside the range.')
        self._Y = i

    @property
    def q_Qw(self):
        '''[float] Ratio between the bottom reactor waste flow and the influent.'''
        return self._q_Qw
    @q_Qw.setter
    def q_Qw(self, i):
        if not 0<=i<=1:
            raise ValueError('`q_Qw` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._q_Qw = i

    @property
    def q_Xw(self):
        '''
        [float] Ratio between the biomass concentration in the reactor and the waste flow,
        only relevant when either of the "lumped" method is used.
        '''
        if self.method == 'lumped-q_Qw':
            return self.Xw / self.Xe
        return _check_if_relevant(self._q_Xw, self.method, 'lumped')
    @q_Xw.setter
    def q_Xw(self, i):
        _check_if_relevant(self._q_Xw, self.method, 'lumped')
        if not i>=1:
            raise ValueError('`q_Xw` should be >=1, '
                             f'the input value {i} is outside the range.')
        self._q_Xw = i

    @property
    def mu_max(self):
        '''
        [float] Maximum specific growth rate, [/hr],
        only relevant when the "separate" method is used.
        '''
        return _check_if_relevant(self._mu_max, self.method, 'separate')
    @mu_max.setter
    def mu_max(self, i):
        _check_if_relevant(self._mu_max, self.method, 'separate')
        if i < 0:
            raise ValueError('`mu_max` should be >= 0, '
                             f'the input value {i} is outside the range.')
        self._mu_max = i

    @property
    def b(self):
        '''
        [float] Specific endogenous decay coefficient, [/hr],
        only relevant when the "separate" method is used.
        '''
        return _check_if_relevant(self._b, self.method, 'separate')
    @b.setter
    def b(self, i):
        _check_if_relevant(self._b, self.method, 'separate')
        if i < 0:
            raise ValueError('`b` should be >= 0, '
                             f'the input value {i} is outside the range.')
        self._b = i

    @property
    def Fxb(self):
        '''
        [float] Biomass transfer ratio from the bottom reacor to the top reactor,
        should be within [0, 1] (ideal to no retention),
        only relevant when the "separate" method is used.
        '''
        return _check_if_relevant(self._Fxb, self.method, 'separate')
    @Fxb.setter
    def Fxb(self, i):
        _check_if_relevant(self._Fxb, self.method, 'separate')
        if not 0<=i<=1:
            raise ValueError('`Fxb` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._Fxb = i

    @property
    def Fxt(self):
        '''
        [float] Biomass transfer ratio from the top reacor to the effluent,
        should be within [0, 1] (ideal to no retention),
        only relevant when the "separate" method is used.
        '''
        return _check_if_relevant(self._Fxt, self.method, 'separate')
    @Fxt.setter
    def Fxt(self, i):
        _check_if_relevant(self._Fxt, self.method, 'separate')
        if not 0<=i<=1:
            raise ValueError('`Fxt` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._Fxt = i

    @property
    def Vb(self):
        '''
        [float] Volume of the bottom reactor, [m3],
        only relevant when the "separate" method is used.
        '''
        return _check_if_relevant(self._Vb, self.method, 'separate')

    @property
    def Vt(self):
        '''
        [float] Volume of the top reactor, [m3],
        only relevant when the "separate" method is used.
        '''
        return _check_if_relevant(self._Vt, self.method, 'separate')

    @property
    def Vliq(self):
        '''
        [float] Total volume of the liquid
        (gas headspace not included), [m3].
        '''
        return self.Qi*self.Si / self.OLRall

    @property
    def HRT(self):
        '''[float] Hydraulic retention time [hr].'''
        return self.Vliq / self.Qi

    @property
    def tau(self):
        '''
        [float] Reactor residence time, [hr]
        (same as the hydraulic retention time, HRT).
        '''
        return self.HRT

    @property
    def SRT(self):
        '''[float] Solid residence time [hr].'''
        if self.method == 'separate':
            return (self.Xb*self.Vb+self.Xe*self.Vt)/(self.q_Qw*self.Xb+self.Fxt*self.Qe*self.Xe)

        return  self.Xe*self.Vliq / (self.q_Qw*self.q_Xw+self.Qe*self.Xe)

    @property
    def AD_rxns(self):
        '''
        [:class:`tmo.ParallelReaction`] Anaerobic digestion reactions
        (biogas production and biomass growth).
        '''
        if self._AD_rxns is None:
            self._AD_rxns = self.get_AD_rxns()
        return self._AD_rxns

    @property
    def decay_rxn(self):
        '''
        [:class:`tmo.Reaction`] Biomass endogeneous decay,
        only relevant when the "separate" method is used.

        .. note::
            Conversion is defaulted to 100%, and needs to be adjusted when used.
        '''
        return _check_if_relevant(self._decay_rxn, self.method, 'separate')