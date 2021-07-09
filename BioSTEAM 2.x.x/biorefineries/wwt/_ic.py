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
    OLRall : float
        Overall organic loading rate, [kg COD/m3/hr].
    CODrm : float
        COD removal.
    qw : float
        Ratio between the bottom reactor waste flow and the influent.
    Y : float
        Biomass yield, [kg biomass/kg substrate COD].
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

    # Design and operating assumptions
    _Fxb = 0.0032
    _Fxt = 0.0281
    _BD_dct = get_BD_dct(1.)

    # Related to cost algorithm
    _default_vessel_type = 'IC'
    _default_vessel_material = 'Stainless steel'
    purchase_cost_algorithms = IC_purchase_cost_algorithms

    # Heating utilities
    auxiliary_unit_names = ('heat_exchanger',)

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 method='separate', OLRall=1.25, CODrm=0.8, qw=0.05,
                 Y=0.07, mu_max=0.01, b=0.00083, V_wf=0.8,
                 vessel_type='IC', vessel_material='Stainless steel',
                 kW_per_m3=0., T=35+273.15, **kwargs):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.method = method
        self.OLRall = OLRall
        self.CODrm = CODrm
        self.qw = qw
        self.Y = Y
        self.mu_max = mu_max
        self.b = b
        self.V_wf = V_wf or self._default_V_wf
        self.vessel_type = 'IC'
        self.vessel_material = vessel_material
        self.kW_per_m3 = kW_per_m3
        self.T = T
        self.heat_exchanger = hx = bst.HXutility(None, None, None, T=T)
        self.heat_utilities = hx.heat_utilities

        # Initiate the attributes
        self._Vliq = self._Vb = self._Vt = None
        self._HRT = self._SRT = None
        self._AD_rxns = None
        # Note that the reaction conversion need to be adjusted when using
        #!!! Need to see if this reaction is good (esp. conversion of N to N2)
        self._decay_rxn = self.chemicals.WWTsludge.get_combustion_reaction()

        for k, v in kwargs.items():
            setattr(self, k, v)


    def _get_AD_rxns(self):
        '''Refresh the auto-generated biogas and growth reactions.'''
        X_growth = self.Y * self.CODrm
        X_biogas = self.CODrm - X_growth
        rxns = get_AD_rxns(self.ins[0], self.BD_dct,
                           X_biogas, X_growth, 'WWTsludge')
        return rxns



    def _run(self):
        inf = tmo.Stream()
        inf.mix_from(self.ins)
        self._inf = inf
        biogas, eff, waste  = self.outs
        inf.split_to(waste, eff, self.qw)

        biogas.T = inf.T
        biogas.phase = 'g'
        biogas.empty()

        Qi= inf.F_vol # influent volumetric flow rate, m3/hr
        Si = compute_stream_COD(inf, self.BD_dct) # influent COD, kg/m3
        Xi = inf.imass['WWTsludge'] # influent biomass concentration, kg/hr

        Qw = waste.F_vol # waste volumetric flow rate, m3/hr
        Qe = Qi - Qw # effluent volumetric flow rate, m3/hr
        St = Qi*Si*(1-self.CODrm) / Qe # effluent COD, kg/m3

        OLRall = self.OLRall # overall organic loading rate, kg COD/m3/hr
        Vliq = self._Vliq = Qi*Si / OLRall # total volume (bottom and top rxs), m3

        # Reactor design using different method
        mu_max, b, Y = self.mu_max, self.b, self.Y
        method = self.method.lower()
        run_inputs = (inf, eff, waste, biogas, mu_max, b, Y, Qi, Xi, Si, St, Vliq, Qe, Qw)
        # Returns substrate concentration in the waste flow
        if method == 'separate':
            Sw = self._run_separate(run_inputs)
        elif method == 'lumped':
            Sw = self._run_lumped(run_inputs)
        else:
            raise ValueError('`method` can only be "separated" or "lumped", '
                             f'not "{method}".')

        # Calculated reactant conversion based on the design results
        rxn_dct = {
            'waste': (Sw, waste),
            'eff': (St, eff)
            }
        for k, v in rxn_dct.items():
            X_S = (Si-v[0]) / Si # substrate conversion
            rxns = self._get_AD_rxns()
            rxns._X *= X_S
            stream_copy = v[1].copy() # to calculate decay conversion
            rxns(v[1].mol)
            # This is biomass change for mu_max-b
            net_sludge = v[1].imol['WWTsludge'] - stream_copy.imol['WWTsludge']

            stream_copy.empty()
            stream_copy.imol['WWTsludge'] = net_sludge/(mu_max-b)*b # the amount undergoes decay

            self.decay_rxn.force_reaction(stream_copy.mol)
            stream_copy.imol['O2'] = 0
            v[1].mix_from((v[1], stream_copy))
            biogas.receive_vent(v[1], accumulate=True)


    def _run_separate(self, run_inputs):
        inf, eff, waste, biogas, mu_max, b, Y, Qi, Xi, Si, St, Vliq, Qe, Qw = run_inputs
        mu_max, b, Y, Fxb, Fxt = \
            self.mu_max, self.b, self.Y, self.Fxb, self.Fxt

        Xb, Xt, Sb, Vb = sp.symbols('Xb, Xt, Sb, Vb', real=True)

        # Mass balances based on biomass/substrate changes in the bottom/top rx,
        # (0 at steady state)
        biomass_b = Qi*Xi - (Qe*Xb*Fxb+Qw*Xb) + Xb*Vb*(mu_max-b)
        biomass_t = Qe*(Fxb*Xb-Fxt*Xt) + Xt*(Vliq-Vb)*(mu_max-b)
        substrate_b = Qi*(Si-Sb) - mu_max*(Xb*Vb/Y)
        substrate_t = Qe*(Sb-St) - mu_max*((Vliq-Vb)*Xt/Y)

        parameters = (Qi, Qe, Si, St, Vliq)

        results = sp.solve(
            (sp.Eq(biomass_b, 0),
             sp.Eq(biomass_t, 0),
             sp.Eq(substrate_b, 0),
             sp.Eq(substrate_t, 0)), (Xb, Xt, Sb, Vb))

        Xb, Xt, Sb, Vb = self._filter_results('separate', parameters, results)

        Vt = Vliq - Vb # volume of the top rx, m3
        self._Vb, self._Vt = Vb, Vt
        self._HRT = Vliq / Qi # hydraulic retention time, hr
        self._SRT = (Xb*Vb+Xt*Vt) / (Qw*Xb+Fxt*Qe*Xt) # solid residence time, hr

        return Sb


    @staticmethod
    def _filter_results(method, parameters, results):
        '''Check if the solution satisfies the design constraints.'''
        Qi, Qe, Si, St, Vliq = parameters
        solutions = []
        for result in results:
            Xb, Xt, Sb, Vb = result
            Vt = Vliq - Vb
            OLRt = Qe*Sb / Vt
            OLRb = Qi*Si / Vb

            if (
                    0 <= OLRt <= OLRb and
                    0 <= St <= Sb <= Si and
                    0 <= Xt <= Xb and
                    0 <= Vb <= Vliq
                ):
                solutions.append(result)

        if len(solutions) == 0 :
            raise DesignError('No feasible design found for the given parameters.')

        elif len(solutions) >1: # find more than one solution
            Xbs = [i[1] for i in solutions]
            index = Xbs.index(min(Xbs)) # choose the one with lowest effluent biomass
            return solutions[index]


    def _run_lumped(self, run_inputs):
        inf, eff, waste, biogas, mu_max, b, Y, Qi, Xi, Si, St, Vliq, Qe, Qw = run_inputs
        rS = Qi*(Si-St) / Vliq # substrate conversion
        Xall = -rS*Y / mu_max # overall biomass concentration in the rx, kg/m3
        Xw = (Qe*Xall-Qi*Xi-Vliq*(Xall*mu_max-b)) / Qw

        self._HRT = Vliq / Qi # hydraulic retention time, hr
        self._SRT = Xall*Vliq / (Qw*Xw+Qe*Xall) # solid residence time, hr

        return St

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
    def CODrm(self):
        '''
        [float] COD removal.
        '''
        return self._CODrm
    @CODrm.setter
    def CODrm(self, i):
        if not 0<=i<=1:
            raise ValueError('`CODrm` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._CODrm = i

    @property
    def qw(self):
        '''[float] Ratio between the bottom reactor waste flow and the influent.'''
        return self._qw
    @qw.setter
    def qw(self, i):
        if not 0<=i<=1:
            raise ValueError('`qw` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._qw = i

    @property
    def Y(self):
        '''[float] Biomass yield, [kg biomass/kg substrate COD].'''
        return self._Y
    @Y.setter
    def Y(self, i):
        if i < 0:
            raise ValueError('`Y` should be >= 0, '
                             f'the input value {i} is outside the range.')
        self._Y = i

    @property
    def mu_max(self):
        '''[float] Maximum specific growth rate, [/hr].'''
        return self._mu_max
    @mu_max.setter
    def mu_max(self, i):
        if i < 0:
            raise ValueError('`mu_max` should be >= 0, '
                             f'the input value {i} is outside the range.')
        self._mu_max = i

    @property
    def b(self):
        '''[float] Specific endogenous decay coefficient, [/hr].'''
        return self._b
    @b.setter
    def b(self, i):
        if i < 0:
            raise ValueError('`b` should be >= 0, '
                             f'the input value {i} is outside the range.')
        self._b = i

    @property
    def Fxb(self):
        '''
        [float] Biomass transfer ratio from the bottom reacor to the top reactor,
        should be within [0, 1] (ideal to no retention).
        '''
        return self._Fxb
    @Fxb.setter
    def Fxb(self, i):
        if not 0<=i<=1:
            raise ValueError('`Fxb` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._Fxb = i

    @property
    def Fxt(self):
        '''
        [float] Biomass transfer ratio from the top reacor to the effluent,
        should be within [0, 1] (ideal to no retention).
        '''
        return self._Fxt
    @Fxt.setter
    def Fxt(self, i):
        if not 0<=i<=1:
            raise ValueError('`Fxt` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._Fxt = i

    @property
    def BD_dct(self):
        '''
        [dict] Biodegradability of chemicals.
        '''
        return self._BD_dct
    @BD_dct.setter
    def BD_dct(self, i):
        for k, v in i.items():
            if not 0<=v<=1:
                raise ValueError('Biodegradability should be within [0, 1], '
                                 f'the input value for chemical "{k}" is '
                                 'outside the range.')
        self._BD_dct = i

    @property
    def Vb(self):
        '''[float] Volume of the bottom reactor, [m3].'''
        if self.method == 'lumped':
            raise ValueError('Volume of the bottom reactor cannot be calcualted using the lumped method.')
        return self._Vb

    @property
    def Vt(self):
        '''[float] Volume of the top reactor, [m3].'''
        if self.method == 'lumped':
            raise ValueError('Volume of the bottom reactor cannot be calcualted using the lumped method.')
        return self._Vt

    @property
    def Vliq(self):
        '''
        [float] Total volume of the bottom and top reactor
        (gas chamber not included), [m3].
        '''
        return self._Vliq

    @property
    def HRT(self):
        '''[float] Hydraulic retention time [hr].'''
        return self._HRT

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
        return self._SRT

    @property
    def AD_rxns(self):
        '''
        [:class:`tmo.ParallelReaction`] Anaerobic digestion reactions
        (biogas production and biomass growth).
        '''
        if self._AD_rxns is None:
            self._AD_rxns = self._get_AD_rxns()
        return self._AD_rxns

    @property
    def decay_rxn(self):
        '''
        [:class:`tmo.Reaction`] Biomass endogeneous decay.

        .. note::
            Conversion is defaulted to 100%, and needs to be adjusted when used.
        '''
        return self._decay_rxn