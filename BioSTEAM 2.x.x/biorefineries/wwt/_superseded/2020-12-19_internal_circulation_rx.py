#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 10:24:11 2020

@author: yalinli_cabbi
"""



# %%

import numpy as np

# Upstream inputs
Q_inf = 325 # influent volumetric flow rate, [m3/hr]
SS_inf = 5040/1e3 # influent soluble substrate conc., [kg/m3 as COD]
X_inf = 230/1e3 # influent biomass conc., [kg/m3]
upstream_inputs = (Q_inf, SS_inf, X_inf)


# Assumptions based on literature
mu_max = 0.01 # maximum specific growth rate, [/hr]
Y = 0.07 # biomass yield, [kg biomass/kg substrate as COD]
F_xb = 0.0032 # biomass transfer ratio from bottom to top, 0-1 (ideal to no retention)
F_xt = 0.0281 # biomass transfer ratio from the top rx to effluent
lit_assumptions = (mu_max, Y, F_xb, F_xt)


# Design parameters
OLR_b = 2100/1e3 # organic loading rate of the bottom rx, [kg/hr/m3]
OLR_t = 1000/1e3 # organic loading rate of the top rx, [kg/hr/m3]
HtoD_b = 2 # heigh-to-diameter ratio of the bottom rx
waste_ratio = 0.05 # sludge-to-influent flow rate ratio
design_param = (OLR_b, OLR_t, HtoD_b, waste_ratio)
# COD_rm = 0.821230158730159 # COD removal
# design_param = (COD_rm, OLR_b, OLR_t, HtoD_b)

rx_param = (
    'Qi', # influent liquid flow
    'Qe', # effluent liquid flow
    'Qw', # waste sludge, [m3/hr]
    'Qg', # gas flow, [m3/hr]
    'OLR', # organic loading rate, [kg/m3/hr]
    'SS', # soluble substrate concentration, [kg/m3 as COD]
    'rSS', # soluble substrate consumption rate (should<=0), [kg/m3/hr as COD]
    'X', # biomass concentration (should>=0), [kg/m3]
    'rX', # biomass growth rate, [kg/m3/hr]
    'V', # reactor volume, [m3]
    'A', # reactor area, [m2]
    'H', # reactor height, [H]
    'D', # reactor diameter, [m]
    'HtoD', # aspect ratio
    'v', # upflow velocity, [m/hr]
    'HRT', # hydraulic retention time, [hr]
    'SRT', # solid retention time, [hr]
    )



all_rx = dict.fromkeys(rx_param)
b_rx = dict.fromkeys(rx_param)
t_rx = dict.fromkeys(rx_param)

design_dcts = {
    'all_rx': all_rx,
    'b_rx': b_rx,
    't_rx': t_rx
    }



# =============================================================================
# Calculations, should be iteratively run to solve for waste_ratio
# =============================================================================

def sol_cod_rm(upstream_inputs, lit_assumptions,
               OLR_b, OLR_t, HtoD_b, waste_ratio):
    
    Q_inf, SS_inf, X_inf = upstream_inputs
    mu_max, Y, F_xb, F_xt = lit_assumptions
    
    design = {}
    # Bottom rx
    design['Q_w'] = Q_w = waste_ratio * Q_inf # waste sludge flow rate, [m3/hr]
    design['Q_b'] = Q_b = Q_inf - Q_w # effluent flow rate, [m3/hr]
    design['V_b'] = V_b = Q_inf * SS_inf / OLR_b # volume of the bottom rx, [m3]
    design['D_b'] = D_b = (V_b/(np.pi/4*HtoD_b))**(1/3) # diameter of the bottom rx, [m]
    design['A'] = A = np.pi/4 * (D_b**2) # rx area (same for bottom and top rxs), [m2]
    design['v'] = Q_b / A # upflow velocity, [m]
    design['H_b'] = D_b * HtoD_b # height of the bottom rx, [m]
    # HRT?
    # SRT?
    
    # Solve bottom rx biomass conc. based on biomass mass balance of the rx, [kg/m3]
    design['X_b'] = X_b = (Q_inf*X_inf)/(F_xb*Q_b+Q_w-mu_max*V_b)
    design['r_xb'] = r_xb = mu_max * X_b # biomass growth rate in the bottom rx, [kg/m3]
    design['r_ssb'] = r_ssb = - r_xb/Y # substrate consumption rate in the bottom rx, [kg/m3 as COD]
    
    # Solve bottom rx soluble substrate conc. based on substrate mass balance of the rx, [kg/m3 as COD]
    design['SS_b'] = SS_b = SS_inf + V_b*r_ssb/Q_inf
    
    # Top rx
    design['Q_t'] = Q_t = Q_b
    design['V_t'] = V_t = Q_b * SS_b / OLR_t
    design['D_t'] = D_t = D_b
    design['H_t'] = H_t = V_t / A
    design['HtoD_t'] = H_t / D_t
    
    # Solve top rx biomass conc. based on biomass mass balance of the rx, [kg/m3]
    design['X_t'] = X_t = (F_xb*Q_b*X_b)/(F_xt*Q_t-mu_max*V_t)
    design['r_xt'] = r_xt = mu_max * X_t
    design['r_sst'] = r_sst = - r_xt/Y
    
    # Solve top rx soluble substrate conc. based on substrate mass balance of the rx, [kg/m3 as COD]
    design['SS_t'] = SS_t = SS_b + V_t*r_sst/Q_b
    
    # Calculate COD removal
    design['COD_rm'] = 1 - (Q_t*SS_t)/(Q_inf*SS_inf)
    
    return design

# design = get_cod_removal(upstream_inputs, lit_assumptions, *design_param)













design = sol_cod_rm(upstream_inputs, lit_assumptions, *design_param)


COD_rm = design['COD_rm']

print(f'COD removal is {design["COD_rm"]:.1%}.')

# %%

import numpy as np
from sympy import symbols, Eq, solve, Symbol

# Assumptions based on literature
mu_max = 0.01 # maximum specific growth rate, [/hr]
Y = 0.07 # biomass yield, [kg biomass/kg substrate as COD]
F_xb = 0.0032 # biomass transfer ratio from bottom to top, 0-1 (ideal to no retention)
F_xt = 0.0281 # biomass transfer ratio from the top rx to effluent
lit_assumptions = (mu_max, Y, F_xb, F_xt)


rx_param = (
    'Qi', # influent liquid flow
    'Qe', # effluent liquid flow
    'Qw', # waste sludge, [m3/hr]
    'Qg', # gas flow, [m3/hr]
    'OLR', # organic loading rate, [kg/m3/hr]
    'SS', # soluble substrate concentration, [kg/m3 as COD]
    'rSS', # soluble substrate consumption rate (should<=0), [kg/m3/hr as COD]
    'X', # biomass concentration (should>=0), [kg/m3]
    'rX', # biomass growth rate, [kg/m3/hr]
    'V', # reactor volume, [m3]
    'A', # reactor area, [m2]
    'H', # reactor height, [H]
    'D', # reactor diameter, [m]
    'HtoD', # aspect ratio
    'v', # upflow velocity, [m/hr]
    'HRT', # hydraulic retention time, [hr]
    'SRT', # solid retention time, [hr]
    )



all_rx = dict.fromkeys(rx_param)
b_rx = dict.fromkeys(rx_param)
t_rx = dict.fromkeys(rx_param)

design_dcts = {
    'all_rx': all_rx,
    'b_rx': b_rx,
    't_rx': t_rx
    }


for i in all_rx.keys(): all_rx[i] = symbols(i+'_all')
for i in b_rx.keys(): b_rx[i] = symbols(i+'_b')
for i in t_rx.keys(): t_rx[i] = symbols(i+'_t')

eqs = []
# Equations within a reactor
for i in all_rx, b_rx, t_rx:
    eqs.extend([
        Eq(i['Qi'], i['Qe']+i['Qw']),
        Eq(i['Qg'], 0), # gas production updated later
        Eq(i['OLR'], i['Qi']*i['SS']/i['V']),
        Eq(i['rX'], mu_max*i['X']),
        Eq(i['rSS'], Y*i['rX']),
        Eq(i['A'], np.pi/4*(i['D']**2)),
        Eq(i['V'], i['A']*i['H']),
        Eq(i['v'], i['Qi']/i['A']),
        Eq(i['HtoD'], i['H']/i['D']),
        Eq(i['HRT'], i['V']/i['Qi']),
        ])
    i['Qg'] = 0
    

# Connections and mass balances between reactors
eqs.extend([
    Eq(all_rx['Qw'], b_rx['Qw']),
    Eq(t_rx['Qw'], 0),
    Eq(all_rx['Qe'], b_rx['Qe']),
    Eq(b_rx['Qe'], t_rx['Qe']),
    Eq(all_rx['V'], b_rx['V']+t_rx['V']),
    Eq(all_rx['A'], b_rx['A']),
    Eq(all_rx['A'], t_rx['A']),
    Eq(all_rx['SRT'],
       (b_rx['X']*b_rx['V']+t_rx['X']*t_rx['V']) / \
           (b_rx['Qw']*b_rx['V']+F_xt*t_rx['X']*t_rx['Qe'])),
    Eq(b_rx['SRT'],
       b_rx['X']*b_rx['V']/(b_rx['X']*b_rx['Qw']+F_xb*b_rx['X']*b_rx['Qe'])),
    Eq(t_rx['SRT'],
       t_rx['X']*t_rx['V']/(t_rx['X']*t_rx['Qw']+F_xt*t_rx['X']*t_rx['Qe'])),
    # Bottom rx SS
    Eq(all_rx['SS']*b_rx['Qi']+b_rx['rSS']*b_rx['V'], b_rx['SS']*b_rx['Qi']),
    # Top rx SS
    Eq(b_rx['SS']*t_rx['Qi']+t_rx['rSS']*t_rx['V'], t_rx['SS']*t_rx['Qi']),
    # Bottm rx biomass
    Eq(all_rx['X']*b_rx['Qi']+t_rx['rX']*b_rx['V'],
       F_xb*b_rx['X']*b_rx['Qe']+b_rx['X']*b_rx['Qw']),
    # Top rx biomass
    Eq(F_xb*b_rx['X']*b_rx['Qe']+t_rx['rX']*t_rx['V'], F_xt*b_rx['X']*t_rx['Qe'])
    ])


    

# COD removal
COD_rm = symbols('COD_rm')
eqs.append(
    Eq(COD_rm, 1-(t_rx['Qe']*t_rx['SS'])/(all_rx['Qi']*all_rx['SS']))
    )



# Upstream inputs
Q_inf = 325 # influent volumetric flow rate, [m3/hr]
SS_inf = 5040/1e3 # influent soluble substrate conc., [kg/m3 as COD]
X_inf = 230/1e3 # influent biomass conc., [kg/m3]
upstream_inputs = (Q_inf, SS_inf, X_inf)


def set_upstream_inputs(Q, SS, X):
    eqs.extend([
        Eq(all_rx['Qi'], Q),
        Eq(all_rx['SS'], SS),
        Eq(all_rx['X'], X)
        ])

set_upstream_inputs(*upstream_inputs)




# Design parameters
OLR_b = 2100/1e3 # organic loading rate of the bottom rx, [kg/hr/m3]
OLR_t = 1000/1e3 # organic loading rate of the top rx, [kg/hr/m3]
HtoD_b = 2 # heigh-to-diameter ratio of the bottom rx
Qw_inf = 0.05*Q_inf # sludge-to-influent flow rate ratio

design_param = {
    'OLR_b': 2100/1e3,
    'OLR_t': 1000/1e3,
    'HtoD_b': 2,
    'Qw_all': 0.05*Q_inf
    }

for i in design_param.keys():
    name = i.split('_')
    param, rx = name
    var = design_dcts[rx+'_rx'][param]
    eqs.append(
        Eq(var, design_param[i])
        )




unknowns = [*[i for i in all_rx.values() if isinstance(i, Symbol)],
            *[i for i in b_rx.values() if isinstance(i, Symbol)],
            *[i for i in t_rx.values() if isinstance(i, Symbol)],
            COD_rm]

solve(eqs, unknowns)





























