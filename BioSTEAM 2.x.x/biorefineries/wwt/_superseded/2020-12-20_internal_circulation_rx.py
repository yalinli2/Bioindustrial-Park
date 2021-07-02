#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 10:24:11 2020

@author: yalinli_cabbi
"""



# %%

import numpy as np
from flexsolve import aitken_secant, IQ_interpolation

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



all_rx = dict.fromkeys(rx_param, 0)
b_rx = dict.fromkeys(rx_param, 0)
t_rx = dict.fromkeys(rx_param, 0)

design_dcts = {
    'all_rx': all_rx,
    'b_rx': b_rx,
    't_rx': t_rx
    }



# Upstream inputs
all_rx['Qi'] = 325 # influent volumetric flow rate, [m3/hr]
all_rx['SS'] = 5040/1e3 # influent soluble substrate conc., [kg/m3 as COD]
all_rx['X'] = 230/1e3 # influent biomass conc., [kg/m3]

def update_design():
    all_rx['Qe'] = all_rx['Qi'] - all_rx['Qw']
    for i in ('Qi', 'Qw', 'Qe'):
        b_rx[i] = all_rx[i]
    t_rx['Qi'] = t_rx['Qe'] = all_rx['Qe']
    
    all_rx['SRT'] = (b_rx['X']*b_rx['V']+t_rx['X']*t_rx['V']) / \
        (b_rx['Qw']*b_rx['V']+F_xt*t_rx['X']*t_rx['Qe'])
    b_rx['SRT'] = b_rx['X']*b_rx['V']/(b_rx['X']*b_rx['Qw']+F_xb*b_rx['X']*b_rx['Qe'])
    t_rx['SRT'] = t_rx['X']*t_rx['V']/(t_rx['X']*t_rx['Qw']+F_xt*t_rx['X']*t_rx['Qe'])
    
    for i in (all_rx, b_rx, t_rx):
        if i is all_rx:
            i['A'] = i['Qe']/i['v'] #!!! Qe or Qi?
            i['D'] = (i['A']/(np.pi/4))**(1/2)
            i['rX'] = mu_max * i['X']
            i['rSS'] = - i['rX']/Y
        else:
            i['OLR'] = i['Qi']*i['SS']/i['V']
            i['A'] = all_rx['A']
            i['D'] = all_rx['D']
            i['v'] = all_rx['v']
        i['H'] = i['V']/i['A']
        i['HtoD'] = i['H']/i['D']
        i['HRT'] = i['V']/i['Qi']
        



def get_cod_rm(mu_max=0.01, Y=0.07, F_xb=0.0032, F_xt=0.0281,
               OLR_all=1.094, V_b_to_V_t=1.086805429, v=6.27, waste_ratio=0.05):
    # Overall
    all_rx['OLR'] = OLR_all
    all_rx['v'] = v
    all_rx['Qw'] = waste_ratio*all_rx['Qi']
    all_rx['Qe'] = all_rx['Qi'] - all_rx['Qw']
    all_rx['V'] = all_rx['Qi']*all_rx['SS']/OLR_all
    
    # Bottom rx
    b_rx['V'] = all_rx['V']/(1+1/V_b_to_V_t)
    # Biomass mass balance of the bottom rx
    b_rx['X'] = (all_rx['Qi']*all_rx['X'])/(F_xb*all_rx['Qe']+all_rx['Qw']-mu_max*b_rx['V'])
    b_rx['rX'] = mu_max * b_rx['X']
    b_rx['rSS'] = - b_rx['rX']/Y
    # Substrate mass balance of the bottom rx
    b_rx['SS'] = all_rx['SS'] + b_rx['V']*b_rx['rSS']/all_rx['Qi']
    
    # Top rx
    t_rx['V'] = b_rx['V']/V_b_to_V_t
    # Biomass mass balance of the top rx
    t_rx['X'] = (F_xb*all_rx['Qe']*b_rx['X'])/(F_xt*all_rx['Qe']-mu_max*t_rx['V'])
    t_rx['rX'] = mu_max * t_rx['X']
    t_rx['rSS'] = - t_rx['rX']/Y
    # Substrate mass balance of the top rx
    t_rx['SS'] = b_rx['SS'] + t_rx['V']*t_rx['rSS']/all_rx['Qe']
    
    all_rx['COD_rm'] = 1 - (t_rx['SS']*all_rx['Qe'])/(all_rx['SS']*all_rx['Qi'])
    
    return all_rx['COD_rm']


# print(f'COD removal is {get_cod_rm():.1%}.')

# update_design()


def rm_at_OLR(OLR_all, V_b_to_V_t, v, waste_ratio, COD_rm):
    rm = get_cod_rm(OLR_all=OLR_all,
                   V_b_to_V_t=V_b_to_V_t,
                   v=v, waste_ratio=waste_ratio)
    print(rm)
    return rm-COD_rm


def solve_param_from_rm(OLR_all=None, V_b_to_V_t=None, v=None, waste_ratio=None, COD_rm=None):
    if not OLR_all:
        # all_rx['OLR'] = IQ_interpolation(f=rm_at_OLR, x0=1e-4, x1=35000/24,
        #                               xtol=1e-4, ytol=1e-4, maxiter=50,
        #                               args=(V_b_to_V_t, v, waste_ratio, COD_rm),
        #                               checkbounds=False)

        all_rx['OLR'] = aitken_secant(f=rm_at_OLR, x0=1e-4, 
                                      args=(V_b_to_V_t, v, waste_ratio, COD_rm))

        update_design()
        return all_rx['OLR']
#!!! Maybe not use solver? Just try things out...
OLR_all = solve_param_from_rm(OLR_all=None,
                    V_b_to_V_t=1.086805429,
                    v=6.27,
                    waste_ratio=0.05,
                    COD_rm=0.8877091665747813)
print(OLR_all)






# %%

# =============================================================================
# OUTDATED
# =============================================================================

# def check_constraints(max_OLR=None, min_COD_rm=None):
#     update_design()
#     max_OLR = max_OLR or 35/24
#     min_COD_rm = min_COD_rm or 0
#     if not 0<=all_rx['OLR']<=max_OLR:
#         print('Check #1 failed.')
#         raise DesignError('No feasible design for the given specifications.')
#     elif not 0<=t_rx['SS']<=b_rx['SS']<=all_rx['SS']:
#         print('Check #2 failed.')
#         raise DesignError('No feasible design for the given specifications.')
#     elif not 0<=t_rx['X']<=b_rx['X']:
#         print('Check #3 failed.')
#         raise DesignError('No feasible design for the given specifications.')
#     elif not min_COD_rm<=b_rx['COD_rm']<=all_rx['COD_rm']<= 1:
#         print(f"Check #4 failed: min_COD_rm is {min_COD_rm}, b_rx is {b_rx['COD_rm']}, all_rx['COD_rm'] is {all_rx['COD_rm']}.")
#         raise DesignError('No feasible design for the given specifications.')




# def solve_param_from_rm(OLR_all=None, V_b_to_V_t=None, v=None, waste_ratio=None, COD_rm=None):
#     if not OLR_all:
#         OLR_all0 = 1e-4
#         while OLR_all0 <=35/24:
#             print(f'\nOLR_all0 is {OLR_all0}.')
#             OLR_all = fs.IQ_interpolation(f=rm_at_OLR, x0=OLR_all0, x1=35/24,
#                                           xtol=1e-4, ytol=1e-4, maxiter=500,
#                                           args=(V_b_to_V_t, v, waste_ratio, COD_rm),
#                                           checkbounds=False)
#             diff = rm_at_OLR(OLR_all, V_b_to_V_t, v, waste_ratio, COD_rm)
#             try:
#                 assert abs(diff)<1e-4
#                 print(f'GOOD, diff is {diff}.')
#             except:
#                 print(f'diff is {diff}, absolute value larger than 1e-4.')
#                 OLR_all0 += 5/1e3
#                 print(f"DIFF not good, OLR_all is {all_rx['OLR']}, COD_rx is {all_rx['COD_rm']}.")
#                 continue
#             all_rx['OLR'] = OLR_all
#             try:
#                 check_constraints()
#                 print('Feasible design found.')
#                 return OLR_all
#             except:
#                 OLR_all0 += 5/1e3
#                 print(f"CONSTR not met, OLR_all is {all_rx['OLR']}, COD_rx is {all_rx['COD_rm']}.")

#         # all_rx['OLR'] = fs.fixed_point(f=rm_at_OLR, x=1e-4, 
#         #                               args=(V_b_to_V_t, v, waste_ratio, COD_rm))

#         print('Did not find a feasible design')
#         update_design()
        
#         return all_rx['OLR']




# rx_param = (
#     ['Qi', 'influent liquid flow', '[m3/hr]'],
#     ['Qe', 'effluent liquid flow', '[m3/hr]'],
#     'Qe', # effluent liquid flow, [m3/hr]
#     'Qw', # waste sludge, [m3/hr]
#     'Qg', # gas flow, [m3/hr]
#     'OLR', # organic loading rate, [kg/m3/hr]
#     'SS', # soluble substrate concentration, [kg/m3 as COD]
#     'rSS', # soluble substrate consumption rate (should<=0), [kg/m3/hr as COD]
#     'X', # biomass concentration (should>=0), [kg/m3]
#     'rX', # biomass growth rate, [kg/m3/hr]
#     'V', # reactor volume, [m3]
#     'A', # reactor area, [m2]
#     'H', # reactor height, [H]
#     'D', # reactor diameter, [m]
#     'HtoD', # aspect ratio
#     'v', # upflow velocity, [m/hr]
#     'HRT', # hydraulic retention time, [hr]
#     'SRT', # solid retention time, [hr]
#     'COD_rm', # COD removal
#     )


























