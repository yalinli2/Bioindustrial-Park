#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 13:15:56 2021

@author: Yalin Li
"""

import sympy as sp

Qi = 325
Si = 5.04
Xi = 0.23
# COD_rm = 0.8
COD_rm = 0.9
qw = 0.05
OLR_overall = 1.25
# OLR_overall = 1.084
Vtot = 1310
mu_max = 0.01
b = 0.00083
# b = 0
Y = 0.07
F_Xb = 0.0032
F_Xt = 0.0281


Qw = qw * Qi
Qe = Qi - Qw
St = Qi*Si/Qe*(1-COD_rm)

Xb, Xt, Sb, Vb = sp.symbols('Xb, Xt, Sb, Vb', real=True)

biomass_b = Qi*Xi - (Qe*Xb*F_Xb+Qw*Xb) + Xb*Vb*(mu_max-b)
biomass_t = Qe*(F_Xb*Xb-F_Xt*Xt) + Xt*(Vtot-Vb)*(mu_max-b)
substrate_b = Qi*(Si-Sb) - mu_max*(Xb*Vb/Y)
substrate_t = Qe*(Sb-St) - mu_max*((Vtot-Vb)*Xt/Y)

# # With some shorthands
# QiXi = Qi * Xi
# QeFXbQw = Qe*F_Xb + Qw
# mu_net = mu_max - b

# QeFXb = Qe * F_Xb
# QeFXt = Qe * F_Xt
# mu_maxY = mu_max / Y

# # Organized to be
# biomass_b2 = QiXi - QeFXbQw*Xb + mu_net*Vb*Xb
# biomass_t2 = QeFXb*Xb - QeFXt*Xt + mu_net*(Vtot-Vb)*Xt
# substrate_b2 = Qi*(Si-Sb) - mu_maxY*Xb*Vb
# substrate_t2 = Qe*(Sb-St) - mu_maxY*(Vtot-Vb)*Xt




# Qw = 16
# Qe = 309
# St = 1.06





eq1 = sp.Eq(biomass_b, 0)
eq2 = sp.Eq(biomass_t, 0)
eq3 = sp.Eq(substrate_b, 0)
eq4 = sp.Eq(substrate_t, 0)

# results = sp.solve((eq1, eq2, eq3, eq4), (Xb, Xt, Sb, Vb))



# eq5 = sp.Eq(biomass_b2, 0)
# eq6 = sp.Eq(biomass_t2, 0)
# eq7 = sp.Eq(substrate_b2, 0)
# eq8 = sp.Eq(substrate_t2, 0)

sp.solve((eq1, eq2, eq3, eq4), (Xb, Xt, Sb, Vb))

# sp.solve((eq5, eq6, eq7, eq8), (Xb, Xt, Sb, Vb))

# Xb = 8.691
# Xt = 1.514
# Sb = 1.289
# Vt = 328
# Vb = 982

# x1 = (Xb*(Qe*F_Xb+Qw)-Qi*Xi)/(mu_max-b)
# x2 = F_Xb*Xb/(F_Xt/(Vtot-(Qe*F_Xb+Qw)/(mu_max-b)-Qi*Xi/Xb/(mu_max-b))-(mu_max-b)/Qe)
# Qe*(Si-St-mu_max/(Qi*Y)*x1)-mu_max/Y*x2

# # Vb
# ((Qe*Xb*F_Xb+Qw*Xb)-Qi*Xi)/(Xb*(mu_max-b))