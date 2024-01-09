# Testing different correlations for the volumetric/bulk interfacial heat transfer coefficient H_if

using Plots
using CoolProp

P_ref = 1e5
T_ref = 70 + 273.15

μ   = PropsSI("VISCOSITY", "P", P_ref, "T", T_ref, "Water")

kf  = PropsSI("CONDUCTIVITY", "P", P_ref, "Q", 0, "Water")

kg  = PropsSI("CONDUCTIVITY", "P", P_ref, "Q", 1, "Water")
Cpg = PropsSI("CPMASS",       "P", P_ref, "Q", 1, "Water")
ρg  = PropsSI("D",            "P", P_ref, "Q", 1, "Water")
αg  = kg/(ρg*Cpg)

σ   = PropsSI("surface_tension", "P", P_ref, "Q", 1, "Water")

# Bubbly superheated liquid

db = 1e-3
α_bub = maximum([αg, 1e-5])
agf = 3.6*α_bub/db
F1 = minimum([0.001, α_bub])/α_bub
F3 = 1
Re = 

H_interface = kf/db*(2.0 + 0.74*Re^0.5)*agf*F2*F3

# Bubbly Subcooled Liquid

H_interface = F3*F5*hfg*ρg*ρf*α_bub/(ρf - ρg)

# Own analysis

Db = 1e-3
k = 0.67
h_if = k/Db*Nu

Re = (1000*2*Db)/(0.004)
Pr = PropsSI("Prandtl", "P", P_ref, "T", T_ref, "Water")
Nu = 2 + (0.4*Re^(1/2) + 0.06*Re^(2/3)*Pr^0.4)

H_interface = 6*k*Nu/Db^2