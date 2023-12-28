using Plots
using CoolProp

P = (1e5 - 1e4):1e2:1.1e5

U_liq = [PropsSI("U", "P", pi, "Q", 0, "Water") for pi in P]
U_vap = [PropsSI("U", "P", pi, "Q", 1, "Water") for pi in P]

# plot(P/1e2, U_vap/1e3)
# plot(P/1e2, U_liq/1e3)
# plot!(P/1e2, (U_vap - U_liq)/1e3)

# T = (5 + 273.15):5:(90 + 273.15)
# U_liq = [PropsSI("U", "T", Ti, "P", 1e5, "Water") for Ti in T]
# plot(T .- 273.15, U_liq)
# Cv = 4184 
# T0 = 273.16
# plot!(T .- 273.15, Cv*(T .- T0))

ρ_vap = [PropsSI("D", "P", pi, "Q", 1, "Water") for pi in P]
plot(P/1e2, ρ_vap)
T_adjust = 170 + 273.15
R = 8.314/0.018
plot!(P/1e2, P/(R*T_adjust) .+ 0.1)