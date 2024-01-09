using Plots
using CoolProp

P = (1e5 - 1e4):1e2:1.1e5

# U_vap = [PropsSI("U", "P"jpi, "Q", 1, "Water") for pi in P]
# U_liqs = [PropsSI("U", "P", pi, "Q", 0, "Water") for pi in P]

# plot(P/1e2, U_vap/1e3, label="Vapor")
# plot!(P/1e2, iq/1e3, label="Liquid")
# plot!(P/1e2, (U_vap - U_liq)/1e3)

T = (20 + 273.15):0.5:(99.99 + 273.15)
U_liqs = [PropsSI("U", "T", Ti, "P", 1e5, "Water") for Ti in T]
plot(T .- 273.15, U_liqs)
Cv = 4184 
T0 = 273.16
plot!(T .- 273.15, Cv*(T .- T0))
hline!([PropsSI("U","P",1e5,"Q",0,"water")])
vline!([100])


# ρ_vap = [PropsSI("D", "P", pi, "Q", 1, "Water") for pi in P]
# plot(P/1e3, ρ_vap)

# p1 = 90e3
# p2 = 110e3
# U1 = PropsSI("U", "P", p1, "Q", 1, "Water")
# U2 = PropsSI("U", "P", p2, "Q", 1, "Water")
# U_linear(p) = (U2 - U1)/(p2 - p1)*(p - p1) + U1
# plot(P/1e3, U_vap)
# plot!(P/1e3, U_linear.(P))

P = 90e3:0.5e3:110e3
U_liqs = [PropsSI("U", "P", pi, "Q", 0, "Water") for pi in P]
plot(P, U_liqs)
