using Plots
using CoolProp
using Roots


Pc = 22.064e6
Tc = 647.096
ρc = 322.0

a1 = -7.85951783; a2 = 1.84408259; a3 = -11.7866497; a4 = 22.6807411; a5 = -15.9618719; a6 = 1.80122502 
function Psat_IAWPS(T)
    θ = 1 - T/Tc
    return Pc*exp(Tc/T*(a1*θ + a2*θ^1.5 + a3*θ^3 + a4*θ^3.5 + a5*θ^4 + a6*θ^7.5))
end

c1 = -2.03150240; c2= -2.68302940; c3 = -5.38626492; c4 = -17.2991605; c5 = -44.7586581; c6 = -63.920106
function ρsat_IAWPS(T)
    θ = 1 - T/Tc
    return ρc*exp(c1*θ^(2/6) + c2*θ^(4/6) + c3*θ^(8/6) + c4*θ^(18/6) + c5*θ^(37/6) + c6*θ^(71/6))
end

function Psat_IAWPS_rho(ρ,Tsat_function)
    return Psat_IAWPS(Tsat_function(ρ))
end

function Tsat_IAWPS_rho_findzero(ρ)
    return find_zero(T -> ρsat_IAWPS(T) - ρ, 100 + 273.15)
end

function ρprime(T)
    θ = 1 - T/Tc
    return -ρsat_IAWPS(T)/Tc * (1/3*c1*θ^(-2/3) + 2/3*c2*θ^(-1/3) + 8/6*c3*θ^(1/3) + 18/16*c4*θ^2 + 37/6*c5*θ^(31/6) + 71/6*c6*θ^(65/6))
end
function Tsat_IAWPS_rho_newton(ρ)
    x = 100 + 273.15
    err = 1
    while err > 0.1
        xnew = x - (ρsat_IAWPS(x) - ρ)/ρprime(x)
        err = abs(xnew - x)
        x = xnew
    end
    return x
end


rhos = 0.1:0.5:200.0

Tsat1 = [PropsSI("T", "D", rhoi, "Q", 1, "Water") - 273.15 for rhoi in rhos]
# Tsat2 = [Tsat_IAWPS_rho_findzero(rhoi) - 273.15 for rhoi in rhos]
# Tsat3 = [Tsat_IAWPS_rho_newton(rhoi) - 273.15 for rhoi in rhos]

Psat1 = [PropsSI("P", "D", rhoi, "Q", 1, "Water")/1e5 for rhoi in rhos]
# Psat2 = [Psat_IAWPS_rho(rhoi, Tsat_IAWPS_rho_findzero)/1e5 for rhoi in rhos]
# Psat3 = [Psat_IAWPS_rho(rhoi, Tsat_IAWPS_rho_newton)/1e5 for rhoi in rhos]

# p1 = plot(rhos, [Tsat1],#, Tsat3],
#           ylabel="Tsat", xlabel="ρ")
# p2 = plot(rhos, [Psat1],#, Psat3],
#           ylabel="Psat", xlabel="ρ")
# plot(p1, p2, layout=(2,1))

# plot(rhos,Psat1)
# R = 8.314/(18/1000)
# plot!(rhos,rhos*R*(160+273.15)/1e5)

vmid = 5.0
vinf = 0.0001:0.001:0.01
vsup = 0.1:0.5:10.0

ps = 1000:10:1e6

Vsat_gas = [PropsSI("V", "P", pi, "Q", 1, "Water") for pi in ps]
Vsat_liq = [PropsSI("V", "P", pi, "Q", 0, "Water") for pi in ps]

plot(Vsat_liq, ps)
plot!(Vsat_gas, ps)
plot!(R)