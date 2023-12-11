using Plots
using CoolProp
using Roots
using CSV

Molar = 18/1000
R = 8.314
A = 8.14019
B = 1810.94
C = 244.485

function Tsat_antoine(ρ)
    f(Tsat) = ρ*R*(Tsat)/Molar - 101325/760*10^(A - B/(C + Tsat - 273.15))
    fprime(Tsat) = ρ*R/Molar - 101325/760*10^(A - B/(C+Tsat))*log(10)*(B/(C+Tsat)^2)
    err = 1
    x = 100 
    while err > 0.1
        xnew = x - f(x)/fprime(x)
        err = xnew - x
        x = xnew
    end
    return x
end

function Psat_ideal(ρ, Tsat_function)
    return ρ*R*(Tsat_function(ρ))/Molar
end

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

function Psat_IAWPS_rho(ρ)
    return Psat_IAWPS(Tsat_IAWPS_rho(ρ))
end

function Tsat_IAWPS_rho(ρ)
    return find_zero(T -> ρsat_IAWPS(T) - ρ, 100 + 273.15)
end

rhos = 0.5:0.01:2.5
Psat1 = [PropsSI("P", "D", rhoi, "Q", 1, "Water")/1e5 for rhoi in rhos]
# Psat3 = [Psat_ideal(rhoi, define tsat)/1e5 for rhoi in rhos]
@time Psat3 = [Psat_IAWPS_rho(rhoi)/1e5 for rhoi in rhos]

Tsat1 = [PropsSI("T", "D", rhoi, "Q", 1, "Water") - 273.15 for rhoi in rhos]
# Tsat2 = [Tsat_antoine(rhoi) for rhoi in rhos]
@time Tsat3 = [Tsat_IAWPS_rho(rhoi) - 273.15 for rhoi in rhos]

p1 = plot(rhos, [Tsat1, Tsat3],
          ylabel="Tsat", xlabel="ρ")
p2 = plot(rhos, [Psat1, Psat3],
          ylabel="Psat", xlabel="ρ")
plot(p1, p2, layout=(2,1))

## IAWPS Psat verification
# Trange = 1:0.01:99
# @time Psats1 = [PropsSI("P", "T", Ti + 273.15, "Q", 0, "Water")/1e5 for Ti in Trange]
# @time Psats2 = Psat_IAWPS.(Trange)/1e5
# plot(Trange, [Psats1, Psats2])
## (verified)


# IAWPS rhosat verification
Trange = (99 + 273.15):1:(150 + 273.15) 
# @time rho_sats1 = [PropsSI("D", "T", Ti, "Q", 1, "Water") for Ti in Trange]
# @time rho_sats2 = ρsat_IAWPS.(Trange)
# plot(Trange .- 273.15, [rho_sats1, rho_sats2])
## (verified)