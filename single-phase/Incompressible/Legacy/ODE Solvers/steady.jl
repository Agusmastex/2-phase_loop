using Plots
using ToeplitzMatrices
using LinearAlgebra
using DifferentialEquations

Q  = 1000.0       # J/s
ρ0 = 1000.0       # kg/m^3
Cp = 4184.0       # J/kg K
R  = 1.0/1000     # m
A  = π*R^2        # m^2
Vz = 2.0          # m/s
L  = 1.0          # m                
Lh = 0.6          # m

β  = 210e-6	      # K⁻¹ 
β  = 0.01	      # K⁻¹ 
T0 = 25.0		  # ºC
g  = 9.8		  # m/sꜝ
f  = 0.01

Q_tilde(z) = z < Lh ? Q/(ρ0*Cp*A*Lh) : 0

tf = 1.5*L/Vz
# Tmax = T0 + Q/(Vz*ρ0*Cp*A)	

# H[1] : Temperature
# H[2] : Pressure

fun(H,p,z) = [Q_tilde(z)/Vz,
			 -0.25*f*ρ0*Vz^2 - ρ0*(1 - β*(H[1]-T0))*g]

prob = ODEProblem(fun, [T0,0], [0,L])
dz = 0.01
z  = 0:dz:L
sol  = solve(prob, Rodas4())

T  = sol(z)[1,:]
p  = sol(z)[2,:]

p = p/1000
p = p .- p[end]

plot(sol.t, sol[1,:], marker=:circle)

plot(z, [T,p], layout=(2,1))