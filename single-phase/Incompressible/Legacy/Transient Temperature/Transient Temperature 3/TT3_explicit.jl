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
Lh = 4.0          # m
L  = 10.0         # m                

T0 = 25.0
Q_tilde   = Q/(ρ0*Cp*A*Lh)

n = 50
dz = L/n
z  = 0:dz:L
T  = T0*ones(n + 1) 
Q_vec = Q_tilde*[zi < Lh ? 1 : 0 for zi in z]
dt = 0.01
tf = 1.5*L/Vz
t  = 0:dt:tf

D = 1/dz*Tridiagonal(-ones(n), ones(n + 1), zeros(n)) # Diferencia hacia atrás

Tmax = T0 + Q/(Vz*ρ0*Cp*A)	

@gif for ti in t

global T
T = T + dt*(Q_vec - Vz*D*T)
T[1] = T0

plot(z,T, ylim = (T0-0.1,Tmax+0.1), title="t = $ti")
end every 5

