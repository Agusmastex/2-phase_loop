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

T_initial = 25.0
Q_tilde   = Q/(ρ0*Cp*A*Lh)

n_intervals = 500
dz = L/n_intervals
z  = 0:dz:L
T  = T_initial*ones(n_intervals + 1)
Q_on  = [zi < Lh ? 1 : 0 for zi in z]
Q_vec = Q_tilde*Q_on
dt = 0.01
tf = 1.5*L/Vz
t  = 0:dt:tf

D = 1/dz*Tridiagonal(-ones(n_intervals), ones(n_intervals + 1), zeros(n_intervals))

Tmax = T_initial + Q/(Vz*ρ0*Cp*A)	

for ti in 0:dt:tf
global T

A = I + dt*Vz*D
b = T + dt*Q_vec
T = A\b
T[1] = T_initial

# plot(z,T, ylim = (T_initial-0.1,Tmax+0.1),  title="t = $ti", show=true)
plot(z,T,marker=:circle, show=true)
end

