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

n_intervals = 50
dz = L/n_intervals
z = 0:dz:L
T = zeros(n_intervals + 1) .+ T_initial
Q_on = [zi < Lh ? 1 : 0 for zi in z]
dt = 0.01
tf = 7.0
t  = 0:dt:tf

T_inside = @view T[2:end-1]
# Q_on = Q_on[2:end-1]

# D = 0.5/dz*Tridiagonal(-ones(n_intervals - 2), ones(n_intervals - 1), zeros(n_intervals - 2))
D = 1/dz*Tridiagonal(-ones(n_intervals), ones(n_intervals + 1), zeros(n_intervals))
D = copy(Matrix(D))
# D[1,1] = 0
# D[end,:] = D[end-1,:]

Tmax = T_initial + Q/(Vz*ρ0*Cp*A)	

for ti in 0:dt:tf

global T_inside, T
# T_inside = T_inside + dt*(Q_on*Q_tilde - Vz*D*T_inside)
# T[end] = T[end-1]

T = T + dt*(Q_on*Q_tilde - Vz*D*T)
T[1] = T_initial
T[end] = T[end-1]

plot(z,T, ylim = (T_initial-0.1,Tmax+0.1), marker=:circle, title="t = $ti", show=true)
# plot(z, T, marker=:circle, show=true)
end

