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

n = 500
dz = L/n
z  = 0:dz:L
T  = T0*ones(n + 1)
Q_vec = Q_tilde*[zi < Lh ? 1 : 0 for zi in z]
dt = 0.01
tf = 1.5*L/Vz
t  = 0:dt:tf

# D = 1/dz*Tridiagonal(-ones(n), ones(n + 1), zeros(n))   # Diferencia hacia atrás, upwind
D = 1/dz*Tridiagonal(-ones(n), zeros(n + 1), ones(n)) # Diferencia centrada, inestable

Tmax = T0 + Q/(Vz*ρ0*Cp*A)	

A = I + dt*Vz*D
A[1,:]    = [1; zeros(n)]
# A[end, :] = [zeros(n-1);-1;1]

for ti in 0:dt:tf
    global T

    b = T + dt*Q_vec
    b[1] = T0
    # b[end] = 0
    T = A\b

    # plot(z,T, ylim = (T0-0.1,Tmax+0.1),  title="t = $ti", show=true)
    plot(z, T, show=true)
end

