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

n = 500           # Número de intervalos (existen n + 1 nodos)
dz = L/n
z  = 0:dz:L
Q_on  = [zi < Lh ? 1 : 0 for zi in z]
Q_vec = Q_tilde*Q_on
dt = 0.01
tf = 1.5*L/Vz
t  = 0:dt:tf
nt = length(t)     # Número de nodos temporales (existen nt - 1 intervalos de tiempo)
T  = zeros(nt,n+1)
T[1,:] = T_initial*ones(n+1)

T_inside = @view T[2:end-1]
# Q_on = Q_on[2:end-1]

# D = 0.5/dz*Tridiagonal(-ones(n - 2), ones(n - 1), zeros(n - 2))
D = 1/dz*Tridiagonal(-ones(n), ones(n + 1), zeros(n))
D = copy(Matrix(D))
D[1,1] = -dz/dt
# D[end,:] = D[end-1,:]

Tmax = T_initial + Q/(Vz*ρ0*Cp*A)	

for j in 1:nt-1
global T
T[j+1,:] = T[j,:] + dt*(Q_vec - Vz*D*T[j,:])
T[j+1,1] = T_initial
end

for j in 1:1:nt
plot(z,T[j,:], ylim = (T_initial-0.1,Tmax+0.1),  title="t = $(t[j])", show=true)
end

