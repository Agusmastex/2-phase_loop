using Plots
using LinearAlgebra
using DifferentialEquations

Q  = 1000.0       # J/s
ρ0 = 1000.0       # kg/m^3
Cp = 4184.0       # J/kg K
R  = 1.0/1000     # m
A  = π*R^2        # m^2
Vz = 2.0          # m/s
Lh = 0.5          # m
L  = 3.0         # m                

T0 = 25.0
Q_tilde   = Q/(ρ0*Cp*A*Lh)

n = 50
dz = L/n
z  = 0:dz:L
T  = T0*ones(n + 1) 
Q_vec = Q_tilde*[zi < Lh ? 1 : 0 for zi in z]
dt = 0.005
tf = 1.5*L/Vz
t  = 0:dt:tf

# Matriz de diferencia finita hacia atrás
D = 1/dz*Tridiagonal(-ones(n), ones(n + 1), zeros(n))

# Force boundary conditions
D[1,1] = 0
Q_vec[1] = 0

# For plotting
Tmax = T0 + Q/(Vz*ρ0*Cp*A)	

# Solve ODE 
fun(T,p,t) = Q_vec - Vz*D*T
prob = ODEProblem(fun, T, [0,tf])
sol  = solve(prob, saveat=t)
plot(sol)

# Plot
t = sol.t
for j in 1:length(t)
    Tj = sol[:,j]
    plot(z,Tj, ylim = (T0-0.1,Tmax+0.1), marker=:circle, title="t = $(t[j])", show=true)
end
