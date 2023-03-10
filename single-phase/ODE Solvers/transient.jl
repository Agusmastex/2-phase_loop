using Plots
using LinearAlgebra
using DifferentialEquations
using Interpolations

Q  = 1000.0       # J/s
ρ0 = 1000.0       # kg/m^3
Cp = 4184.0       # J/kg K
R  = 1.0/1000     # m
A  = π*R^2        # m^2
Vz = 2.0          # m/s
Lh = 0.5          # m
L  = 3.0          # m                

f  = 0.04
# β  = 210e-6
β  = 0.1
g  = 9.8

T0 = 25.0
Q_tilde   = Q/(ρ0*Cp*A*Lh)

n = 500
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

# Solve for temperature
Tfun(T,p,t) = Q_vec - Vz*D*T
Tprob = ODEProblem(Tfun, T, [0,tf])
Tsol = solve(Tprob, saveat=t)

t  = Tsol.t

Δp_loss   = 0.25*f*(L/R)*ρ0*Vz^2
Δp_static = ρ0*g*L
Δp_max    = Δp_loss + Δp_static

# Solve for pressure and simultaneously plot
for j in 1:length(t)
    Tj = Tsol[:,j]

    Tj_fun = LinearInterpolation(z, Tj)
    p_fun(u,p,z) = - 0.25*f*ρ0*Vz^2/R - ρ0*(1 - β*(Tj_fun(z) - T0))*g
    p_prob = ODEProblem(p_fun, 0, [0,L])
    psol = solve(p_prob, saveat=z)

    p_plot = psol.u
    z_plot = psol.t

    p_plot = p_plot + ρ0*g*z_plot             # Para graficar presión modifcada
    p_plot = p_plot/1000

    if j % 1 == 0
       p1 = plot(z, Tj,
                   ylim = (T0-0.1,Tmax+0.1),
                   title="t = $(t[j])")

       p2 = plot(z_plot, p_plot, 
                   # ylim = (-Δp_loss/1000, 0)      # Para graficar presión modificada
                   # ylim = (-Δp_max/1000, 0)     # Para graficar presión normal
       )

       p3 = plot(p1, p2, layout=(2,1))
       display(p3)
    end
end