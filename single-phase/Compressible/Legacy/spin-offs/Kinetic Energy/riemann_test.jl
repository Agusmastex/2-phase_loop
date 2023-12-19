using Plots
using SparseArrays

L  = 1.0
tf = 0.9

dx = 0.01
dt = 0.005

x  = 0:dx:L
n  = length(x) - 1
t  = 0:dt:tf
nt = length(t) - 1

# Upwind
D = spdiagm(
    -1 => -1*ones(n),
     0 =>    ones(n+1)
   )
D = 1/dx * dropzeros(D)

γ = 1.4

f(ρ,m,E) = m^2/ρ*(2-0.5*γ) + (γ-1)*E
g(ρ,m,E) = γ*E - 0.5*m^2/ρ*(γ-1)

f_vec(ρ,m,E) = [f(ρ[k], m[k], E[k]) for k in 1:n+1]
g_vec(ρ,m,E) = [g(ρ[k], m[k], E[k]) for k in 1:n+1]

T0 = 110 + 273.15
p0 = 1e5
Molar = 18/1000
R_gas = 8.314
Cv = 1.4e3

ρ0 = p0*Molar/(R_gas*T0)
v0 = 1.0
U0 = Cv*T0

R₀ = 1/100
Cd = 0.01

ρ = ρ0*ones(n+1) 
m = [x[i] < L/8 ? ρ0*v0 : 0 for i in 1:n+1]
E = [x[i] < L/8 ? ρ0*U0 + 0.5*ρ0*v0^2 : ρ0*U0 for i in 1:n+1]


# ρ = ρ0*ones(n+1)
# m = ρ0*(v0-0.001)*ones(n+1)
# E = (ρ0*U0 + 0.5*ρ0*(v0-0.001)^2)*ones(n+1)

ρ_save = [ρ]
m_save = [m]
E_save = [E]

for j in 1:nt
    global ρ,m,E
    F = 2/R₀*Cd*[abs(m[k])*m[k]/ρ[k] for k in 1:n+1]

    ρ_new = ρ - dt*D*m    
    m_new = m - dt*(F + D*f_vec(ρ,m,E))
    E_new = E - dt*D*g_vec(ρ,m,E)

    ρ = ρ_new
    m = m_new
    E = E_new

    ρ[1] = ρ0
    m[1] = ρ0*v0
    E[1] = ρ0*U0 + 0.5*ρ0*v0^2

    push!(ρ_save, ρ)
    push!(m_save, m)
    push!(E_save, E)
end

n_save = length(m_save)
v_save = [m_save[i]./ρ_save[i] for i in 1:n_save]
T_save = [(E_save[i]./ρ_save[i] - 0.5*v_save[i].^2)/Cv for i in 1:n_save]
p_save = [ρ_save[i].*T_save[i] * R_gas/Molar/1e5 for i in 1:n_save]


for j in 1:nt
    p1 = plot(x, v_save[j], title="v")
    p2 = plot(x, p_save[j], title="p")
    # p3 = plot(x, E_save[j], title="E",
    p3 = plot(x, T_save[j] .- 273.15, title="T", 
                 xlabel="t = $(t[j])")

    p = plot(p1,p2,p3, layout=(3,1))
    display(p)
end