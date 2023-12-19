using Plots
using SparseArrays

L  = 1.0
tf = 0.1

dx = 0.05
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
D = 1/dx * D

# Centered
Dprime = spdiagm(
    -1 => -1*ones(n),
     1 =>    ones(n)
   )
Dprime[1,1] = -1
Dprime[end,end] = 1
Dprime = 0.5/dx * Dprime


γ = 1.4

f(ρ,m,e) = m^2/ρ + e/Molar*(γ-1)
g(ρ,m,e) = e*m/ρ

f_vec(ρ,m,e) = [f(ρ[k], m[k], e[k]) for k in 1:n+1]
g_vec(ρ,m,e) = [g(ρ[k], m[k], e[k]) for k in 1:n+1]

# T0 = 110 + 273.15
# p0 = 1e5
# Molar = 18/1000
# R_gas = 8.314
# Cv = 1.4

# ρ0 = p0*Molar/(R_gas*T0)
# v0 = 0.001
# U0 = Cv*T0

# R₀ = 1/100
# Cd = 0.00

ρ0 = 1
v0 = 1
U0 = 1
R₀ = 1
Cd = 0
Molar = 1

ρ = ρ0*ones(n+1)
m = zeros(n+1)
e = ρ0*U0*ones(n+1)

ρ_save = []
m_save = []
e_save = []

for j in 1:nt
    global ρ,m,e

    v = m./ρ
    F = Cd/R₀*[abs(m[k])*m[k]/ρ[k] for k in 1:n+1]
    G = e/Molar*(γ-1).*Dprime*v

    ρ_new = ρ - dt*D*m    
    m_new = m - dt*(F + D*f_vec(ρ,m,e))
    e_new = e - dt*(G + D*g_vec(ρ,m,e))

    ρ = ρ_new
    m = m_new
    e = e_new

    ρ[1] = ρ0
    m[1] = ρ0*v0
    e[1] = ρ0*U0

    push!(ρ_save, ρ)
    push!(m_save, m)
    push!(e_save, e)
end

# n_save = length(m_save)
# v_save = [m_save[i]./ρ_save[i] for i in 1:n_save]
# T_save = [(E_save[i]./ρ_save[i] - 0.5*v_save[i].^2)/Cv for i in 1:n_save]
# p_save = [ρ_save[i].*T_save[i] * R_gas/Molar/1e5 for i in 1:n_save]


for j in 1:10:nt
    p1 = plot(x, ρ_save[j], title="ρ")
    p2 = plot(x, m_save[j], title="m")
    p3 = plot(x, e_save[j], title="e", 
                 xlabel="t = $(t[j])")

    p = plot(p1,p2,p3, layout=(3,1))
    display(p)
end