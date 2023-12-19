using Plots
using SparseArrays

T  = 1.0
L  = 5.0
dt = 0.000001
dx = 0.2
x  = 0:dx:L
t  = 0:dt:T
nt = length(t)
N  = length(x) - 1

# Upwind
D = spdiagm(
    -1 => -1*ones(N),
     0 =>    ones(N+1)
     )
D[1,1] = 0
D = 1/dx * dropzeros(D)

R₀ = 1/100
Cd = 0.0001
R_gas = 8.314
T0 = 298
Molar = 18/1000
a = R_gas*T0/Molar
a = 10.0
p0 = 1e5
ρ_uniform = p0*Molar/(R_gas*T0)
v_uniform = 2.0

ρ = zeros(nt, N+1)
m = zeros(nt, N+1)
ρ[1,:] .= ρ_uniform
m[1,:] .= ρ_uniform*v_uniform

ρ_inflow = ρ_uniform
v_inflow = v_uniform

ρ[1,1] = ρ_inflow
m[1,1] = ρ_inflow*v_inflow

f(ρ,m) = m^2/ρ + a*ρ

for n in 1:nt-1
    F = [f(ρ[n,j],m[n,j]) for j in 1:N+1]
    J = Cd/R₀*[m[n,j]*abs(m[n,j])/ρ[n,j] for j in 1:N+1]
    J[1] = 0

    ρ[n+1,:] = ρ[n,:] - dt*D*m[n,:]
    m[n+1,:] = m[n,:] - dt*(D*F + J)
end

v = m./ρ

for n in 1:100:nt
    p1 = plot(x,ρ[n,:])#, ylims = (0.999*ρ_uniform, 1.001*ρ_inflow))
    p2 = plot(x,v[n,:])#, ylims = (0, 1.2*v_inflow))
    p  = plot(p1,p2, layout=(2,1),title="t=$(t[n])")
    display(p)
end