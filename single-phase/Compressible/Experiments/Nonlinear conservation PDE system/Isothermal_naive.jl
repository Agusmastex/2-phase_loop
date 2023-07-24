using Plots

T  = 0.03
L  = 10
dx = 0.05
dt = 0.00005
x  = 0:dx:L
t  = 0:dt:T
N  = length(x) - 1
nt = length(t)
c = dt/dx

R = 8.314
T = 298
Molar = 18/1000
a = R*T/Molar
p0 = 10^5
ρ0 = p0*Molar/(R*T)
dρ = 0.01*ρ0
v0 = 0

ρ  = zeros(nt,N+1)
v  = zeros(nt,N+1)
ρ[1,:] = ρ0*ones(N+1) + dρ*exp.(-(x.-L/2).^2)
ρ[2,:] = ρ[1,:] # Needed for Leapfrog
v[1,:] = v0*ones(N+1)
v[2,:] = v[1,:] # Needed for Leapfrog
m  = ρ.*v
ρ_init = copy(ρ)
m_init = copy(m)

F(ρ,m) = m^2/ρ + a*ρ

# Upwind
for n in 1:nt-1
        ρ[n+1,1] = ρ[n,1] - c*(m[n,1] - m[n,end])
        m[n+1,1] = m[n,1] - c*(F(ρ[n,1], m[n,1]) - F(ρ[n,end], m[n,end]))
    for j in 2:N+1
        ρ[n+1,j] = ρ[n,j] - c*(m[n,j] - m[n,j-1])
        m[n+1,j] = m[n,j] - c*(F(ρ[n,j], m[n,j]) - F(ρ[n,j-1], m[n,j-1]))
    end
end
ρ_upwind = copy(ρ)
v_upwind = m./ρ
ρ = copy(ρ_init)
m = copy(m_init)

# Lax-Wendroff
# for n in 1:nt-1
#         ρ[n+1,1] = ρ[n,1] - 0.5*c*(m[n,2] - m[n,end]) + 0.5*c^2*(m[n,2] - 2*m[n,1] + m[n,end])
#         m[n+1,1] = m[n,1] - 0.5*c*(F(ρ[n,2], m[n,2]) - F(ρ[n,end], m[n,end])) + 
#                           + 0.5*c^2(F(ρ[n,2],m[n,2]) - 2*F(ρ[n,1],m[n,1]) + F(ρ[n,end],m[n,end]))
#     for j in 2:N
#         ρ[n+1,j] = ρ[n,j] - 0.5*c*(m[n,j+1] - m[n,j-1]) + 0.5*c^2*(m[n,j+1] - 2*m[n,j] + m[n,j-1])
#         m[n+1,j] = m[n,j] - 0.5*c*(F(ρ[n,j+1], m[n,j+1]) - F(ρ[n,j-1], m[n,j-1])) + 
#                           + 0.5*c^2(F(ρ[n,j+1],m[n,j+1]) - 2*F(ρ[n,j],m[n,j]) + F(ρ[n,j-1],m[n,j-1]))
#     end
#         # ρ[n+1,end] = ρ[n,end] - 0.5*c*(m[n,1] - m[n,end-1]) + 0.5*c^2*(m[n,1] - 2*m[n,end] + m[n,end-1])
#         # m[n+1,end] = m[n,end] - 0.5*c*(F(ρ[n,1], m[n,1]) - F(ρ[n,end-1], m[n,end-1])) + 
#         #                   + 0.5*c^2(F(ρ[n,1],m[n,1]) - 2*F(ρ[n,end],m[n,end]) + F(ρ[n,end-1],m[n,end-1]))
# end
# ρ_wendroff = copy(ρ)
# v_wendroff = m./ρ
# ρ = copy(ρ_init)
# m = copy(m_init)

# FTCS
for n in 1:nt-1
        ρ[n+1,1] = ρ[n,1] - 0.5*c*(m[n,2] - m[n,end])
        m[n+1,1] = m[n,1] - 0.5*c*(F(ρ[n,2], m[n,2]) - F(ρ[n,end], m[n,end]))
    for j in 2:N
        ρ[n+1,j] = ρ[n,j] - 0.5*c*(m[n,j+1] - m[n,j-1])
        m[n+1,j] = m[n,j] - 0.5*c*(F(ρ[n,j+1], m[n,j+1]) - F(ρ[n,j-1], m[n,j-1]))
    end
        ρ[n+1,end] = ρ[n,end] - 0.5*c*(m[n,1] - m[n,end-1])
        m[n+1,end] = m[n,end] - 0.5*c*(F(ρ[n,1], m[n,1]) - F(ρ[n,end-1], m[n,end-1]))
end
ρ_FTCS = copy(ρ)
v_FTCS = m./ρ
ρ = copy(ρ_init)
m = copy(m_init)

# Leapfrog
for n in 2:nt-1
        ρ[n+1,1] = ρ[n-1,1] - c*(m[n,2] - m[n,end])
        m[n+1,1] = m[n-1,1] - c*(F(ρ[n,2], m[n,2]) - F(ρ[n,end], m[n,end]))
    for j in 2:N      
        ρ[n+1,j] = ρ[n-1,j] - c*(m[n,j+1] - m[n,j-1])
        m[n+1,j] = m[n-1,j] - c*(F(ρ[n,j+1], m[n,j+1]) - F(ρ[n,j-1], m[n,j-1]))
    end
        ρ[n+1,end] = ρ[n-1,end] - c*(m[n,1] - m[n,end-1])
        m[n+1,end] = m[n-1,end] - c*(F(ρ[n,1], m[n,1]) - F(ρ[n,end-1], m[n,end-1]))
end
ρ_leapfrog = copy(ρ)
v_leapfrog = m./ρ
ρ = copy(ρ_init)
m = copy(m_init)

for n in 1:nt
    p1 = plot(title="ρ", legend=false)
    # plot!(x,ρ_upwind[n,:])
    # plot!(x,ρ_wendroff[n,:])
    # plot!(x,ρ_FTCS[n,:])
    plot!(x,ρ_leapfrog[n,:])
    p2 = plot(title="Vz", legend=false)
    # plot!(x,v_upwind[n,:])
    # plot!(x,v_wendroff[n,:])
    # plot!(x,v_FTCS[n,:])
    plot!(x,v_leapfrog[n,:])
    p = plot(p1,p2,layout=(2,1))
    display(p)
end