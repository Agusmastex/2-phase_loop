using Plots

function circulant(a)
    n = length(a)
    M = zeros(n,n)
    for j in 1:n
        for i in 1:n
            if i >= j
                M[i,j] = a[i-j+1]
            else
                M[i,j] = a[n+i-j+1]
            end
        end
    end
    return M
end

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
v0 = 600

ρ  = zeros(nt,N+1)
v  = zeros(nt,N+1)
ρ[1,:] = ρ0*ones(N+1) + dρ*exp.(-(x.-L/2).^2)
v[1,:] = v0*ones(N+1)
m  = ρ.*v
ρ_init = copy(ρ)
m_init = copy(m)

f(ρ,m) = m^2/ρ + a*ρ

D_upwind = circulant([1 -1 zeros(1,N-1)])

# Upwind (explicit)
for n in 1:nt-1
    ρ[n+1,:] = ρ[n,:] - c*D_upwind*m[n,:]
    F = [f(ρ[n,j], m[n,j]) for j in 1:N+1]
    m[n+1,:] = m[n,:] - c*D_upwind*F
end
ρ_upwind = copy(ρ)
v_upwind = m./ρ
ρ = copy(ρ_init)
m = copy(m_init)




for n in 1:nt
    p1 = plot(title="ρ", legend=false)
    plot!(x,ρ_upwind_e[n,:])
    plot!(x,ρ_upwind_i[n,:])
    p2 = plot(title="Vz", legend=false)
    plot!(x,v_upwind_e[n,:])
    plot!(x,v_upwind_i[n,:])
    p = plot(p1,p2,layout=(2,1))
    display(p)
end