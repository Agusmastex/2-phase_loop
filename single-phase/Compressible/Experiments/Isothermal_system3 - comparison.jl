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

T  = 0.04
L  = 10
dx = 0.05
dt = 0.00005
x  = 0:dx:L
t  = 0:dt:T
N  = length(x) - 1
nt = length(t)

R = 8.314
T0 = 298
Molar = 18/1000
a = R*T0/Molar
p0 = 10^5
ρ_baseline = p0*Molar/(R*T)
dρ = 0.01*ρ_baseline
v_uniform = 2.0

# Initial density distribution
ρ0_fun(x) = ρ_baseline + dρ*exp(-(x-L/2)^2)

# Initial velocity distribution
v0_fun(x) = v_uniform

# To compute the solution
# First index: space, Second index: (density, mass flux)
U = zeros(N+1,2)
U[:,1] = ρ0_fun.(x) # Discretized initial density
U[:,2] = ρ0_fun.(x).*v0_fun.(x) # Discretized initial velocity
U_init = copy(U)

# To store the solution
# First index: time, second index: space
ρ = zeros(nt,N+1)
v = zeros(nt,N+1)

f(U) = [U[2], U[2]^2/U[1] + a*U[1]]
D_upwind  = circulant([1 -1 zeros(1,N-1)  ])
D_center1 = circulant([0 -1 zeros(1,N-2) 1])

# Upwind 
for n in 1:nt-1
    global U
    F = mapslices(f,U,dims=2)
    U = U - dt/dx*D_upwind*F
    ρ[n+1,:] = U[:,1]
    v[n+1,:] = U[:,2]./U[:,1]
end

ρ_upwind = copy(ρ)
v_upwind = copy(v)
U = copy(U_init)

# FTCS
for n in 1:nt-1
    global U
    F = mapslices(f,U,dims=2)
    U = U - 0.5*dt/dx*D_center1*F
    ρ[n+1,:] = U[:,1]
    v[n+1,:] = U[:,2]./U[:,1]
end

ρ_FTCS = copy(ρ)
v_FTCS = copy(v)
U = copy(U_init)

slack = 2.5

for n in 1:nt
    p1 = plot(ylims=(ρ_baseline-dρ,ρ_baseline+dρ))
    plot!(x,ρ_upwind[n,:])
    plot!(x,ρ_FTCS[n,:])
    p2 = plot(ylims=(v_uniform-slack, v_uniform+slack))
    plot!(x,v_upwind[n,:])
    plot!(x,v_FTCS[n,:])
    p = plot(p1,p2,layout=(2,1))#,title="t=$(t[n])")
    display(p)
end