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
v0 = 400

ρ  = zeros(nt,N+1)
v  = zeros(nt,N+1)
ρ[1,:] = ρ0*ones(N+1) + dρ*exp.(-(x.-L/2).^2)
v[1,:] = v0*ones(N+1)
m = ρ.*v

U = zeros(2,N+1)
U[1,:] = ρ[1,:]
U[2,:] = m[1,:]

f(U) = [U[2], U[2]^2/U[1] + a*U[1]]

# Upwind (explicit)
for n in 1:nt-1
    U[:,1] = U[:,1] - dt/dx*(f(U[:,1]) - f(U[:,end]))
    for j in 2:N+1
        U[:,j] = U[:,j] - dt/dx*(f(U[:,j]) - f(U[:,j-1]))
    end
    p1 = plot(x, U[1,:], title="ρ", legend=false)
    p2 = plot(x, U[2,:], title="Vz", legend=false)
    p = plot(p1,p2,layout=(2,1))
    display(p)
end

# for n in 1:nt
#     p1 = plot(title="ρ", legend=false)
#     # plot!(x,ρ[n,:])
#     p2 = plot(title="Vz", legend=false)
#     # plot!(x,v[n,:])
#     p = plot(p1,p2,layout=(2,1))
#     display(p)
# end