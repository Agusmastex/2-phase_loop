using Plots
using SparseArrays
using LinearAlgebra


L  = 20.0
dt = 0.00001
dx = 0.05
x  = 0:dx:L
t  = 0:dt:T
nt = length(t)
N  = length(x) - 1


R₀ = 1/100
Cd = 2.5
R_gas = 8.314
T0 = 298
Molar = 18/1000
a = R_gas*T0/Molar
p0 = 1e5
ρ_uniform = p0*Molar/(R_gas*T0)
v_uniform = 300.0

ρ0 = ρ_uniform .+ exp.(-(x .- L/2).^2)

Q_sol = zeros(nt,2N+2)
Q_sol[1,1:N+1]    .= ρ_uniform
Q_sol[1,N+2:2N+2] .= ρ_uniform*v_uniform

Q_sol[1,1:N+1]    = ρ0
Q_sol[1,N+2:2N+2] = ρ0*v_uniform

function K(Q)
    ρ = Q[1:N+1]
    m = Q[N+2:end]
    k = Cd/R₀*[m[i]^2/ρ[i] for i in 1:N+1]
    return [zeros(N+1); k]
end

function g(Q)
    ρ = Q[1:N+1]
    m = Q[N+2:end]
    f = [m[i]^2/ρ[i] + a*ρ[i] for i in 1:N+1]
    return [m; f]
end

D = spdiagm(
    -1 => -1*ones(N),
     0 =>    ones(N+1)
     )
D[1,1] = 0
D = dropzeros(D)

big_D = [D zeros(N+1,N+1);
         zeros(N+1,N+1) D]


F(Qn,Qnp1) = 1/dt*(Qnp1 - Qn) + 1/dx*big_D*g(Qnp1) + K(Qnp1)

J1 = 1/dt*I
J2 = 1/dx*Tridiagonal(-ones(N), ones(N+1), zeros(N))
function Jac(Q)
    ρ = Q[1:N+1]
    m = Q[N+2:2N+2]
    diag3  = [1/dx*(-m[i]^2/ρ[i] + a) - Cd/R₀*(m[i]/ρ[i])^2 for i in 1:N+1]
    ldiag3 = [1/dx*((m[i]/ρ[i])^2 - a) for i in 1:N]
    J3 = Tridiagonal(ldiag3, diag3, zeros(N))
    diag4  = [1/dt + 1/dx*(2*m[i]/ρ[i]) + 2*Cd/R₀*(m[i]/ρ[i]) for i in 1:N+1]
    ldiag4 = [-1/dx*(2*m[i]/ρ[i]) for i in 1:N]
    J4 = Tridiagonal(ldiag4, diag4, zeros(N))
    return sparse([J1 J2; J3 J4])
end

ρ = zeros(nt,N+1)
v = zeros(nt,N+1)

tol = 1e-6
for n in 1:nt-1
    println()
    println(n)
    res = 1
    Qk = Q_sol[n,:]
    k = 0
    while res > tol
        k = k + 1
        print(k)
        A = Jac(Qk)
        b = -F(Q_sol[n,:],Qk)
        xk = A\b
        res = norm(xk)
        Qk = Qk + xk
    end
    Q_sol[n+1,:] = Qk
    Q_sol[n+1,1] = ρ_uniform
    Q_sol[n+1,N+2] = ρ_uniform*v_uniform
    ρ[n+1,:] = Qk[1:N+1]
    v[n+1,:] = Qk[N+2:2N+2]./ρ[n+1,:]
    # p1 = plot(x,ρ[n+1,:])
    # p2 = plot(x,v[n+1,:])
    # p  = plot(p1,p2, layout=(2,1),title="t=$(t[n])")
    # display(p)
end

for n in 1:nt
    p1 = plot(x,ρ[n,:])
    p2 = plot(x,v[n,:])
    p  = plot(p1,p2, layout=(2,1),title="t=$(t[n])")
    display(p)
end