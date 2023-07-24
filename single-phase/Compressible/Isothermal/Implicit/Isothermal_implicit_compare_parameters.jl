using Plots
using SparseArrays
using LinearAlgebra

T  = 5.0
L  = 5.0
dt = 1e-2
dx = 0.2
x_plot  = 0:dx:L
nt = length(0:dt:T)
N  = length(x_plot) - 1

R₀ = 1/100
R_gas = 8.314
T0 = 298
Molar = 18/1000
a = R_gas*T0/Molar
p0 = 1e5

F(Qk,Qn) = 1/dt*(Qk - Qn) + 1/dx*big_D*A(Qk) + K(Qk)

# Difference matrix
D = spdiagm(
    -1 => -1*ones(N),
     0 =>    ones(N+1)
   )
D[1,1] = 0
dropzeros!(D)
big_D = sparse([D zeros(N+1,N+1); zeros(N+1,N+1) D])

# Initialize
ρ_uniform = p0*Molar/(R_gas*T0)
v_uniform = 2.0
Q1 = zeros(2N+2)
Q1[1:N+1]    .= ρ_uniform
Q1[N+2:2N+2] .= ρ_uniform*v_uniform

rhos = []
vs = []

# Tweak parameters
Δn = round(Int, (nt-1)/500)
Δn = 1
# a = 1.0
as = (0, 1.0, 2.0)
Cd = 0.0001
# Cds = (0.0001, 0.0002, 0.0003)

simulation_number = 0 
for a in as
    simulation_number += 1
    ρ_save = []
    v_save = []


    function K(Q)
        ρ = Q[1:N+1]
        m = Q[N+2:2N+2]
        k = [Cd/R₀*m[j]^2/ρ[j] for j in 2:N+1]
        return [zeros(N+2); k]
    end

    f(ρ,m) = m^2/ρ + a*ρ

    function A(Q)
        ρ = Q[1:N+1]
        m = Q[N+2:2N+2]
        fvec = [f(ρ[j], m[j]) for j in 1:N+1]
        return [m; fvec]
    end

    # Jacobian
    J1 = 1/dt*sparse(I,N+1,N+1)
    J2 = 1/dx*spdiagm(
        -1 => -1*ones(N),
         0 =>    ones(N+1)
       )
    J2[1,1] = 0
    function Jac(Q)
        ρ = Q[1:N+1]
        m = Q[N+2:2N+2]
        diagonal_3 =       [       1/dx*( -(m[j]/ρ[j])^2 + a)  -   Cd/R₀*(m[j]/ρ[j])^2 for j in 1:N+1]
        lower_diagonal_3 = [       1/dx*(  (m[j]/ρ[j])^2 - a)                          for j in 2:N+1]
        diagonal_4 =       [1/dt + 1/dx*( 2*m[j]/ρ[j])         + 2*Cd/R₀* m[j]/ρ[j]    for j in 1:N+1]
        lower_diagonal_4 = [       1/dx*(-2*m[j]/ρ[j])                                 for j in 2:N+1]

        J3 = 1/dx*spdiagm(
        -1 => lower_diagonal_3,
         0 => diagonal_3
       )
        J4 = 1/dx*spdiagm(
        -1 => lower_diagonal_4,
         0 => diagonal_4
       )

        J3[1,1] = 0
        J4[1,1] = 1/dt

        return sparse([J1 J2; J3 J4])
    end

    # Solve
    tol = 1e-6
        Qn = copy(Q1)
        Qk = copy(Q1)
        for n in 1:nt-1
            res = 1.0
            k = 0
            while res > tol
                k = k + 1
                x = Jac(Qk)\-F(Qk,Qn)
                Qk = Qk + x
                res = norm(x)
                print(k)
            end
            println()
            if n % Δn == 0
                ρ = Qk[1:N+1]
                m = Qk[N+2:2N+2]
                push!(ρ_save, ρ)
                push!(v_save, m./ρ)
                println()
                println("$simulation_number $(round(n/nt*100, digits=1))% ")
            end
            Qn = copy(Qk)
        end

        ρ_save = stack(ρ_save, dims=1)
        v_save = stack(v_save, dims=1)

        push!(rhos, ρ_save)
        push!(vs, v_save)
end


n_save = length(rhos[1][:,1])
t_save = LinRange(0,T,n_save)
number_graphs = length(rhos)
println("Memory $((sizeof(rhos) + sizeof(vs))/1000) MB")

for n in 1:n_save
    p1 = plot(title="ρ", xlabel=("t= $(round(t_save[n], digits=1))"))#, ylims=(0.72, 0.80))
    p2 = plot(title="Vz")#, ylims=(1.80,2))
    for i in 1:number_graphs
    plot!(p1, x_plot, rhos[i][n,:], label="a=$(as[i])")
    plot!(p2, x_plot, vs[i][n,:], label="a=$(as[i])")
    end
    p  = plot(p1,p2, layout=(2,1))
    display(p)
end