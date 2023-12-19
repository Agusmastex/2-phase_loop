using Plots
using SparseArrays

T  = 20.0
L  = 30.0
dt = 1e-2
dx = 0.2
x_plot  = 0:dx:L
nt = length(0:dt:T)
N  = length(x_plot) - 1

# Upwind
D = spdiagm(
    -1 => -1*ones(N),
     0 =>    ones(N+1)
   )
D[1,1] = 0
D = 1/dx * dropzeros(D)

R₀ = 1/100
R_gas = 8.314
T0 = 100 + 273.15
Molar = 18/1000
a = R_gas*T0/Molar
p0 = 1e5
ρ_uniform = p0*Molar/(R_gas*T0)
v_uniform = 2.0

# Tweak paramteres
Δn = round(Int, (nt-1)/500)
Δn = 1
a = 0.0
as = (0.0, 1.0, 2.0)
# as = (0.0)
Cd = 0.0001
# Cds = (0.0001, 0.0002, 0.0003)

rhos = []
vs = []

simulation_number = 0
for a in as
    simulation_number += 1

    f(ρ,m) = m^2/ρ + a*ρ

    ρ_save = []
    m_save = []

    ρ = zeros(N+1)
    m = zeros(N+1)
    ρ .= ρ_uniform
    m .= ρ_uniform*v_uniform

    ρ[1] = ρ_uniform
    m[1] = ρ_uniform*v_uniform

    for n in 1:nt-1
        global ρ, m
        F = [f(ρ[j],m[j]) for j in 1:N+1]
        J = Cd/R₀*[m[j]*abs(m[j])/ρ[j] for j in 1:N+1]
        J[1] = 0

        ρ = ρ - dt*D*m
        m = m - dt*(D*F + J)

        if n % Δn == 0
            push!(ρ_save, ρ)
            push!(m_save, m)
            println("$simulation_number $(round(n/nt*100, digits=1))% ")
        end
    end

    print("End. ")

    ρ_save = stack(ρ_save, dims=1)
    m_save = stack(m_save, dims=1)
    v_save = m_save./ρ_save

    push!(rhos, ρ_save)
    push!(vs, v_save)

end

n_save = length(rhos[1][:,1])
t_save = LinRange(0,T,n_save)
number_graphs = length(rhos)
println("Memory $((sizeof(rhos) + sizeof(vs))/1000) MB")

@gif for n in 1:n_save
    p1 = plot(title="ρ", xlabel=("t= $(round(t_save[n], digits=1))"))#, ylims=(0.72, 0.80))
    p2 = plot(title="Vz")#, ylims=(1.80,2))
    for i in 1:number_graphs
    plot!(p1, x_plot, rhos[i][n,:], label="a=$(as[i])")
    plot!(p2, x_plot, vs[i][n,:], label="a=$(as[i])")
    end
    p  = plot(p1,p2, layout=(2,1))
    # display(p)
end