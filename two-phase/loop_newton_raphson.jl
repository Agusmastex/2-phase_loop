using LinearAlgebra
using SparseArrays
using CoolProp
using Plots

""" 
Solves the 2-phase HEM 3-equation loop problem 
Making the difference matrix cyclical fails
"""

@time begin

# Auxiliary functions 

  function plot_this(list, names)
    plots = [plot(s, item, title=name) for (item, name) in zip(list, names)]
    p = plot(plots...)
    return p
  end

# Geometry
    D = 9.1e-3
    H = 2.3275
    W = 2.02
    L = 2(H + W)
    Lh = 0.575
    Lc = 0.2
    s0_heater = W + 0.05
    s0_cooler = 2W + H + 0.05

# Grid
    N  = 20
    ds = L/N
    s  = 0:ds:(L-ds)

# Constants
    Qh = 3e3
    Qc = Qh
    A  = 0.25*π*D^2
    g = zeros(N)
    g[W .< s .< W + H] .= -9.8
    g[2W + H .< s .< 2W + 2H] .= 9.8
    S = zeros(N)    
    S[s0_heater .< s .< s0_heater + Lh] .=  Qh/(A*Lh)
    S[s0_cooler .< s .< s0_cooler + Lc] .= -Qc/(A*Lc)
    f = 0.02*ones(N)

# Periodic central difference matrix
    # D = spdiagm(
    #     -1 => -ones(N-1),
    #      1 =>  ones(N-1)
    # )
    # D[1,N] = -1
    # D[N,1] =  1
    # D = 0.5/ds*D

# Periodic upwind difference matrix
    D = spdiagm(
        -1 => -ones(N-1),
         0 =>  ones(N)
    )
    # D[1,N] = -1
    D = D/ds

# Equation of state
    f_hat(h,p) = PropsSI("D", "H", h, "P", p, "Water")

# Root function
    function F(Q)
        matrix = reshape(Q,N,4)
        ρ,v,h,p = eachcol(matrix)

        F_mass = D*(ρ.*v);
        F_momentum = v.*D*v + 1 ./ ρ .* D*p + 0.5*f.*v.*abs.(v) - g;
        F_energy = D*(ρ.*h.*v) - v.*D*p - S;
        F_eos = ρ - f_hat.(h,p)

        F_mass[1]     = ρ[1] - ρ0
        F_momentum[1] = v[1] - v0
        F_energy[1]   = h[1] - h0
        F_eos[1]      = p[1] - p0

        return [F_mass; F_momentum; F_energy; F_eos]
    end

# Jacobian
    e(j) = I(4N)[:,j]
    δ = 1e-6
    J(Q) = hcat(((F(Q + δ*e(j)) - F(Q))/δ for j in 1:4N)...)


# Initialize
    v0 = 0.5        
    p0 = 57.98e5    

    T_sat = PropsSI("T", "P", p0, "Q", 0, "Water")
    T_subcooling = 12
    T0 = T_sat - T_subcooling
    h0 = PropsSI("H", "P", p0, "T", T0, "Water")
    ρ0 = PropsSI("D", "P", p0, "T", T0, "Water")

    ρ_init = ρ0*ones(N)
    v_init = v0*ones(N)
    h_init = h0*ones(N)
    p_init = p0*ones(N)

    Qk = [ρ_init; v_init; h_init; p_init]

# Main loop
    tol = 1e-2
    res = 1.0
    k = 0
    while res > tol
        global Qk, res

        k = k + 1 
        print(k)

        # matrix = reshape(Qk,N,4)
        # ρ,v,h,p = eachcol(matrix)
        # display(plot_this([ρ,v,h,p], ["ρ","v","h","p"]))

        x = -J(Qk)\F(Qk)
        Qk = Qk + x
        res = norm(x)
    end

println()

matrix = reshape(Qk,N,4)
ρ,v,h,p = eachcol(matrix)
p = p/1e5
h = h/1e3
display(plot_this([ρ,v,h,p], ["ρ","v","h","p"]))

end