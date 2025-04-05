using LinearAlgebra
using SparseArrays
using CoolProp
using Plots

""" 
Solves the 2-phase steady state vertical pipe flow with constant upstream boundary conditions
using the 3-equation Homogeneous Equilibrium Model (HEM)
upwind finite difference
and manual implementation of Newton-Raphson employing numerical jacobian
"""

# Auxiliary functions 

  function plot_this(list, names)
    plots = [plot(z, item, 
                title=name, 
                marker=:circle,
                markersize=1.5,
                linewidth=0) 
            for (item, name) in zip(list, names)]
    return plots
  end

# Geometry
    D = 9.1e-3
    L = 0.5
    Lh = 0.2*L
    z0_heater = (L - Lh)/2

# Grid
    N  = 20
    dz = L/N
    z  = 0:dz:(L-dz)

# Constants
    Qh = 1225
    A  = 0.25*π*D^2
    g  = 9.8
    S  = zeros(N)    
    S[z0_heater .< z .< z0_heater + Lh] .=  Qh/(A*Lh)
    f  = 0.02*ones(N)
    ΔT_subcooling = 10

# Upwind difference matrix
    D = spdiagm(
        -1 => -ones(N-1),
         0 =>  ones(N)
    )
    D = D/dz

# Equation of state
    f_hat(h,p) = PropsSI("D", "H", h, "P", p, "Water")

# Root function
    function F(Q)
        matrix = reshape(Q,N,4)
        ρ,v,h,p = eachcol(matrix)

        F_mass = D*(ρ.*v);
        F_momentum = v.*D*v + 1 ./ ρ .* D*p + 0.5*f.*v.*abs.(v) .+ g;
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
    p0 = 1e5

    T_sat = PropsSI("T", "P", p0, "Q", 0, "Water")
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
        global Qk, res, k

        k = k + 1 
        print(k)

        x = -J(Qk)\F(Qk)
        Qk = Qk + x
        res = norm(x)
    end

# Derived fields

    matrix = reshape(Qk,N,4)
    ρ,v,h,p = eachcol(matrix)

    T = PropsSI.("T", "P", p, "H", h, "Water") .- 273.15
    Q = PropsSI.("Q", "P", p, "H", h, "Water")
    Q[Q .== -1] .= 0
    ρg = PropsSI.("D", "P", p, "Q", 1, "Water")
    ρl = PropsSI.("D", "P", p, "Q", 0, "Water")
    α = (Q./ρg)./(Q./ρg + (1 .- Q)./ρl)
    p = p .- p[end]
    h = h/1e3

# Plotting 

    field_dict = Dict(
    "ρ" => ρ,
    "v" => v,
    "h" => h,
    "p" => p,
    "T" => T,
    "Q" => Q,
    "α" => α,
    "ρg" => ρg,
    "ρl" => ρl
    )

    select = ["ρ","α","T","p"]

    fields = [field_dict[key] for key in select]
    plots = plot_this(fields, select)
    for p in plots
        vline!(p, [z0_heater,z0_heater+Lh])
    end
    plot(plots...)