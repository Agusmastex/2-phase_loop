using LinearAlgebra
using SparseArrays
using CoolProp
using Plots
using NonlinearSolve

""" 
Solves the 2-phase steady state vertical pipe flow with constant upstream boundary conditions
using the 3-equation Homogeneous Equilibrium Model (HEM)
upwind finite difference
and package NonlinearSolve.jl to solve the resulting system of 4N equations
(currently failing)
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
    d = 9.1e-3
    L = 0.5
    Lh = 0.2*L
    z0_heater = (L - Lh)/2

# Grid
    N  = 20
    dz = L/N
    z  = 0:dz:(L-dz)

# Constants
    # Upstream conditions
    v0 = 0.5        
    p0 = 1e5
    T_sat = PropsSI("T", "P", p0, "Q", 0, "Water")
    ΔT_subcooling = 10
    T0 = T_sat - ΔT_subcooling
    h0 = PropsSI("H", "P", p0, "T", T0, "Water")
    ρ0 = PropsSI("D", "P", p0, "T", T0, "Water")
    # True constants
    A = 0.25*π*d^2
    g = 9.8
    rugosity = 0
    # Heat input
    Qh = 1225
    S  = zeros(N)    
    S[z0_heater .< z .< z0_heater + Lh] .=  Qh/(A*Lh)

# Upwind difference matrix
    D = spdiagm(
        -1 => -ones(N-1),
         0 =>  ones(N)
    )
    D = D/dz

# Equation of state
    f_hat(h,p) = PropsSI("D", "H", h, "P", p, "Water")

# Friction factor

  function f(Re)
    f1 = (-2.457*log((7 / Re)^0.9 + 0.27*rugosity/d))^16
    f2 = (37530 / Re)^12
    return 8*((8 / Re)^12 + (f1 + f2)^(-1.5))^(1/12)
  end


# Root function
    function F(Q,p,t)
        matrix = reshape(Q,N,4)
        ρ,v,h,p = eachcol(matrix)

        μ  = PropsSI.("V","H",h,"P",p,"Water")
        Re = ρ.*v*d./μ

        F_mass = D*(ρ.*v);
        F_momentum = v.*D*v + 1 ./ ρ .* D*p + 0.5*f.(Re)/d.*v.*abs.(v) .+ g;
        F_energy = D*(ρ.*h.*v) - v.*D*p - S;
        F_eos = ρ - f_hat.(h,p)

        F_mass[1]     = ρ[1] - ρ0
        F_momentum[1] = v[1] - v0
        F_energy[1]   = h[1] - h0
        F_eos[1]      = p[1] - p0

        return [F_mass; F_momentum; F_energy; F_eos]
    end

# Initialize
    ρ_init = ρ0*ones(N)
    v_init = v0*ones(N)
    h_init = h0*ones(N)
    p_init = p0*ones(N)

    Q0 = [ρ_init; v_init; h_init; p_init]

# Solve using NonlinearSolve.jl

problem = NonlinearProblem(F, Q0, nothing)
#problem = IntervalNonlinearProblem(F,Q0-)
solution = solve(problem, FastShortcutNonlinearPolyalg(), show_trace = Val(true))
Q = solution.u

# Derived fields

    matrix = reshape(Q,N,4)
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