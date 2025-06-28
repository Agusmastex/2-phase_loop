using LinearAlgebra
using SparseArrays
using CoolProp
using Plots

""" 
Solves the 2-phase steady state vertical pipe flow with constant upstream boundary conditions
using the 3-equation Homogeneous Equilibrium Model (HEM)
upwind finite difference
primitive dimensional variables
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
    d_pipe = 38.1e-3
    d_rod  = 19.05e-3
    Dh = d_pipe - d_rod
    L = 5.19
    Lh = 3
    z0_heater = 0.203

# Grid
    N  = 20
    dz = L/N
    z  = 0:dz:(L-dz)

# Constants
    # Upstream conditions
    v0 = 0.2        
    p0 = 1e5
    T_sat = PropsSI("T", "P", p0, "Q", 0, "Water")
    ΔT_subcooling = 10
    T0 = T_sat - ΔT_subcooling
    h0 = PropsSI("H", "P", p0, "T", T0, "Water")
    ρ0 = PropsSI("D", "P", p0, "T", T0, "Water")
    # True constants
    A = 0.25*π*(d_pipe^2 - d_rod^2)
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
    # D = collect(D)

# Equation of state
    f_hat(h,p) = PropsSI("D", "H", h, "P", p, "Water")

# Friction factor

  function f(Re)
    f1 = (-2.457*log((7 / Re)^0.9 + 0.27*rugosity/Dh))^16
    f2 = (37530 / Re)^12
    return 8*((8 / Re)^12 + (f1 + f2)^(-1.5))^(1/12)
  end


# Root function
    function F(Q)
        matrix = reshape(Q,N,4)
        ρ,v,h,p = eachcol(matrix)

        μ  = PropsSI.("V","H",h,"P",p,"Water")
        Re = ρ.*v*Dh./μ

        F_mass = D*(ρ.*v);
        F_momentum = v.*D*v + 1 ./ ρ .* D*p + 0.5*f.(Re)/Dh.*v.*abs.(v) .+ g;
        F_energy = D*(ρ.*h.*v) - v.*D*p - S;
        F_eos = ρ - f_hat.(h,p)

        F_mass[1]     = ρ[1] - ρ0
        F_momentum[1] = v[1] - v0
        F_energy[1]   = h[1] - h0
        F_eos[1]      = p[1] - p0

        return [F_mass; F_momentum; F_energy; F_eos]
    end

# Jacobian
    function J(Q)
        e(j) = I(4N)[:,j]
        δ = 1e-6
        FQk = F(Q)
        Jac = hcat(((F(Q + δ*e(j)) - FQk)/δ for j in 1:4N)...)
        return Jac
    end


# Initialize
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
    p = p/1e2 #kPa
    # h = h/1e3 #kJ

# Experimental data
    ports = [0, 79.3, 128.7, 168, 214.3, 247.3]
    z_exp = ports*Dh
    T_exp = [
        124.77099193628348,
        131.9847321842359 ,
        141.03053432113606,
        148.8167927627915 ,
        149.04580213690016,
        148.7595417296565 ,
    ]
    P_exp = [
        496.7938914759427 ,
        482.8244262338762 ,
        474.12213641094945,
        467.4809152302949 ,
        460.38167879580203,
        455.49618042527464,
    ]


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
 p = plot(plots...)
#  display(p)

# p1 = plot(z,T, title="T")
# plot!(z_exp, T_exp, marker=:circle, linewidth=0)
# p2 = plot(z,p, title="P")
# plot!(z_exp, P_exp, marker=:circle, linewidth=0)
# plot(p1, p2)
