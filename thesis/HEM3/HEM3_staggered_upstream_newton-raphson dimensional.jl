using LinearAlgebra
using SparseArrays
using CoolProp
using Plots

""" 
Solves the 2-phase steady state vertical pipe flow with constant upstream boundary conditions
using the 3-equation Homogeneous Equilibrium Model (HEM)
upwind staggered finite volume
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
    z_scalar = -dz/2:dz:(L - dz/2)
    z_vect   = 0:dz:L

# Constants
    # Upstream conditions
    v_half = 0.2        
    p_half = 1e5
    T_sat = PropsSI("T", "P", p_half, "Q", 0, "Water")
    ΔT_subcooling = 10
    T_half = T_sat - ΔT_subcooling
    h_half = PropsSI("H", "P", p_half, "T", T_half, "Water")
    ρ_half = PropsSI("D", "P", p_half, "T", T_half, "Water")
    # True constants
    A = 0.25*π*(d_pipe^2 - d_rod^2)
    g = 9.8
    rugosity = 0
    # Heat input
    Qh = 1225
    S  = zeros(N+1)    
    S[z0_heater .< z_scalar .< z0_heater + Lh] .=  Qh/(A*Lh)

# Upwind difference matrix
    D_upw = spdiagm(
        -1 => -ones(N),
         0 =>  ones(N+1)
    )
    D_upw = D_upw/dz

# Centered difference matrix

    D_center = spdiagm(
        -1 => -ones(N),
         1 =>  ones(N)
    )
    D_center[N+1,:] = [zeros(N-1); -2; 2]
    D_center = D_center/(2*dz)

# Mean matrix

    M = spdiagm(
        0 => ones(N+1),
        1 => ones(N)
    )
    M[N+1,:] = [zeros(N-1); -1; 3]
    M = M/2


# Equation of state
    f_hat(h,p) = PropsSI("D", "H", h, "P", p, "Water")

# Friction factor

  function f(Re)
    f1 = (-2.457*log((7 / Re)^0.9 + 0.27*rugosity/Dh))^16
    f2 = (37530 / Re)^12
    # return 8*((8 / Re)^12 + (f1 + f2)^(-1.5))^(1/12)
    return 0.2
  end


# Root function
    function F(Q)
        matrix = reshape(Q,N+1,4)
        ρ,v,h,p = eachcol(matrix)

        μ  = PropsSI.("V","H",h,"P",p,"Water")
        Re = ρ.*v*Dh./μ

        F_mass = D_upw*(ρ.*v);
        F_momentum = v.*D_upw*v + (D_center*p)./(M*ρ) + 0.5*f.(Re)/Dh.*v.*abs.(v) .+ g;
        F_energy = D_upw*(ρ.*h.*v) - (M*v).*(D_center*p) - S;
        F_eos = ρ - f_hat.(h,p)

        F_mass[1]     = ρ[1] + ρ[2] - 2*ρ_half
        F_momentum[1] = v[1] - v_half 
        F_energy[1]   = h[1] + h[2] - 2*h_half
        F_eos[1]      = p[1] + p[2] - 2*p_half

        return [F_mass; F_momentum; F_energy; F_eos]
    end

# Jacobian
    function J(Q)
        e(j) = I(4(N+1))[:,j]
        δ = 1e-6
        FQk = F(Q)
        Jac = hcat(((F(Q + δ*e(j)) - FQk)/δ for j in 1:4(N+1))...)
        return Jac
    end


# Initialize
    ρ_init = ρ_half*ones(N+1)
    v_init = v_half*ones(N+1)
    h_init = h_half*ones(N+1)
    p_init = p_half*ones(N+1)

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

    matrix = reshape(Qk,N+1,4)
    ρ,v,h,p = eachcol(matrix)

    T = PropsSI.("T", "P", p, "H", h, "Water") .- 273.15
    Q = PropsSI.("Q", "P", p, "H", h, "Water")
    Q[Q .== -1] .= 0
    ρg = PropsSI.("D", "P", p, "Q", 1, "Water")
    ρl = PropsSI.("D", "P", p, "Q", 0, "Water")
    α = (Q./ρg)./(Q./ρg + (1 .- Q)./ρl)
    p = p .- p[end]
    p = p/1e3 #kPa
    h = h/1e3 #kJ


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


p1 = plot(z_scalar, ρ)
p2 = plot(z_vect,   v)
p3 = plot(z_scalar, h)
p4 = plot(z_scalar, p)

plot(p1,p2,p3,p4)

#     select = ["ρ","α","T","p"]

#     fields = [field_dict[key] for key in select]
#     plots = plot_this(fields, select)
#     for p in plots
#         vline!(p, [z0_heater,z0_heater+Lh])
#     end
#  p = plot(plots...)