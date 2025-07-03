using LinearAlgebra
using SparseArrays
using CoolProp
using Plots

""" 
Solves the 2-phase steady state vertical pipe flow with constant upstream boundary conditions
using the 3-equation Homogeneous Equilibrium Model (HEM)
upwind staggered finite volume
nondimensional version of the equations
and manual implementation of Newton-Raphson employing numerical jacobian
"""

# Conditions input
    P_in = 750e3
    v_in = 1.0
    ΔT_sub = 11
    q_flux = 241e3

# Geometry
    L = 2.8 + 1.7
    d_inner = 0.0191
    d_outer = 0.0381
    Dh = d_outer - d_inner
    Lh = 2.8/L
    z0_heater = 0.0/L

# Grid
    N  = 50
    dz = 1/N
    z_scalar = -dz/2:dz:(1 - dz/2)
    z_vect   = 0:dz:1

# Constants
    # Upstream conditions
    v_half = v_in
    p_half = P_in
    T_sat = PropsSI("T", "P", p_half, "Q", 0, "Water")
    T_half = T_sat - ΔT_sub
    h_half = PropsSI("H", "P", p_half, "T", T_half, "Water")
    ρ_half = PropsSI("D", "P", p_half, "T", T_half, "Water")
    # True constants
    A = 0.25*π*(d_outer^2 - d_inner^2)
    g = 9.8
    rugosity = 0
    # Heat input
    Qh = q_flux*Lh*L*π*d_inner
    S  = zeros(N+1)    
    S[z0_heater .< z_scalar .< z0_heater + Lh] .=  1
    # Derived quantities
    m_dot = ρ_half*v_half*A
    h_Lh = h_half + Qh/m_dot
    # Dimensionless groups
    Fr = v_half^2/(g*L)
    Ec = v_half^2/(h_Lh - h_half)

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

# Mean matrices

    M_ρ = spdiagm(
        0 => ones(N+1),
        1 => ones(N))
    M_ρ[N+1,:] = [zeros(N-1); -1; 3]
    M_ρ = M_ρ/2

    M_v = spdiagm(
        0 => ones(N+1),
        1 => ones(N))
    M_v = M_v/2
    


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
        matrix = reshape(Q,N+1,4)
        ρ,v,h,p = eachcol(matrix)

        h_dim = h_half .+ (h_Lh - h_half)*h
        p_dim = p_half .+ ρ_half*g*L*p

        μ  = PropsSI.("V","H",h_dim,"P",p_dim,"Water")
        Re = abs.(ρ_half*v_half*ρ.*v*Dh./μ)

        F_mass = D_upw*(ρ.*v)
        F_momentum = v.*D_upw*v + 1/Fr * (D_center*p)./(M_ρ*ρ) + 0.5*f.(Re)*L/Dh.*v.*abs.(v) .+ 1/Fr
        F_energy = D_upw*(ρ.*h.*v) - Ec/Fr * (M_v*v).*(D_center*p) - S/Lh
        F_eos = ρ_half*ρ - f_hat.(h_dim, p_dim)

        F_mass[1]     = ρ[1] + ρ[2] - 2
        F_momentum[1] = v[1] - 1
        F_energy[1]   = h[1] + h[2]
        F_eos[1]      = p[1] + p[2]

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
    ρ_init = ones(N+1)
    v_init = ones(N+1)
    h_init = zeros(N+1)
    p_init = zeros(N+1)

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

    ρ = ρ_half*ρ
    v = v_half*v
    h = h_half .+ (h_Lh - h_half)*h
    p = p_half .+ ρ_half*g*L*p

    T = PropsSI.("T", "P", p, "H", h, "Water")
    Q = PropsSI.("Q", "P", p, "H", h, "Water")
    Q[Q .== -1] .= 0
    ρg = PropsSI.("D", "P", p, "Q", 1, "Water")
    ρl = PropsSI.("D", "P", p, "Q", 0, "Water")
    α = (Q./ρg)./(Q./ρg + (1 .- Q)./ρl)

    hl_sat = PropsSI.("H", "P", p, "Q", 0, "Water")
    T_sat = PropsSI.("T", "P", p, "Q", 1, "Water")

# Plotting 

    using Plots
    using LaTeXStrings


    enthalpies = [h, hl_sat]
    enthalpies = map(x -> x/1e3, enthalpies)

    temperatures = [T, T_sat]
    temperatures = map(x -> x .- 273.15, temperatures)


    enthalpy_plot = plot(
        z_scalar, enthalpies,
        label = [L"H" L"H^{sat}_l"],
        title = L"H", 
        marker = [2 0],
        legend = :bottomright,
    )

    temperature_plot = plot(
        z_scalar, temperatures,
        label = [L"T" L"T_{sat}"],
        title = L"T", 
        marker = [2 0],
        legend = :bottomright,
    )

    void_plot = plot(
        z_scalar, α,
        title = L"\alpha", 
        label = nothing,
        marker = 2,
        ylims = (0, 0.8),
    )

    pressure_plot = plot(
        z_scalar, p,
        title = L"P", 
        label = nothing,
        marker = 0,
    )

plots = [
    enthalpy_plot,
    void_plot,
    temperature_plot,
    pressure_plot,
]

for plot in plots
    vline!(plot,  [z0_heater + Lh], color=:black, label=nothing)
end

plot(plots...)