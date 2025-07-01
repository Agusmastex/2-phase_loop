module HEM3

using LinearAlgebra
using SparseArrays
using CoolProp
using Plots

""" 
Solves the 2-phase steady state vertical pipe flow with constant upstream boundary conditions
using the 3-equation Homogeneous Equilibrium Model (HEM)
collocated upwind finite difference
nondimensional version of the equations
and manual implementation of Newton-Raphson employing numerical jacobian
"""

# function run(P_in, G, ΔT_sub, q_flux) # mass flux given
function run(P_in, j_in, ΔT_sub, q_flux; n_nodes=20) # velocity given

# Geometry
    L = 5.03
    d_inner = 19.05e-3
    d_outer = 38.10e-3
    Dh = d_outer - d_inner
    Lh = 3.0/L
    z0_heater = 0.0/L

# Grid
    N  = n_nodes
    dz = 1/N
    z  = 0:dz:(1-dz)

# Constants
    # Upstream conditions
    p0 = P_in
    T_sat = PropsSI("T", "P", p0, "Q", 0, "Water")
    T0 = T_sat - ΔT_sub
    h0 = PropsSI("H", "P", p0, "T", T0, "Water")
    ρ0 = PropsSI("D", "P", p0, "T", T0, "Water")
    # v0 = G/ρ0 # mass flux given
    v0 = j_in # velocity given
    # True constants
    A = 0.25*π*(d_outer^2 - d_inner^2)
    g = 9.8
    rugosity = 0
    # Heat input
    Qh = q_flux*Lh*L*π*d_inner
    S  = zeros(N)    
    S[z0_heater .< z .< z0_heater + Lh] .=  1
    # Derived quantities
    m_dot = ρ0*v0*A
    h_Lh = h0 + Qh/m_dot
    # Dimensionless groups
    Fr = v0^2/(g*L)
    Ec = v0^2/(h_Lh - h0)

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
    f1 = (-2.457*log((7 / Re)^0.9 + 0.27*rugosity/Dh))^16
    f2 = (37530 / Re)^12
    return 8*((8 / Re)^12 + (f1 + f2)^(-1.5))^(1/12)
  end


# Root function
    function F(Q)
        matrix = reshape(Q,N,4)
        ρ,v,h,p = eachcol(matrix)

        μ  = PropsSI.("V", "H", h0 .+ (h_Lh - h0)*h,"P", p0 .+ ρ0*g*L*p, "Water")
        Re = ρ0*v0 * ρ.*v*Dh./μ

        F_mass = D*(ρ.*v)
        F_momentum = v.*D*v + 1/Fr * 1 ./ ρ .* D*p + 0.5*f.(Re)*L/Dh.*v.*abs.(v) .+ 1/Fr
        F_energy = D*(ρ.*h.*v) - Ec/Fr * v.*D*p - S/Lh
        F_eos = ρ0*ρ - f_hat.(h0 .+ (h_Lh - h0)*h, p0 .+ ρ0*g*L*p)

        F_mass[1]     = ρ[1] - 1
        F_momentum[1] = v[1] - 1
        F_energy[1]   = h[1] - 0
        F_eos[1]      = p[1] - 0

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
    ρ_init = ones(N)
    v_init = ones(N)
    h_init = zeros(N)
    p_init = zeros(N)

    Qk = [ρ_init; v_init; h_init; p_init]

# Main loop
    tol = 1e-2
    res = 1.0
    k = 0
    while res > tol
        # global Qk, res, k

        k = k + 1 
        print(k)

        x = -J(Qk)\F(Qk)
        Qk = Qk + x
        res = norm(x)
    end

# Derived fields

    matrix = reshape(Qk,N,4)
    ρ,v,h,p = eachcol(matrix)

    ρ = ρ0*ρ
    v = v0*v
    h = h0 .+ (h_Lh - h0)*h
    p = p0 .+ ρ0*g*L*p

    T = PropsSI.("T", "P", p, "H", h, "Water") .- 273.15
    Q = PropsSI.("Q", "P", p, "H", h, "Water")
    Q[Q .== -1] .= 0
    ρg = PropsSI.("D", "P", p, "Q", 1, "Water")
    ρl = PropsSI.("D", "P", p, "Q", 0, "Water")
    α = (Q./ρg)./(Q./ρg + (1 .- Q)./ρl)
    # p = p .- p[end]
    p = p/1e3
    h = h/1e3

    dz = L/N
    z  = 0:dz:(L-dz)
    Lh = Lh*L
    z0_heater = z0_heater*L

    fields = Dict(
       "z" => z,
       "alpha" => α,
       "T" => T,
       "P" => p,
    )

return fields
end

end