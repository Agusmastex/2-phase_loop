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
    Lh = L
    z0_heater = (L - Lh)/2

# Grid
    N  = 20
    dz = L/N
    z  = 0:dz:(L-dz)

# Constants
    # Upstream conditions
    v0 = 0.5        
    p0 = 1e5
    ΔT_subcooling = 0.1
    T0 = T_sat - ΔT_subcooling
    h0 = PropsSI("H", "P", p0, "T", T0, "Water")
    ρ0 = PropsSI("D", "P", p0, "T", T0, "Water")
    # Thermodynamic constants
    T_sat = PropsSI("T", "P", p0, "Q", 0, "Water")
    hfg = 
        PropsSI("H", "P", p0, "Q", 1, "Water") - 
        PropsSI("H", "P", p0, "Q", 0, "Water")
    Cp = PropsSI("CPMASS", "P", p0, "Q", 1, "Water")
    # True constants
    A  = 0.25*π*D^2
    g  = 0
    # Heat input
    x_target = 0.01
    m_dot = ρ0*v0*A
    Qh = m_dot*(x_target*hfg + Cp*ΔT_subcooling)
    S  = zeros(N)    
    S[z0_heater .< z .< z0_heater + Lh] .=  Qh/(A*Lh)
    # Friction factor
    f  = 0.0*ones(N)

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
        F_energy = D*(ρ.*h.*v) - 0*v.*D*p - S;
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
    # p = p .- p[end]
    # h = h/1e3

# Analytical solution

H = h0 .+ Qh/m_dot*z/Lh
rhog = PropsSI("D", "P", p0, "Q", 1, "Water")
rhol = PropsSI("D", "P", p0, "Q", 0, "Water")
slope = Qh/m_dot/hfg*(rhol/rhog - 1)
V = v0 .+ v0*slope*z/Lh
P = p0 .- ρ0*v0^2*slope*z/Lh
RHO = m_dot./(A*V)

function myplot(field, name)
    return plot(z, field, title=name, marker=:circle, markersize=1.5)
end

p1 = myplot(ρ,"ρ")
plot!(z, RHO)
p2 = myplot(v,"v")
plot!(z, V)
p3 = myplot(h,"h")
plot!(z, H)
p4 = myplot(p,"p")
plot!(z, P)

plots = [p1, p2, p3, p4]
plot(plots...)