include("C:\\Users\\mateo\\Files\\Investigación\\2-phase_loop\\thesis\\run\\models\\HFM4SZ.jl")


""" 
Solves the 2-phase steady state 1D vertical flow with constant upstream boundary conditions
using the 4-equation Homogeneous Flow Model (HFM)
the Saha-Zuber correlation for critical enthalpy of the point of Net Vapor Generation
the Lahey method for subcooled boiling source term in the vapor mass equation
upwind staggered finite volume
nondimensional version of the equations
and manual implementation of Newton-Raphson employing numerical jacobian
"""

# Conditions input
    P_in = 750e3
    v_in = 1.0
    ΔT_sub = 11
    q_flux = 241e3
    n_nodes = 20

# Geometry
    struct Geom
        L
        Dh
        Lh
        z0_heater
        A_flow
        A_wall
    end

    L = 2.8 + 1.7
    d_inner = 0.0191
    d_outer = 0.0381
    Dh = d_outer - d_inner
    Lh = 2.8
    z0_heater = 0.0
    A_flow = 0.25*π*(d_outer^2 - d_inner^2)
    A_wall = π*d_inner*Lh
    geom = Geom(L, Dh, Lh, z0_heater, A_flow, A_wall)


# Run and unpack

   field_dict =  HFM4SZ.run(P_in, v_in, ΔT_sub, q_flux, geom, n_nodes=n_nodes)
   z_scalar = field_dict["z"]
   α = field_dict["alpha"]
   T = field_dict["T"]
   T_sat = field_dict["T_sat"]
   T_cr = field_dict["T_cr"]
   h = field_dict["h"]
   hl_sat = field_dict["hl_sat"]
   h_cr = field_dict["h_cr"]
   p = field_dict["p"]
   Pe = field_dict["Pe"]

# Report

    println()
    println("h_cr/hl_sat = $(maximum(h_cr./hl_sat))")
    println("Pe = $(Int(round(sum(Pe)/(n_nodes+1))))")

# Plotting 

    using Plots
    using LaTeXStrings

    enthalpies = [h, hl_sat, h_cr]

    temperatures = [T, T_sat, T_cr]


    enthalpy_plot = plot(
        z_scalar, enthalpies,
        label = [L"H" L"H^{sat}_l" L"H_{cr}"],
        title = L"H", 
        marker = [1 0 0],
        legend = :bottomright,
    )

    temperature_plot = plot(
        z_scalar, temperatures,
        label = [L"T" L"T_{sat}" L"T_{cr}"],
        title = L"T", 
        marker = [1 0 0],
        legend = :bottomright,
    )

    void_plot = plot(
        z_scalar, α,
        title = L"\alpha", 
        label = nothing,
        marker = 1,
        ylims = (0, 0.8),
    )

    pressure_plot = plot(
        z_scalar, p/1e3,
        title = L"P", 
        label = nothing,
        marker = 1,
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

    using DelimitedFiles

    database_path = "C:\\Users\\mateo\\Files\\Investigación\\2-phase_loop\\thesis\\run\\data\\relap-2016\\"
    cd(database_path * "fig-2")
    (z_α, α_ref) = eachcol(readdlm("alpha.csv", ','))
    (z_T, T_ref) = eachcol(readdlm("T.csv", ','))
    (z_p, p_ref) = eachcol(readdlm("alpha.csv", ','))

    z_α = z_α*Dh
    z_T = z_T*Dh
    z_p = z_p*Dh

    plot!(void_plot,        z_α, α_ref, label="ref")
    plot!(temperature_plot, z_T, T_ref, label="ref")
    plot!(pressure_plot,    z_p, p_ref, label="ref")

    plot(plots...)