include("C:\\Users\\mateo\\Files\\Investigación\\2-phase_loop\\thesis\\run\\models\\HFM4SZ.jl")
using Plots
using LaTeXStrings
using TOML

""" 
Solves the 2-phase steady state 1D vertical flow with constant upstream boundary conditions
using the 4-equation Homogeneous Flow Model (HFM)
the Saha-Zuber correlation for critical enthalpy of the point of Net Vapor Generation
the Lahey method for subcooled boiling source term in the vapor mass equation
upwind staggered finite volume
nondimensional version of the equations
and manual implementation of Newton-Raphson employing numerical jacobian
"""

for scenario in [
    # "fig-2",
    # "fig-3", 
    # "fig-4", 
    # "fig-7",
    # "fig-8",
    "fig-9"
    ]

n_nodes = 20

# Conditions and geometry input
    database_path = "C:\\Users\\mateo\\Files\\Investigación\\2-phase_loop\\thesis\\run\\data\\relap-2016\\"
    cd(database_path * scenario)

    conditions = TOML.parsefile("conditions.toml")
    P_in = conditions["P_in"]
    v_in = conditions["v_in"]
    ΔT_sub = conditions["T_sub"]
    q_flux = conditions["q_flux"]

    geom = TOML.parsefile("geometry.toml")
    L  = geom["L"]
    Lh = geom["Lh"]
    d_inner = geom["d_inner"]
    d_outer = geom["d_outer"]
    Dh = d_outer - d_inner


# Run and unpack

   field_dict =  HFM4SZ.run(conditions, geom, n_nodes=n_nodes)
   z      = field_dict["z"]
   α      = field_dict["alpha"]
   T      = field_dict["T"]
   T_sat  = field_dict["T_sat"]
   T_cr   = field_dict["T_cr"]
   h      = field_dict["h"]
   hl_sat = field_dict["hl_sat"]
   h_cr   = field_dict["h_cr"]
   p      = field_dict["p"]
   Pe     = field_dict["Pe"]
   Γ      = field_dict["Gamma"]

    println()

# Plotting 

    enthalpies = [h, hl_sat, h_cr]
    enthalpies = map(x -> x/1e3, enthalpies)
    temperatures = [T, T_sat, T_cr]
    temperatures = map(x -> x .- 273.15, temperatures)

    enthalpy_plot = plot(
        z, enthalpies,
        label = [L"H" L"H^{sat}_l" L"H_{cr}"],
        title = L"H", 
        marker = [1 0 0],
        legend = :bottomright,
    )

    temperature_plot = plot(
        z, temperatures,
        label = [L"T" L"T_{sat}" L"T_{cr}"],
        title = L"T", 
        marker = [1 0 0],
        legend = :bottomright,
    )

    void_plot = plot(
        z, α,
        title = L"\alpha", 
        label = nothing,
        marker = 1,
        ylims = (0, 0.8),
    )

    pressure_plot = plot(
        z, p/1e3,
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

    # Mark saturated boiling start
    using Roots
    using Dierckx
    root_f(zi) = Spline1D(z, hl_sat)(zi) - Spline1D(z, h)(zi)
    try
        z_satboil = find_zero(root_f, (0,L))
        for plot in plots
            vline!(plot,  [z_satboil], color=:black, label=nothing)
        end
    catch
        println("No saturated boiling")
    end

    # Mark heater end
    for plot in plots
        vline!(plot,  [Lh], color=:black, label=nothing)
    end

    using DelimitedFiles

    cd(database_path * scenario)
    (z_α, α_ref) = eachcol(readdlm("alpha.csv", ','))
    (z_T, T_ref) = eachcol(readdlm("T.csv", ','))
    (z_p, p_ref) = eachcol(readdlm("P.csv", ','))

    z_α = z_α*Dh
    z_T = z_T*Dh
    z_p = z_p*Dh

    plot!(void_plot,        z_α, α_ref, label = nothing)
    plot!(temperature_plot, z_T, T_ref, label = nothing)
    plot!(pressure_plot,    z_p, p_ref, label = nothing)

    plot_title = "P_{\\mathrm{in}} = $(P_in/1e3) \\; \\mathrm{kPa} \\quad v_{\\mathrm{in}} = $v_in \\; \\mathrm{m/s} \\quad ΔT_{\\mathrm{sub}} = $ΔT_sub \\mathrm{ºC} \\quad q^{''} = $(q_flux/1e3) \\; \\mathrm{kW/m}^2"
    plot_title = latexstring(plot_title)
    title = plot(title = plot_title, grid = false, showaxis = false, bottom_margin = -50Plots.px, titlefontsize=11)
    l = @layout [
     grid(2,2)
     a{0.01h}
    ]   
    final_plot = plot(plots..., title, layout=l, size = (600,500))
    display(final_plot)
end