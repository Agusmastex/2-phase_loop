using Plots
using LaTeXStrings
using TOML
using DelimitedFiles

# model = "HFM4SZ"
model = "HEM3"
n_nodes = 15

root = "C:\\Users\\mateo\\Files\\Investigación\\2-phase_loop\\thesis\\run"
include("C:\\Users\\mateo\\Files\\Investigación\\2-phase_loop\\thesis\\run\\models\\" * model * ".jl")

for scenario in readdir(root * "\\data")

    println(scenario)

# Conditions and geometry input
    cd(root * "\\data\\" * scenario)

    conditions = TOML.parsefile("conditions.toml")
    P_in = conditions["P_in"]
    # v_in = conditions["v_in"]
    ΔT_sub = conditions["T_sub"]
    q_flux = conditions["q_flux"]

    geom = TOML.parsefile("geometry.toml")
    L  = geom["L"]
    Lh = geom["Lh"]
    d_inner = geom["d_inner"]
    d_outer = geom["d_outer"]
    Dh = d_outer - d_inner

# Run and unpack

   field_dict =  HEM3.run(conditions, geom, n_nodes=n_nodes)
   z      = field_dict["z"]
   α      = field_dict["alpha"]
   T      = field_dict["T"]
   T_sat  = field_dict["T_sat"]
#    T_cr   = field_dict["T_cr"]
   h      = field_dict["h"]
   hl_sat = field_dict["hl_sat"]
#    h_cr   = field_dict["h_cr"]
   p      = field_dict["p"]
#    Pe     = field_dict["Pe"]
#    Γ      = field_dict["Gamma"]
   v_in   = round(field_dict["v"][1], digits=2)

# Plotting 

    enthalpies = [h, hl_sat]#, h_cr]
    enthalpies = map(x -> x/1e3, enthalpies)
    temperatures = [T, T_sat]#, T_cr]
    temperatures = map(x -> x .- 273.15, temperatures)
    p = p/1e3

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
        z, p,
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

    cd(root * "\\data\\" * scenario)
    (z_α, α_ref) = eachcol(readdlm("alpha.csv", ','))
    (z_T, T_ref) = eachcol(readdlm("T.csv", ','))
    (z_p, p_ref) = eachcol(readdlm("P.csv", ','))

    z_α = z_α*Dh
    z_T = z_T*Dh
    z_p = z_p*Dh

    scatter!(void_plot,        z_α, α_ref, label = "ref", marker=:star5)#, markersize=2)
    scatter!(temperature_plot, z_T, T_ref, label = "ref", marker=:star5)#, markersize=2)
    scatter!(pressure_plot,    z_p, p_ref, label = "ref", marker=:star5)#, markersize=2)

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