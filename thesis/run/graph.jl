using Plots
using LaTeXStrings
using TOML
using DelimitedFiles
using CSV

root = "C:\\Users\\mateo\\Files\\Investigación\\2-phase_loop\\thesis\\run"

function prepare_plot(scenario)
    enthalpy_plot = plot(
        title = L"H", 
        legend = :bottomright,
    )

    temperature_plot = plot(
        title = L"T", 
        legend = :bottomright,
    )

    void_plot = plot(
        title = L"\alpha", 
        ylims = (0, 0.8),
    )

    pressure_plot = plot(
        title = L"P", 
    )

    cd(root * "\\data\\" * scenario)

    conditions = TOML.parsefile("conditions.toml")
    P_in = conditions["P_in"]
    v_in = round(conditions["G"]/1000, digits=2)
    ΔT_sub = conditions["T_sub"]
    q_flux = conditions["q_flux"]

    geom = TOML.parsefile("geometry.toml")
    Dh = geom["d_outer"] - geom["d_inner"]

    (z_α, α_ref) = eachcol(readdlm("alpha.csv", ','))
    (z_T, T_ref) = eachcol(readdlm("T.csv", ','))
    (z_p, p_ref) = eachcol(readdlm("P.csv", ','))

    z_α = z_α*Dh
    z_T = z_T*Dh
    z_p = z_p*Dh

    scatter!(void_plot,        z_α, α_ref, label = "ref", marker=:star5)#, markersize=1.5)
    scatter!(temperature_plot, z_T, T_ref, label = "ref", marker=:star5)#, markersize=1.5)
    scatter!(pressure_plot,    z_p, p_ref, label = "ref", marker=:star5)#, markersize=1.5)

    plot_title = "P_{\\mathrm{in}} = $(P_in/1e3) \\; \\mathrm{kPa} \\quad v_{\\mathrm{in}} = $v_in \\; \\mathrm{m/s} \\quad ΔT_{\\mathrm{sub}} = $ΔT_sub \\mathrm{ºC} \\quad q^{''} = $(q_flux/1e3) \\; \\mathrm{kW/m}^2"
    plot_title = latexstring(plot_title)

    plots = Dict(
        "H" => enthalpy_plot,
        "T" => temperature_plot,
        "α" => void_plot,
        "P" => pressure_plot,
    )

    return plots, plot_title
end

function add_plot(plots, result, model_name)
    data = CSV.File(root * "\\results\\" * result)

    z = data.z

    p = data.p
    α = data.alpha

    h = data.h
    hl_sat = data.hl_sat
    h_cr = data.h_cr

    T = data.T
    T_sat = data.T_sat

    enthalpies = [h, hl_sat, h_cr]
    enthalpies = map(x -> x/1e3, enthalpies)
    temperatures = [T, T_sat]
    temperatures = map(x -> x .- 273.15, temperatures)
    p = p/1e3

    plot!(plots["H"], z, enthalpies  , marker=:circle, markersize=1.5, markerstrokewidth=0.1, label=nothing)#, label=[L"H" L"H_{sat}" L"H_{cr}"])
    plot!(plots["T"], z, temperatures, marker=:circle, markersize=1.5, markerstrokewidth=0.1, label=nothing)#, label=[L"T" L"T_{sat}"])
    plot!(plots["α"], z, α           , marker=:circle, markersize=1.5, label=model_name, markerstrokewidth=0.1)
    plot!(plots["P"], z, p           , marker=:circle, markersize=1.5, label=model_name, markerstrokewidth=0.1)

    return plots
end
function finalize_plot(plots, plot_title)

    title = plot(title = plot_title, grid = false, showaxis = false, bottom_margin = -50Plots.px, titlefontsize=11)
    l = @layout [
     grid(2,2)
     a{0.01h}
    ]   
    final_plot = plot(
        plots["H"], plots["α"], plots["T"], plots["P"], 
        title, layout=l, size = (600,500))
    display(final_plot)
end


for scenario in readdir(root * "\\data")
    try
    plots, plot_title = prepare_plot(scenario)
    plots = add_plot(plots, scenario * "_HFM4SZ.csv", "SZ")
    plots = add_plot(plots, scenario * "_HFM4SZ_C.csv", "C")
    finalize_plot(plots, plot_title)
    catch
        println(scenario)
    end
end