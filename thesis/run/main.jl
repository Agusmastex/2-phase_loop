using DelimitedFiles
using Plots
using LaTeXStrings

models = ["HEM3", "HFM4SG", "HFM4SZ"]
databases = ["relap-2016"]

for model_name in models
    include("models\\" * model_name * ".jl")
end

root = "C:\\Users\\mateo\\Files\\Investigación\\2-phase_loop\\thesis\\run"

function read_conditions(conditions_filepath)
    input_file = open(conditions_filepath, "r")
    lines = readlines(input_file)
    conditions = [parse(Float64, split(line, "\t")[1]) for line in lines]
    conditions[1] = conditions[1]*1e3   # kPa to Pa
    conditions[4] = conditions[4]*1e3   # kW/m2 to W/m2
    return conditions
end

function run_model(model_name, conditions, n_nodes)
    model = getfield(Main, Symbol(model_name))
    fields = model.run(conditions..., n_nodes=n_nodes)

    header = hcat(keys(fields)...)
    matrix = hcat(values(fields)...)

    [header; matrix]
end

function produce_data(models, databases, n_nodes)
    for database in databases
        database_folder = root * "\\data\\" * database
        cd(database_folder)
        for scenario in readdir()
            cd(scenario)

            conditions = read_conditions("conditions.txt")

            for model_name in models
                results_path = root * "\\results\\" * "\\" * database * "\\" * scenario
                try 
                    cd(results_path)
                catch
                    mkpath(results_path)
                    cd(results_path)
                end

                results_name = database * "_" * model_name * "_" * scenario
                println(results_name)
                csv = run_model(model_name, conditions, n_nodes)
                writedlm(results_name * ".csv", csv, ',')
                println()
            end
            cd(database_folder)
            println()
        end
    end
end

function plot_data3(database)
    database_folder = root * "\\data\\" * database
    database_results_folder = root * "\\results\\" * database
    cd(database_folder)
    for scenario in readdir()
        α_plot = plot(title = L"\alpha", ylims=(-0.03,0.8))
        T_plot = plot(title = L"T")
        p_plot = plot(title = L"P")

        cd(database_folder * "\\" *scenario)

        input_file = open("conditions.txt", "r")
        lines = readlines(input_file)
        conditions = [split(line, "\t")[1] for line in lines]

        plot_title = "P_{\\mathrm{in}} = $(conditions[1]) \\; \\mathrm{kPa} \\quad v_{\\mathrm{in}} = $(conditions[2]) \\; \\mathrm{m/s} \\quad ΔT_{\\mathrm{sub}} = $(conditions[3]) \\mathrm{ºC} \\quad q^{''} = $(conditions[4]) \\; \\mathrm{kW/m}^2"
        plot_title = latexstring(plot_title)

        (z_α, α) = eachcol(readdlm("alpha.csv", ','))
        (z_T, T) = eachcol(readdlm("T.csv", ','))
        (z_p, p) = eachcol(readdlm("P.csv", ','))

        d_inner = 0.0191
        d_outer = 0.0381
        Dh = d_outer - d_inner

        z_α = z_α*Dh
        z_T = z_T*Dh
        z_p = z_p*Dh

        plot!(α_plot, z_α, α, label = "ref")
        plot!(T_plot, z_T, T, label = "ref")
        plot!(p_plot, z_p, p, label = "ref")

        Lh = 2.8
        vline!(α_plot, [Lh], color=:black, label=nothing)
        vline!(T_plot, [Lh], color=:black, label=nothing)
        vline!(p_plot, [Lh], color=:black, label=nothing)

        cd(database_results_folder * "\\" * scenario)

        for result_name in readdir()
            matrix = readdlm(result_name, ',')
            header = matrix[1,:]
            fields = eachcol(matrix[2:end,:])    
            field_dict = Dict(zip(header, fields))

            model = split(result_name, "_")[2]

            plot!(α_plot, field_dict["z"], field_dict["alpha"], label = model)
            plot!(T_plot, field_dict["z"], field_dict["T"], label = model)
            plot!(p_plot, field_dict["z"], field_dict["p"], label = model)
            
        end
        
        title = plot(title = plot_title, grid = false, showaxis = false, bottom_margin = -50Plots.px, titlefontsize=11)
        l = @layout [
         grid(1,3)
         a{0.01h}
        ]   
        this_plot = plot(α_plot, T_plot, p_plot, title, layout=l)
        plot_name = root * "\\plots\\" * "\\" * database * "\\" * scenario
        try
            savefig(this_plot, plot_name)
        catch
            mkpath(root * "\\plots\\" * "\\" * database)
            savefig(this_plot, plot_name)
        end
    end

end

function test_mesh_independence(model_name, conditions, nodes_array)
        println(model_name)
    for n_nodes in nodes_array
        println("n = $n_nodes")
        @time csv = run_model(model_name, conditions, n_nodes)
        println()
        results_filepath = root * "\\results\\mesh_independence\\" * "\\" * model_name * "\\"  * model_name * "_n=$n_nodes" * ".csv"
        try
            writedlm(results_filepath, csv, ',')
        catch
            mkpath(root * "\\results\\mesh_independence\\" * model_name)
            writedlm(results_filepath, csv, ',')
        end
    end
end

function plot_single(plot, csv_path, plotted_field; label = nothing)
    matrix = readdlm(csv_path, ',')
    header = matrix[1,:]
    fields = eachcol(matrix[2:end,:])    
    field_dict = Dict(zip(header, fields))
    plot!(
        plot, 
        field_dict["z"],
        field_dict[plotted_field],
        label = label,
    )
end

# test_mesh_independence("HFM4_SG", [750e3, 1.0, 11, 241e3], [5,10,20,40,80])

# cd(root * "\\results\\mesh_independence\\HFM4_SG")
# p = plot()
# for item in readdir()
#     plot_single(p, item, "P", label=item)
# end
# display(p)

# produce_data(models, databases, 50)

# plot_data3("relap-2016")

function four_plots(database, models)
    database_folder = root * "\\data\\" * database
    database_results_folder = root * "\\results\\" * database
    cd(database_folder)
    for scenario in readdir()
        α_plot = plot(title = L"\alpha", ylims=(-0.03,0.8))
        T_plot = plot(title = L"T")
        p_plot = plot(title = L"P")
        h_plot = plot(title = L"H")

        cd(database_folder * "\\" *scenario)

        input_file = open("conditions.txt", "r")
        lines = readlines(input_file)
        conditions = [split(line, "\t")[1] for line in lines]

        plot_title = "P_{\\mathrm{in}} = $(conditions[1]) \\; \\mathrm{kPa} \\quad v_{\\mathrm{in}} = $(conditions[2]) \\; \\mathrm{m/s} \\quad ΔT_{\\mathrm{sub}} = $(conditions[3]) \\mathrm{ºC} \\quad q^{''} = $(conditions[4]) \\; \\mathrm{kW/m}^2"
        plot_title = latexstring(plot_title)

        (z_α, α_ref) = eachcol(readdlm("alpha.csv", ','))
        (z_T, T_ref) = eachcol(readdlm("T.csv", ','))
        (z_p, p_ref) = eachcol(readdlm("P.csv", ','))

        d_inner = 0.0191
        d_outer = 0.0381
        Dh = d_outer - d_inner

        z_α = z_α*Dh
        z_T = z_T*Dh
        z_p = z_p*Dh

        plot!(α_plot, z_α, α_ref, label = "ref")
        plot!(T_plot, z_T, T_ref, label = "ref")
        plot!(p_plot, z_p, p_ref, label = "ref")

        Lh = 2.8
        vline!(α_plot, [Lh], color=:black, label=nothing)
        vline!(T_plot, [Lh], color=:black, label=nothing)
        vline!(p_plot, [Lh], color=:black, label=nothing)

        cd(database_results_folder * "\\" * scenario)

        for result_name in readdir()
            model = split(result_name, "_")[2]
            if model in models
            matrix = readdlm(result_name, ',')
            header = matrix[1,:]
            fields = eachcol(matrix[2:end,:])    
            field_dict = Dict(zip(header, fields))


            enthalpies = []
            temperatures = []
            enthalpy_labels = []
            temperature_labels = []

            try
                enthalpies = [field_dict["h"], field_dict["hl_sat"], field_dict["h_cr"]]
                temperatures = [field_dict["T"]]#, field_dict["T_sat"], field_dict["T_cr"]]
                enthalpy_labels = [L"H" L"H^{sat}_l" L"H_{cr}"]
                temperature_labels = [L"T" L"T_{sat}" L"T_{cr}"]
            catch
                enthalpies = [field_dict["h"], field_dict["hl_sat"]]
                temperatures = [field_dict["T"]]#, field_dict["T_sat"]]
                enthalpy_labels = [L"H" L"H^{sat}_l"]
                temperature_labels = [L"T" L"T_{sat}"]
            end

            enthalpy_labels = map(x -> x * " $model", enthalpy_labels)
            temperature_labels = map(x -> x * " $model", temperature_labels)

            z_scalar = field_dict["z"]
            α = field_dict["alpha"]
            p = field_dict["p"]

            plot!(
                h_plot,
                z_scalar, enthalpies,
                label = enthalpy_labels,
                title = L"H", 
                marker = 0,
                legend = :bottomright,
            )

            plot!(
                T_plot,
                z_scalar, temperatures,
                label = temperature_labels,
                title = L"T", 
                marker = 0,
                legend = :bottomright,
            )

            plot!(
                α_plot,
                z_scalar, α,
                title = L"\alpha", 
                label = model,
                marker = 0,
                ylims = (0, 0.8),
            )

            plot!(
                p_plot,
                z_scalar, p,
                title = L"P", 
                label = model,
                marker = 0,
            )
        end
            
        end
        
        title = plot(title = plot_title, grid = false, showaxis = false, bottom_margin = -50Plots.px, titlefontsize=11)
        l = @layout [
         grid(2,2)
         a{0.01h}
        ]   
        this_plot = plot(h_plot, α_plot, T_plot, p_plot, title, layout=l ,size = (1000, 1000))
        plot_name = root * "\\plots\\" * "\\" * database * "\\" * scenario
        try
            savefig(this_plot, plot_name)
        catch
            mkpath(root * "\\plots\\" * "\\" * database)
            savefig(this_plot, plot_name)
        end
    end

end

# produce_data(["HFM4SZ"], ["relap-2016"], 20)
# four_plots("relap-2016", ["HFM4SZ"])