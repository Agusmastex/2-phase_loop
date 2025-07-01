using DelimitedFiles
using Plots
using LaTeXStrings

models = ["HEM3", "HFM4_SG"]
databases = ["relap-2016"]
nodes = 100

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

function produce_data(models, databases)
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
                csv = run_model(model_name, conditions, nodes)
                writedlm(results_name * ".csv", csv, ',')
                println()
            end
            cd(database_folder)
            println()
        end
    end
end

function plot_data(database)
    database_folder = root * "\\data\\" * database
    database_results_folder = root * "\\results\\" * database
    cd(database_folder)
    for scenario in readdir()
        α_plot = plot(title = L"\alpha", ylims=(-0.03,0.8))
        T_plot = plot(title = L"T")
        p_plot = plot(title = L"P")

        cd(database_folder * "\\" *scenario)

        (z_α, α) = eachcol(readdlm("alpha.csv", ','))
        (z_T, T) = eachcol(readdlm("T.csv", ','))
        (z_p, p) = eachcol(readdlm("P.csv", ','))

        d_inner = 0.0191
        d_outer = 0.0381
        Dh = d_outer - d_inner

        z_α = z_α*Dh
        z_T = z_T*Dh
        z_p = z_p*Dh

        plot!(α_plot, z_α, α, label = "exp")
        plot!(T_plot, z_T, T, label = "exp")
        plot!(p_plot, z_p, p, label = "exp")

        cd(database_results_folder * "\\" * scenario)

        for result_name in readdir()
            matrix = readdlm(result_name, ',')
            header = matrix[1,:]
            fields = eachcol(matrix[2:end,:])    
            field_dict = Dict(zip(header, fields))

            model = split(result_name, "_")[2]

            plot!(α_plot, field_dict["z"], field_dict["alpha"], label = model)
            plot!(T_plot, field_dict["z"], field_dict["T"], label = model)
            plot!(p_plot, field_dict["z"], field_dict["P"], label = model)
            
        end
        this_plot = plot(α_plot, T_plot, p_plot, layout=(1,3))
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

test_mesh_independence("HFM4_SG", [750e3, 1.0, 11, 241e3], [5,10,20,40,80])

cd(root * "\\results\\mesh_independence\\HFM4_SG")
p = plot()
for item in readdir()
    plot_single(p, item, "P", label=item)
end
display(p)
