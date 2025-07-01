using DelimitedFiles
using Plots
using LaTeXStrings

models = ["HEM3", "HFM4_SG"]
databases = ["relap-2016"]
nodes = 20

for model_name in models
    include("models\\" * model_name * ".jl")
end

root = "C:\\Users\\mateo\\Files\\Investigación\\2-phase_loop\\thesis\\run"

d_inner = 19.05e-3
d_outer = 38.10e-3
Dh = d_outer - d_inner

function produce_data(models, databases)
    for database in databases
        database_folder = root * "\\data\\" * database
        cd(database_folder)
        for scenario in readdir()
            cd(scenario)
            input_file = open("conditions.txt", "r")
            lines = readlines(input_file)
            conditions = [parse(Float64, split(line, "\t")[1]) for line in lines]
            conditions[1] = conditions[1]*1e3   # kPa to Pa
            conditions[4] = conditions[4]*1e3   # kW/m2 to W/m2

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

                model = getfield(Main, Symbol(model_name))
                fields = model.run(conditions..., n_nodes=nodes)

                header = hcat(keys(fields)...)
                matrix = hcat(values(fields)...)

                writedlm(results_name * ".csv", [header; matrix], ',')
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
        α_plot = plot(title = L"\alpha")
        T_plot = plot(title = L"T")
        p_plot = plot(title = L"P")

        cd(database_folder * "\\" *scenario)

        (z_α, α) = eachcol(readdlm("alpha.csv", ','))
        (z_T, T) = eachcol(readdlm("T.csv", ','))
        (z_p, p) = eachcol(readdlm("P.csv", ','))

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
        print()
    end

end

# produce_data(models, databases)
plot_data("relap-2016")