using DelimitedFiles
using Plots
using LaTeXStrings
using TOML

models = [
    # "HEM3",
    "HFM4SZ",
    ]
databases = [
    "relap-2016",
    ]

for model_name in models
    include("models\\" * model_name * ".jl")
end

root = "C:\\Users\\mateo\\Files\\Investigación\\2-phase_loop\\thesis\\run"

function run_model(model_name, conditions, geom, n_nodes)
    model = getfield(Main, Symbol(model_name))
    fields = model.run(conditions, geom, n_nodes=n_nodes)

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

            conditions = TOML.parsefile("conditions.toml")
            geom = TOML.parsefile("geometry.toml")

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
                csv = run_model(model_name, conditions, geom, n_nodes)
                writedlm(results_name * ".csv", csv, ',')
                println()
            end
            cd(database_folder)
            println()
        end
    end
end

produce_data(["HFM4SZ"], ["relap-2016"], 40)