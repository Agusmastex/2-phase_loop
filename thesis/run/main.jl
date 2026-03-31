using DelimitedFiles
using Plots
using LaTeXStrings
using TOML

root = "C:\\Users\\mateo\\Files\\Investigación\\2-phase_loop\\thesis\\run"

for model_name in readdir(root * "\\models")
    include("models\\" * model_name)
end


function run_model(model_name, conditions, geom, n_nodes)
    model = getfield(Main, Symbol(model_name))
    fields = model.run(conditions, geom, n_nodes=n_nodes)

    header = hcat(keys(fields)...)
    matrix = hcat(values(fields)...)

    [header; matrix]
end

function produce_data(models, n_nodes)
    cd(root * "\\data")
    for scenario in readdir()
        cd(root * "\\data\\" * scenario)

        conditions = TOML.parsefile("conditions.toml")
        geom = TOML.parsefile("geometry.toml")

        for model_name in models
            results_path = root * "\\results\\"

            results_name = scenario * "_" * model_name
            println(results_name)
            csv = run_model(model_name, conditions, geom, n_nodes)
            writedlm(results_path * results_name * ".csv", csv, ',')
            println()
        end
        println()
    end
end

# produce_data(["HFM4SZ", "HFM4SZ_C"], 20)