include("run_experiment.jl")

root = "C:\\Users\\mateo\\Files\\Investigación\\2-phase_loop\\JJI+i"
data_folder = "relap_2016_data/"
results_folder ="relap_2016_results/"
cd(root * "\\" * data_folder)

for folder in readdir()
    cd(folder)
    for subfolder in readdir()
        cd(subfolder)
        plot_name = folder * "_" * subfolder
        println(plot_name)
        experiment(
            root * "\\" * results_folder * "\\" * plot_name
            )
        cd("..")
        println()
        println()
    end
    cd("..")
end