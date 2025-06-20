using DelimitedFiles
include("HEM3.jl")

function experiment(plot_name)
input_file = open("conditions.txt", "r")
lines = readlines(input_file)
data = [parse(Float64, split(line, "\t")[1]) for line in lines]

data[1] = data[1]*1e3   # kPa to Pa
data[4] = data[4]*1e3   # kW/m2 to W/m2

z, fields = run(data...)

α = readdlm("alpha.csv", ',')[:,2]
T = readdlm("T.csv", ',')[:,2]
P = readdlm("P.csv", ',')[:,2]

exp_fields = (α, T, P)

z_exp = [58.5, 95.3, 116.4, 136.5, 156.7]
d_inner = 19.05e-3
d_outer = 38.10e-3
Dh = d_outer - d_inner
z_exp = z_exp*Dh

plots = [plot(z, item, title=key) for (key, item) in fields]

for (plot, exp_data) in zip(plots, exp_fields)
    plot!(plot, z_exp, exp_data, marker=:circle)
end

plot(plots..., layout=(1,3), legend=false)
savefig(plot_name)
end