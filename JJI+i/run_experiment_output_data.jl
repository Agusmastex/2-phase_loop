using DelimitedFiles
include("HEM3.jl")

function experiment(data_file_name)
input_file = open("conditions.txt", "r")
lines = readlines(input_file)
data = [parse(Float64, split(line, "\t")[1]) for line in lines]

data[1] = data[1]*1e3   # kPa to Pa
data[4] = data[4]*1e3   # kW/m2 to W/m2

z, fields = run(data...)

A = [z fields[1][2] fields[2][2] fields[3][2]]

writedlm(data_file_name * ".csv", A, ',')

end