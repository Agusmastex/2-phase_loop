# Calculates thermal expansion coefficient β using the same scheme as Leon Matos
# Data from NIST

# Get data from textfile
using Plots

file = open("fluid.txt")
lines = readlines(file)
close(file)

table = []
for line in lines
		push!(table,split(line,"\t"))
end

density = []
temperature = []

for element in table
		push!(density,element[3])
		push!(temperature,element[1])
end

popfirst!(density)
popfirst!(temperature)

ρ = parse.(Float64,density[1:end-2])
T = parse.(Float64,temperature[1:end-2])


# Calculate β

β = []
for i in 1:length(ρ)-1
		beta = -(ρ[i+1]/ρ[i] - 1)/0.1666
		push!(β,beta)
end

h = plot(T[2:end],β)
display(h)

# Alternative β

β1 = diff(ρ)./diff(T)
β1 = -β1./ρ[2:end]

h = plot(T[2:end],β1)
display(h)
