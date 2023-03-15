using Plots

f(x) = exp(-x^2)
x = LinRange(-10,10,100)
println("Yesh")
h = plot(x,f.(x))
display(h)
println("Non")


