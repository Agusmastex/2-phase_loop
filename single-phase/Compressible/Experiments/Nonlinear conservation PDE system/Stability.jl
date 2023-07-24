using Plots

T  = 100
plot(t,exp.(-t))

for k = 1:0.2:1.8
    # nt = ceil(Int,T/k)
    nt = length(0:k:T)

    y0 = 1
    y  = [y0 zeros(1,nt-1)]

    for i = 1:nt-1
        y[i+1] = y[i]*(1-k)
    end

    t = 0:k:T
    plot!(t,y',show=true)
end
