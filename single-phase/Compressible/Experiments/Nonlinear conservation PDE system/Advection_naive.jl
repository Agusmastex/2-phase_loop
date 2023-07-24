using Plots

T  = 30
L  = 10
dx = 0.1
dt = 0.05
x  = 0:dx:L
t  = 0:dt:T
N  = length(x) - 1
nt = length(t)
v = dt/dx

u  = zeros(nt,N+1)
u[1,:] = exp.(-(x.-2).^2)
u[2,:] = exp.(-(x.-2).^2)
u0 = copy(u)


# Upwind
for n in 1:nt-1
    u[n+1,1] = u[n,1] - v*(u[n,1] - u[n,end])
    for j in 2:N+1
        u[n+1,j] = u[n,j] - v*(u[n,j] - u[n,j-1])
    end
end

u_upwind = copy(u)
u = copy(u0)

# Lax-Friedrichs
for n in 1:nt-1
    u[n+1,1] = 0.5*(u[n,2] + u[n,end]) - 0.5*v*(u[n,2] - u[n,end])
    for j in 2:N
        u[n+1,j] = 0.5*(u[n,j+1] + u[n,j-1]) - 0.5*v*(u[n,j+1] - u[n,j-1])
    end
    u[n+1,end] = 0.5*(u[n,1] + u[n,end-1]) - 0.5*v*(u[n,1] - u[n,end-1])
end

u_friedrichs = copy(u)
u = copy(u0)

# Lax-Wendroff
for n in 1:nt-1
    u[n+1,1] = u[n,1] - 0.5*v*(u[n,2] - u[n,end]) + 0.5*v^2*(u[n,2] - 2*u[n,1] + u[n,end])
    for j in 2:N
        u[n+1,j] = u[n,j] - 0.5*v*(u[n,j+1] - u[n,j-1]) + 0.5*v^2*(u[n,j+1] - 2*u[n,j] + u[n,j-1])
    end
    u[n+1,end] = u[n,end] - 0.5*v*(u[n,1] - u[n,end-1]) + 0.5*v^2*(u[n,1] - 2*u[n,end] + u[n,end-1])
end

u_wendroff = copy(u)
u = copy(u0)

# Leapfrog
for n in 2:nt-1
    u[n+1,1] = u[n-1,1] - v*(u[n,2] - u[n,end])
    for j in 2:N
        u[n+1,j] = u[n-1,j] - v*(u[n,j+1] - u[n,j-1])
    end
    u[n+1,end] = u[n-1,end] - v*(u[n,1] - u[n,end-1])
end

u_leapfrog = copy(u)
u = copy(u0)

# Beam-Warming
for n in 1:nt-1
    u[n+1,1] = u[n,1] - 0.5*v*(3*u[n,1] - 4*u[n,end] + u[n,end-1]) +
                       0.5*v^2*(u[n,1] - 2*u[n,end] + u[n,end-1])

    u[n+1,2] = u[n,2] - 0.5*v*(3*u[n,2] - 4*u[n,1] + u[n,end]) +
                       0.5*v^2*(u[n,2] - 2*u[n,1] + u[n,end])
    for j in 3:N+1
        u[n+1,j] = u[n,j] - 0.5*v*(3*u[n,j] - 4*u[n,j-1] + u[n,j-2]) +
                            0.5*v^2*(u[n,j] - 2*u[n,j-1] + u[n,j-2])
    end
end

u_warming = copy(u)
u = copy(u0)

# Forward-Time Centered-Difference
for n in 2:nt-1
    u[n+1,1] = u[n,1] - 0.5*v*(u[n,2] - u[n,end])
    for j in 2:N
        u[n+1,j] = u[n-1,j] - 0.5*v*(u[n,j+1] - u[n,j-1])
    end
    u[n+1,end] = u[n,end] - 0.5*v*(u[n,1] - u[n,end-1])
end

u_FTCS = copy(u)
u = copy(u0)

u_true = copy(u0)
for n in 1:nt
    for j in 1:N+1
        u_true[n,j] = sum(exp(-(x[j]-t[n]-2+i*L)^2) for i in (0,1,2,3))
    end
end

for n in 1:nt
    p = plot(ylims=(-0.5,1))
    plot!(x,u_upwind[n,:],       label="Upwind")
    plot!(x,u_friedrichs[n,:],   label="Lax-Friedrichs")
    plot!(x,u_wendroff[n,:],     label="Lax-Wendroff")
    plot!(x,u_leapfrog[n,:],     label="Leapfrog")
    plot!(x,u_warming[n,:],      label="Beam-Warming")
    plot!(x,u_true[n,:],         label="Analytical", linewidth=3) 
    plot!(x,u_FTCS[n,:],         label="FTCS")
    display(p)
end