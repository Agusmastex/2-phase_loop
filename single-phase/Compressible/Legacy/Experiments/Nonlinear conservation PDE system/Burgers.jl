using Plots

T  = 10
L  = 10
dx = 0.1
dt = 0.05
x  = 0:dx:L
t  = 0:dt:T
N  = length(x) - 1
nt = length(t)
v  = dt/dx

u = zeros(nt,N+1)
u[1,:] = exp.(-(x.-2).^2)
u[1,:] = sin.(-x).^2
u[1,:] = 0.05*x
u0 = copy(u)

# Upwind
for n in 1:nt-1
    u[n+1,1] = u[n,1] - 0.5*v*(u[n,1]^2 - u[n,end]^2)
    for j in 2:N+1
        u[n+1,j] = u[n,j] - 0.5*v*(u[n,j]^2 - u[n,j-1]^2)
    end
end
u_upwind = copy(u)
u = copy(u0)

# Lax-Wendroff
for n in 1:nt-1
    u[n+1,1] = u[n,1] - 0.25*v*(u[n,2]^2 - u[n,end]^2) + 0.25*v^2*(u[n,2]^2 - 2*u[n,1]^2 + u[n,end]^2)
    for j in 2:N
        u[n+1,j] = u[n,j] - 0.25*v*(u[n,j+1]^2 - u[n,j-1]^2) + 0.25*v^2*(u[n,j+1]^2 - 2*u[n,j]^2 + u[n,j-1]^2)
    end
    u[n+1,end] = u[n,end] - 0.25*v*(u[n,1]^2 - u[n,end-1]^2) + 0.25*v^2*(u[n,1]^2 - 2*u[n,end]^2 + u[n,end-1]^2)
end

u_wendroff = copy(u)
u = copy(u0)

# Beam-Warming
for n in 1:nt-1
    u[n+1,1] = u[n,1] - 0.25*v*(3*u[n,1]^2 - 4*u[n,end]^2 + u[n,end-1]^2) +
                       0.25*v^2*(u[n,1]^2 - 2*u[n,end]^2 + u[n,end-1]^2)

    u[n+1,2] = u[n,2] - 0.25*v*(3*u[n,2]^2 - 4*u[n,1]^2 + u[n,end]^2) +
                       0.25*v^2*(u[n,2]^2 - 2*u[n,1]^2 + u[n,end]^2)
    for j in 3:N+1
        u[n+1,j] = u[n,j] - 0.25*v*(3*u[n,j]^2 - 4*u[n,j-1]^2 + u[n,j-2]^2) +
                            0.25*v^2*(u[n,j]^2 - 2*u[n,j-1]^2 + u[n,j-2]^2)
    end
end
u_warming = copy(u)
u = copy(u0)

for n in 1:nt
    p = plot(ylims=(-1,1))
    plot!(x,u_upwind[n,:],      label="Upwind")
    plot!(x,u_wendroff[n,:],    label="Lax-Wendroff")
    plot!(x,u_warming[n,:],     label="Beam-Warming")
    display(p)
end