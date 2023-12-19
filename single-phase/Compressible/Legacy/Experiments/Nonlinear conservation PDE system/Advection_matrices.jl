using Plots
using LinearAlgebra

function circulant(a)
    n = length(a)
    M = zeros(n,n)
    for j in 1:n
        for i in 1:n
            if i >= j
                M[i,j] = a[i-j+1]
            else
                M[i,j] = a[n+i-j+1]
            end
        end
    end
    return M
end

T  = 30
L  = 10
dx = 0.1
dt = 0.05
x  = 0:dx:L
t  = 0:dt:T
N  = length(x) - 1
nt = length(t)
v  = dt/dx

D_center1   = circulant([0 -1 zeros(1,N-2) 1])
D_center2   = circulant([-2 1 zeros(1,N-2) 1])
D_upwind    = circulant([1 -1 zeros(1,N-1)])
D_backward1 = circulant([3 -4 1 zeros(1,N-2)])
D_backward2 = circulant([1 -2 1 zeros(1,N-2)])

u = zeros(nt,N+1)
u[1,:] = exp.(-(x.-2).^2)
u0 = copy(u)

# Upwind (explicit)
for n in 1:nt-1
    u[n+1,:] = u[n,:] - v*D_upwind*u[n,:]
end
u_upwind_e = copy(u)
u = copy(u0)

# Upwind (implicit)
A = I + v*D_upwind
for n in 1:nt-1
    b = u[n,:]
    u[n+1,:] = A\b
end
u_upwind_i = copy(u)
u = copy(u0)

# Lax-Wendroff
for n in 1:nt-1
    u[n+1,:] = u[n,:] - 0.5*v*D_center1*u[n,:] + 0.5*v^2*D_center2*u[n,:]
end
u_wendroff = copy(u)
u = copy(u0)

# Beam-Warming
for n in 1:nt-1
    u[n+1,:] = u[n,:] - 0.5*v*D_backward1*u[n,:] + 0.5*v^2*D_backward2*u[n,:]
end
u_warming = copy(u)
u = copy(u0)

u_true = copy(u0)
for n in 1:nt
    for j in 1:N+1
        u_true[n,j] = sum(exp(-(x[j]-t[n]-2+i*L)^2) for i in (0,1,2,3))
    end
end



for n in 1:nt
    p = plot(ylims=(-0.2,1))
    plot!(x,u_upwind_e[n,:],      label="Upwind explicit")
    plot!(x,u_upwind_i[n,:],      label="Upwind implicit")
    # plot!(x,u_wendroff[n,:],    label="Lax-Wendroff")
    # plot!(x,u_warming[n,:],     label="Beam-Warming")
    # plot!(x,u_true[n,:],        label="Analytical", linewidth=3) 
    display(p)
end