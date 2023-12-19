using Plots
using NumericalIntegration
using LinearAlgebra

Q  = 100.0       # J/s
ρ0 = 1000.0       # kg/m^3
Cp = 4184.0       # J/kg K
R₀ = 1.0/1000     # m
A  = π*R₀^2        # m^2
# Lh = L
Cd = 0.0001
β  = 1e-4  
T0 = 25.0 
g  = 9.8

tf = 60.0
dt = 0.1
L  = 5.0
Lh = 0.2*L        # m
dz = 0.2
nt = length(0:dt:tf)
z  = 0:dz:L
N  = length(z) - 1

D_back = 1/dz*Tridiagonal(-ones(N), ones(N + 1), zeros(N)) 				# Diferencia finita hacia atrás
D_back[1,1] = 0
D_fwd  = 1/dz*Tridiagonal(zeros(N), -ones(N + 1), ones(N)) 				# Diferencia finita hacia adelante
D_fwd[end,end] = 0

function D(vz)
	if vz >= 0
		 D_m = D_back 
	else
		 D_m = D_fwd
    end
	return D_m
end

Q_tilde  = Q/(ρ0*Cp*A*Lh)
Q_vec    = [zi < Lh ? Q_tilde : 0 for zi in z]
Q_vec[1] = 0

v_init = 0.01
T_init = T0*ones(N+1)

Δn = 1
v_save = [v_init]
T_save = [T_init]

# TAVGS  = [T0]

v_old = copy(v_init)
T_old = copy(T_init)

for n in 1:nt-1
    T_new = T_old + Q_vec - v_old*dt*D(v_old)*T_old
    T_avg = integrate(z,T_old)/L
    v_new = v_old + dt*(β*(T_avg - T0)*g - Cd/R₀*v_old*abs(v_old))

    if n % Δn == 0
        push!(v_save,v_new)
        push!(T_save,T_new)
        # push!(TAVGS, T_avg)
    end

    v_old = copy(v_new)
    T_old = copy(T_new)
end

n_save = length(v_save)
T_save = stack(T_save, dims=1)
t_save = LinRange(0,tf,n_save)

print("$((sizeof(v_save)+sizeof(T_save))/1e6) MB")

for n in 1:nt
    p1 = plot(z,T_save[n,:], title="Temperatura", xlabel="Espacio")
    p2 = plot(t_save[1:n], v_save[1:n], title="Velocidad", xlabel="Tiempo")
    # p3 = plot(t_save[1:n], TAVGS[1:n], title="Temperatura promedio")
    p  = plot(p1,p2, layout=(3,1))
    display(p)
end

T_avg_steady = T0 + Q_tilde*Lh
v_steady = 1