using Plots
using ToeplitzMatrices
using LinearAlgebra
using DifferentialEquations
using NumericalIntegration

Q  = 1000.0       # J/s
ρ0 = 1000.0       # kg/m^3
Cp = 4184.0       # J/kg K
R  = 1.0/1000     # m
A  = π*R^2        # m^2
Vz = 2.0          # m/s
L  = 10.0          # m                
Lh = 0.2*L        # m

# β  = 210e-6	      # K⁻¹ 
β  = 0.01         # K⁻¹ 
T0 = 25.0		  # ºC
g  = 9.8		  # m/sꜝ
f  = 0.01         # 1
Δp = 50e3         # Pa

# n  = 100
# dz = L/n

dz = 0.05
z  = 0:dz:L
n  = length(z) - 1

Q_tilde = Q/(ρ0*Cp*A*Lh)
Q_vec   = [zi < Lh ? Q_tilde : 0 for zi in z]

tf = 8.0
th = 0.5
# th = 2*tf
dt = 0.005
t  = 0:dt:tf
nt = length(t)
Δp_tilde = Δp/(ρ0*L)
Δp_vec   = [tj < th ? Δp_tilde : 0 for tj in t]

# D = 1/dz*Tridiagonal(-ones(n), ones(n + 1), zeros(n))

function D_f(vz)
	if vz > 0
		 D_m = 1/dz*Tridiagonal(-ones(n), ones(n + 1), zeros(n)) 				# Diferencia finita hacia atrás
	else
		 D_m = 1/dz*Tridiagonal(zeros(n), -ones(n + 1), ones(n)) 				# Diferencia finita hacia adelante
	end
	return D_m
end
v_old =  0.0
v  = [v_old]
T  = T0*ones(n+1)

Tlast = T0
Tmax = 120

for j in 1:nt-1
	global Tlast, v_old, T
	if v_old > 0
		Tlast = T[end]
	end
	T = T + dt*(Q_vec - v_old*D_f(v_old)*T)

	if v_old < 0
		T[end] = Tlast
	else
		T[1] = T0
	end
	
	T_avg = integrate(z,T)/L
	v_new = v_old + dt*(Δp_vec[j] - 0.25*f*v_old*abs(v_old) - g*(1 - β*(T_avg - T0)))
	push!(v,v_new)
	v_old = v_new

	if j % 5 == 0 # && t[j] > 10
		p1 = plot(z,T, title="t = $(t[j])", ylim = (T0 - 1, Tmax + 1))
		p2 = plot(t[1:j+1], v)
	    p3 = plot(p1, p2, layout=(2,1))
		display(p3)
	end

end

# sleep(2)
# plot(t,v[1:end-1])
