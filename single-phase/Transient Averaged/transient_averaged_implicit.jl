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
L  = 3.0          # m                
Lh = 0.2*L        # m

# β  = 210e-6	      # K⁻¹ 
β  = 0.0	      # K⁻¹ 
T0 = 25.0		  # ºC
g  = 9.8		  # m/sꜝ
f  = 0.01         # 1
Δp = 50e3         # Pa

n  = 50
dz = L/n
z  = 0:dz:L
Q_tilde = Q/(ρ0*Cp*A*Lh)
Q_vec   = [zi < Lh ? Q_tilde : 0 for zi in z]

tf = 18.0
th = 10.0
# th = 2*tf
dt = 0.01
t  = 0:dt:tf
nt = length(t)
Δp_tilde = Δp/(ρ0*L)
Δp_vec   = [tj < th ? Δp_tilde : 0 for tj in t]

D = 1/dz*Tridiagonal(-ones(n), ones(n + 1), zeros(n))

v_old =  0.0
v  = [v_old]
T  = T0*ones(n+1)

for j in 1:nt
	A = I + dt*v_old*D
    b = T + dt*Q_vec
    b[1] = T0
    b[end] = 0
    T = A\b
	
	T_avg = integrate(z,T)/L
	v_new = v_old + dt*(Δp_vec[j] - 0.25*f*v_old*abs(v_old) - g*(1 - β*(T_avg - T0)))
	push!(v,v_new)
	v_old = v_new

	if j % 1 == 0 && t[j] > 9
		p1 = plot(z,T, title="t = $(t[j])")
		p2 = plot(t[1:j+1], v)
		p3 = plot(p1, p2, layout=(2,1))
		display(p3)
	end

end

# sleep(2)
# plot(t,v[1:end-1])
