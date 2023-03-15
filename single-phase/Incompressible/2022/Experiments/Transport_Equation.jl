using Plots

# Resuelve la ecuación de transporte u_t + c*u_x = 0
# Con una condición inicial u(x,0) = f(x) donde f es una función definida por partes con forma triangular

Q  = 1000.0       # J/s
ρ0 = 1000.0       # kg/m^3
Cp = 4184.0       # J/kg K
R  = 1.0/1000     # m
A  = π*R^2        # m^2
Vz = 2.0          # m/s
Lh = 4.0          # m
L  = 10.0         # m                

T_initial = 25.0
Q_tilde   = Q/(ρ0*Cp*A*Lh)

n_nodes = 100
dz = L/n_nodes
z = -0.5*dz:dz:L+dz
# T = zeros(n_nodes + 2) .+ T_initial
# Q_on = [zi < Lh ? 1 : 0 for zi in z]
T = zeros(n_nodes + 2)

h0 = 2
a  = 2
b  = 1

for i in 1:n_nodes+2
		if z[i] < h0 - b
				T[i] = T_initial

		elseif z[i] < h0
				T[i] = a/b * (z[i] - h0) + h0 + T_initial

		elseif z[i] < h0 + b
				T[i] = -a/b * (z[i] - h0) + h0 + T_initial

		elseif z[i] > h0 + b
				T[i] = T_initial
		end
end

# plot(z,T, marker=:circle)
# readline()

Q_on = [0 for zi in z] # Término de source nulo
dt = 0.01
tf = 7.0
t  = 0:dt:tf

Tmax = T_initial + a

Tnew = copy(T)
for ti in t
		Told = copy(Tnew)
		for i in 2:n_nodes+1
				Tnew[i] = Told[i] + dt*(Q_on[i]*Q_tilde - Vz*(Told[i] - Told[i-1])/dz)
		end
		Tnew[end] = Tnew[end-1]
		Tnew[1] = 2*T_initial - Tnew[2]

		plot(z,Tnew, ylim = (T_initial-0.1,Tmax+0.1), title="t = $ti", show=true)

		#plot(z,Tnew, marker=:circle, title="t = $ti", show=true)
end 


