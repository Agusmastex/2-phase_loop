using Plots

# Calcula el perfil de presiones en estado estacionario y sin fricción del problema de tubería abierta
# Se empieza a notar el carácter parbólico del perfil sólo para valores grandes de β

Q  = 1000.0       # J/s
ρ0 = 1000.0       # kg/m^3
Cp = 4184.0       # J/kg K
R  = 1.0/1000     # m
A  = π*R^2        # m^2
Vz = 2.0          # m/s
Lh = 4.0          # m
L  = 10.0         # m                

T0 = 25.0         # C
f  = 0.0
β  = 0.03         # K⁻¹ 
#β  = 210e-6      # K⁻¹ 
g  = 9.8          # m/s²

T_initial = 25.0
Q_tilde   = Q/(ρ0*Cp*A*Lh)

n_nodes = 50
dz = L/n_nodes
z = -0.5*dz:dz:L+dz
T = zeros(n_nodes + 2) .+ T_initial
Q_on = [zi < Lh ? 1 : 0 for zi in z]
dt = 0.01
tf = 7.0
t  = 0:dt:tf

Tmax = T_initial + Q/(Vz*ρ0*Cp*A)	

Tnew = copy(T)
for ti in t
		Told = copy(Tnew)
		for i in 2:n_nodes+1
				Tnew[i] = Told[i] + dt*(Q_on[i]*Q_tilde - Vz*(Told[i] - Told[i-1])/dz)
		end
		Tnew[end] = Tnew[end-1]
		Tnew[1] = 2*T_initial - Tnew[2]

		#plot(z,Tnew, ylim = (T_initial-0.1,Tmax+0.1), marker=:circle, title="t = $ti", show=true)
end 
T = copy(Tnew)

p = zeros(n_nodes + 2)

for i in 2:n_nodes+1
		p[i+1] = p[i] - dz*ρ0*g*(1 - β*(T[i]-T0))
end

plot(z,p, marker=:circle)

