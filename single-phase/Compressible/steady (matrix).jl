using Plots
using LinearAlgebra
using DifferentialEquations

Q  = 1000.0         # J/s
Cp = 1996.0  		# J/kg K
R  = 5.0/100		# m
A  = π*R^2   		# m^2
L  = 10.0    	 	# m                
Lh = 0.5*L   		# m

T0 = 110.0 + 273.15	# ºC
g  = 9.8		 	# m/sꜝ
# f  = 0.05   	    # 1
f  = 3.0    	    # 1
Δp = 0.0 	       	# Pa

# Initial conditions
P0 = 1e5
v0 = 10.0			# m/s

M  = 18e-3          # kg/mol
Rg = 8.314
Cp = Cp*M           # J/mol K
Cv = Cp - Rg        # J/mol K
ρ0 = P0*M/(Rg*T0) 			# kg/m³
m  = ρ0*v0			# kg/s

# Setup
Q_hat(z) = z < Lh ? Q/(A*Lh) : 0
S = [0 0; -0.025f*m/R 0]

function F(X,p,z)
    M = [Cp/Rg*X[2] Cv/Rg*X[1]; m 1]
    b = [Q_hat(z); 0]
    return M\(b+S*X)
end

# Solve
prob = ODEProblem(F,[v0,P0],[0,L])
sol  = solve(prob, Rodas4())

# Post-process
z = 0:0.01:L
v = sol(z)[1,:]
P = sol(z)[2,:]
T = M/(m*Rg) * (P.*v)

P = 1e-5*P
T = T .- 273.15

plot_1 = plot(z,T, title="T")
plot_2 = plot(z,P, title="P")
plot_3 = plot(z,v, title="Vz")


final_plot = plot(plot_1, plot_2, plot_3, layout=(3,1))