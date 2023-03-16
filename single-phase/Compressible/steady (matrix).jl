using Plots
# using LinearAlgebra
using DifferentialEquations

Q  = 1.0  		    # J/s
Cp = 1996.0  		# J/kg K
R  = 1.0/1000		# m
A  = π*R^2   		# m^2
L  = 10.0    	 	# m                
Lh = 0.5*L   		# m

T0 = 110.0 + 273.15	# ºC
g  = 9.8		 	# m/sꜝ
# f  = 0.05   	    # 1
f  = 3.0    	    # 1
Δp = 0.0 	       	# Pa

ρ0 = 0.60 			# kg/m³
v0 = 3.0			# m/s
m  = ρ0*v0			# kg/s
P0 = 1e5
M  = 18e-3
Rg = 8.314
Cp = Cp*M           # J/mol K
Cv = Cp - Rg        # J/mol K

Q_hat(z) = z < Lh ? Q/(A*Lh) : 0

S = [0 0; -0.025f*m/R 0]

function F(X,p,z)
    M = [Cp/Rg*X[2] Cv/Rg*X[1]; m 1]
    b = [Q_hat(z); 0]
    return M/(b+S*X)
end

prob = ODEProblem(F,[v0,P0],[0,L])
sol  = solve(prob)
plot(sol, layout=(2,1))