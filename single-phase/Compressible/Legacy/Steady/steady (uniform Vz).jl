using Plots
# using LinearAlgebra
using DifferentialEquations

Q  = 1.0    		# J/s
Cp = 1996.0  		# J/kg K
R  = 1.0/1000		# m
A  = π*R^2   		# m^2
L  = 10.0    	 	# m                
Lh = 0.5*L   		# m

T0 = 100.0		  	# ºC
g  = 9.8		 	# m/sꜝ
# f  = 0.05   	    # 1
f  = 3.0    	    # 1
Δp = 0.0 	       	# Pa

ρ0 = 0.60 			# kg/m³
v0 = 3.0			# m/s
m  = ρ0*v0			# kg/s
dpdz = -Δp/L

Q_hat(z) = z < Lh ? Q/(m*Cp*A*Lh) : 0

# H[1]: Velocity
# H[2]: Temperature

F(H,p,z) = Q_hat(z) + v0/(m*Cp)*dpdz

prob = ODEProblem(F,T0,[0,L])
sol = solve(prob, Rodas4())

z = 0:0.01:L
T = sol(z)

plot(z,T)
