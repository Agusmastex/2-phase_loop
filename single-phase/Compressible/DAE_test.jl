using Plots
# using LinearAlgebra
using DifferentialEquations

Q  = 1.0  		# J/s
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
p0 = 1e5
M  = 18e-3
Rg = 8.314
Cv = Cp - Rg/M

Q_hat(z) = z < Lh ? Q/(A*Lh) : 0

function system(out,du,u,p,z)
    global
    v,T,P = u
    out[1] = - du[3] - 0.025*f*m*v - m*du[1]
    out[2] = Q_hat(z) - P*du[1] - m*Cv*du[2]
    out[3] = P*M - m/v * Rg*T
end

u0  = [v0, T0, p0]
du0 = [0.1,0.1,0.1]
zspan = [0,L]

prob = DAEProblem(system, du0, u0, zspan, differential_vars=[true,true,false])
sol  = solve(prob)
plot(sol, layout=(3,1))

# function f2(out,du,u,p,t)
#   out[1] = - 0.04u[1]              + 1e4*u[2]*u[3] - du[1]
#   out[2] = + 0.04u[1] - 3e7*u[2]^2 - 1e4*u[2]*u[3] - du[2]
#   out[3] = u[1] + u[2] + u[3] - 1.0
# end

# u₀ = [1.0, 0, 0]
# du₀ = [-0.04, 0.04, 0.0]
# tspan = (0.0,100000.0)

# (0.0, 100000.0)

# prob = DAEProblem(f2, du₀, u₀, tspan, differential_vars=[true,true,false])

# # using Sundials
# sol = solve(prob)
# plot(sol, xscale=:log10, tspan=(1e-6, 1e5), layout=(3,1))