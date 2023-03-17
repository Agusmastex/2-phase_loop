using Plots
using LinearAlgebra
using DifferentialEquations

Q  = 1000.0         # J/s       Potencia disipada en la resistencia
Cp = 1996.0  		# J/kg K    Capacidad calorífica a presión constante específica
R₀ = 3.0/100		# m         Radio de la tubería
A  = π*R₀^2   		# m^2       Área de sección transversal de la tubería
L  = 10.0    	 	# m         Longitud de la tubería 
Lh = 0.5*L   		# m         Longitud del tramo del calentador

T0 = 110.0 + 273.15	# K         Temperatura de entrada, y de referencia para Boussinesq
g  = 9.8		 	# m/s       Aceleración de la gravedad
f  = 0.05   	    # 1         Factor de fricción de Darcy
# f  = 3.0    	    # 1

# Initial conditions
P0 = 5e5            # Pa        Presión de entrada
v0 = 10.0			# m/s       Velocidad de entrada

M  = 18e-3          # kg/mol    Masa molar del fluido
Rg = 8.314          # J/mol k   Constante de los gases ideales
Cp = Cp*M           # J/mol K   Capacidad calorífica a presión constante molar
Cv = Cp - Rg        # J/mol K   Capacidad calorífica a volumen constante molar
ρ0 = P0*M/(Rg*T0) 	# kg/m³     Densidad en la entrada
m  = ρ0*v0			# kg/m²s    Densidad de flujo de masa 

# Setup
Q_hat(z) = z < Lh ? Q/(A*Lh) : 0

function system(out,du,u,p,z)
    global
    v,T,P = u
    out[1] = du[3] + 0.025*f*m*v + m*du[1]
    out[2] = Q_hat(z) - P*du[1] - m*Cv*du[2]
    out[3] = P*M - m/v * Rg*T
end

u0  = [v0, T0, P0]
du0 = [0,0,0] # This should be unnecesary
zspan = [0,L]

prob = DAEProblem(system, du0, u0, zspan, differential_vars=[true,true,false])
sol  = solve(prob)

# Post-process
z = 0:0.01:L
v = sol(z)[1,:]
T = sol(z)[2,:]
P = sol(z)[3,:]
ρ = m./v

P = P/100000
T = T .- 273.15

variables = [T,P,v,ρ]
variable_names = ["T", "P", "Vz", "ρ"]

println("Flujo másico")
display(m*A*1000)

plots = [plot(z, quantity, title = quantity_title) for (quantity, quantity_title) in zip(variables, variable_names)]
final_plot = plot(plots..., layout=(2,2), legend=false)
