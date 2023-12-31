using Plots
using SparseArrays
using LinearAlgebra
using CoolProp

tf = 1.0
dt = 0.01
dz = 0.05
  
# Change this
  f = 0.5
  Qval = 100e3
  g = 0
  Δn = 1
  T_adjust = 160 + 273.15

# Grid
  L  = 0.5
  z  = 0:dz:L
  t  = 0:dt:tf
  nt = length(t) - 1
  N  = length(z) - 1

# Constants
  D = 0.05
  Lh = 0.5*L
  A_flow = 0.25*π*D^2
  ρ1 = 1000
  # hfg = PropsSI("H","P",1e5,"Q",1,"Water") - PropsSI("H","P",1e5,"Q",0,"Water")
  hfg = 2e6
  R_gas = 8.314/0.018
  Cp_liq = 4184.0
  h_liqsat = PropsSI("H", "P", 1e5, "Q", 0, "Water") 
  h_critboil = 0.8*h_liqsat # Net void critical enthalpy

# Nonlinear root function
  centered_D = spdiagm(
    -1 => -1*ones(N),
     1 =>    ones(N)
  )
  centered_D[1,:] = zeros(N+1)
  centered_D[end,:] = [zeros(N-2); 2; -8; 6]
  centered_D = 0.5/dz*dropzeros(centered_D)

  upwind_D = spdiagm(
    -1 => -ones(N),
     0 =>  ones(N+1)
  )
  upwind_D[1,:] = zeros(N+1)
  upwind_D = 1/dz*dropzeros(upwind_D)

  function safe(vector)
    return [0; vector[2:end]]
  end

  p1 = 90e3
  p2 = 110e3
  rho1 = PropsSI("D", "P", p1, "Q", 1, "Water")
  rho2 = PropsSI("D", "P", p2, "Q", 1, "Water")
  function ρ_vapsat(P)
    return (rho2 - rho1)/(p2 - p1)*(P - p1) + rho1
  end
  U1 = PropsSI("U", "P", p1, "Q", 1, "Water")
  U2 = PropsSI("U", "P", p2, "Q", 1, "Water")
  function U_vapsat(P)
    return (U2 - U1)/(p2 - p1)*(P - p1) + U1
  end
  Tsat1 = PropsSI("T", "P", p1, "Q", 1, "Water")
  Tsat2 = PropsSI("T", "P", p2, "Q", 1, "Water")
  function T_vapsat(P)
    return (Tsat2 - Tsat1)/(p2 - p1)*(P - p1) + Tsat1
  end
    
  T_triple = 273.16
  function U_liq(T)
    return Cp_liq*(T - T_triple)
  end

  function boiling(T,P,α)
    Γ_w = zeros(N+1)
    for i in 1:N+1
      if U_liq(T[i]) > h_liqsat
        Γ_w[i] = q[i]/hfg
      elseif U_liq(T[i]) > h_critboil
        Γ_w[i] = q[i]*(U_liq(T[i]) - h_critboil)/(h_liqsat - h_critboil)/hfg
      end
    end
    H_interface = 5e5*α
    T_vapor = T_vapsat.(P)
    Γ_ig = H_interface.*(T - T_vapor)/hfg
    # Γ_w = zeros(N+1)
    return Γ_ig + Γ_w
    # return [ U_liq(T[i]) > h_liqsat ? q[i]/hfg : 0 for i in 1:N+1 ]
    # return [ T[i] > 100 + 273.15 ? q[i]/hfg : 0 for i in 1:N+1 ]
  end

  function F(Qk, Qn)
    p = Qk[1:N+1]
    α = Qk[N+2:2N+2]
    vm = Qk[2N+3:3N+3]
    T = Qk[3N+4:4N+4]

    p_n = Qn[1:N+1]
    α_n = Qn[N+2:2N+2]
    vm_n = Qn[2N+3:3N+3]
    T_n = Qn[3N+4:4N+4]

    ρ2   = ρ_vapsat.(p)
    ρ2_n = ρ_vapsat.(p_n)

    ρm   = (1 .- α)*ρ1 + α.*ρ2
    ρm_n = (1 .- α_n)*ρ1 + α_n.*ρ2_n

    U2   = U_vapsat.(p)
    U2_n = U_vapsat.(p_n)

    U1   = U_liq.(T)
    U1_n = U_liq.(T_n)

    ρmUm   = (1 .- α)*ρ1.*U1 + α.*ρ2.*U2
    ρmUm_n = (1 .- α_n)*ρ1.*U1_n + α_n.*ρ2_n.*U2_n

    Γ = boiling(T,p,α)

    return [
      (ρm     - ρm_n      )/dt + upwind_D*(ρm.*vm);
      (α.*ρ2  - α_n.*ρ2_n )/dt + upwind_D*(α.*ρ2.*vm) - safe(Γ);
      (vm     - vm_n      )/dt + vm.*upwind_D*vm + (1 ./ ρm).*upwind_D*p + safe(0.5*f/D*vm.*vm .+ g);
      (ρmUm - ρmUm_n)/dt + upwind_D*(ρmUm.*vm) + safe(-q + p.*centered_D*vm)
    ]
  end

# Jacobian
  e(i) = [ k == i ? 1 : 0 for k in 1:4N+4 ]
  δ = 1e-4
  function Jac(Q, Qn)
    J = zeros(4N+4,4N+4)
    for j in 1:4N+4
      J[:,j] = (F(Q + δ*e(j), Qn) - F(Q, Qn))/δ
    end
    return J
  end

# Initialize fields
  p₀ = 1e5
  α₀ = 0.01
  vm₀ = 1.0
  T₀ = 70 + 273.15
  
  q = [zi <= Lh ? Qval/(Lh*A_flow) : 0 for zi in z]

  p_initial  = p₀*ones(N+1)
  α_initial  = α₀*ones(N+1)
  vm_initial = vm₀*ones(N+1)
  T_initial  = T₀*ones(N+1)

  Qn = [p_initial; α_initial; vm_initial; T_initial]

  p_save  = [p_initial]
  α_save  = [α_initial]
  vm_save = [vm_initial]
  T_save  = [T_initial]
  t_save  = [0.0]

# Main loop
  tol = 1e-6
  Qk = copy(Qn)
  for n in 1:nt+1
      global Qk, Qn
      res = 1.0
      k = 0
      while res > tol
          k = k + 1
          x = Jac(Qk,Qn)\-F(Qk,Qn)
          Qk = Qk + x
          res = norm(x)
          print(k)
      end
      println()
      if n % Δn == 0
          p = Qk[1:N+1]
          α  = Qk[N+2:2N+2]
          vm = Qk[2N+3:3N+3]
          T = Qk[3N+4:4N+4]
          push!(p_save, p)
          push!(α_save, α)
          push!(vm_save, vm)
          push!(T_save, T)
          push!(t_save, t[n])
          println()
          println("$(round(n/(nt+1)*100, digits=1))% ")
      end
      Qn = copy(Qk)
  end

# Calculate derived fields
  n_save  = length(p_save)
  T_save  = [T_save[j] .- 273.15 for j in 1:n_save]

  field_dict = Dict(
    "p"  => p_save/1e3,
    "α"  => α_save,
    "vm" => vm_save,
    "T"  => T_save,

  )

# Plotting

  select = ["p","α", "vm","T"]
 
  ## Dynamical limits
  for n in 1:n_save
    global p
    plots = [plot(z,field_dict[name][n], title=name, formatter=:plain) for name in select]
    xlabel!("t = $(t_save[n])")
    p = plot(plots...)
    display(p)
  end

  ## Static limits
  # field_min = Dict()
  # field_max = Dict()
  
  # for item in field_dict
  #   name, field_save = item
  #   matrix = hcat(field_save...)
  #   field_min[name] = minimum(matrix)
  #   field_max[name] = maximum(matrix)
  # end
  
  # for n in 1:n_save
  #   global p 
  #   plots = [plot(z,field_dict[name][n], ylims=(0.995*field_min[name], 1.005*field_max[name]), title=name, formatter=:plain) for name in select]
  #   xlabel!("t = $(t_save[n])")
  #   p = plot(plots...)
  #   display(p)
  #   # sleep(1)
  # end

  # simulation = "Backward Upwind Nonconserved NumJac [p α vm]"
  # p[:plot_title] = simulation
  # plot(p)
  # savefig(simulation)