using Plots
using SparseArrays
using LinearAlgebra

tf = 5.0
dt = 0.02
dz = 0.01
  
# Change this
  f = 0.1
  Qval = 1000
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
  hfg = 2.2e6
  R_gas = 8.314/(18/1000)

# Nonlinear root function
  centered_D = spdiagm(
    -1 => -1*ones(N),
     1 =>    ones(N)
  )
  centered_D[1,:] = zeros(N+1)
  centered_D[end,:] = [zeros(N-1); -2; 2]
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

  function F(Qk, Qn)
    ρm = Qk[1:N+1]
    α = Qk[N+2:2N+2]
    vm = Qk[2N+3:3N+3]
    ρm_n = Qn[1:N+1]
    α_n = Qn[N+2:2N+2]
    vm_n = Qn[2N+3:3N+3]
    ρ2 = (ρm - (1 .- α)*ρ1)./α
    p = R_gas*T_adjust*ρ2
    return [
      (ρm - ρm_n)/dt + upwind_D*(ρm.*vm);
      (α  - α_n )/dt - upwind_D*((1 .- α).*vm) - safe(Γ/ρ1);
      (vm - vm_n)/dt + vm.*upwind_D*vm + upwind_D*p./ρm + safe(0.5*f/D*vm.*vm .+ g)
    ]
  end

# Jacobian
  e(i) = [ k == i ? 1 : 0 for k in 1:3N+3 ]
  δ = 1e-4
  function Jac(Q, Qn)
    J = zeros(3N+3,3N+3)
    for j in 1:3N+3
      J[:,j] = (F(Q + δ*e(j), Qn) - F(Q, Qn))/δ
    end
    return J
  end

# Initialize fields
  ρ2 = 0.5
  α₀ = 0.1
  ρm₀ = (1-α₀)*ρ1 + α₀*ρ2
  vm₀ = 1.0
  
  q = [zi <= Lh ? Qval/(Lh*A_flow) : 0 for zi in z]
  Γ = q/hfg

  ρm_initial = ρm₀*ones(N+1)
  α_initial  = α₀*ones(N+1)
  vm_initial = vm₀*ones(N+1)

  Qn = [ρm_initial; α_initial; vm_initial]

  ρm_save = [ρm_initial]
  α_save  = [α_initial]
  vm_save = [vm_initial]
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
          ρm = Qk[1:N+1]
          α  = Qk[N+2:2N+2]
          vm = Qk[2N+3:3N+3]
          push!(ρm_save, ρm)
          push!(α_save, α)
          push!(vm_save, vm)
          push!(t_save, t[n])
          println()
          println("$(round(n/(nt+1)*100, digits=1))% ")
      end
      Qn = copy(Qk)
  end

# Calculate derived fields
  n_save  = length(ρm_save)
  ρ2_save = [(ρm_save[j] - (1 .- α_save[j])*ρ1)./α_save[j] for j in 1:n_save]
  p_save  = [R_gas*T_adjust*ρ2_save[j]/1e2 for j in 1:n_save]
  Gm_save = [ρm_save[j].*vm_save[j] for j in 1:n_save]

  field_dict = Dict(
    "ρm" => ρm_save,
    "α"  => α_save,
    "vm" => vm_save,

    "ρ2" => ρ2_save,
    "p"  => p_save,
    "Gm" => Gm_save
  )

# Plotting

  select = ["ρm","α", "vm","p"]
 
  ## Dynamical limits
  for n in 1:n_save
    plots = [plot(z,field_dict[name][n], title=name, formatter=:plain) for name in select]
    xlabel!("t = $(t_save[n])")
    p = plot(plots...)#, layout=(length(select), 1))
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
  #   p = plot(plots..., layout=(length(select), 1))
  #   display(p)
  #   # sleep(1)
  # end

  # simulation = "Backward Upwind Primitive NumJac"
  # p[:plot_title] = simulation
  # plot(p)
  # savefig(simulation)