using Plots
using SparseArrays
using LinearAlgebra

tf = 1.0
dt = 0.02
dz = 0.01
  
# Change this
  f = 0
  Qval = 100
  g = 0
  Δn = 1

# Grid
  L  = 0.5
  z  = 0:dz:L
  t  = 0:dt:tf
  nt = length(t) - 1
  N  = length(z) - 1

# Constants
  γ = 1.4
  D = 0.05
  Lh = 0.5*L
  A_flow = 0.25*π*D^2
  Cv = 1.3e3

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
    ρ = Qk[1:N+1]
    v = Qk[N+2:2N+2]
    U = Qk[2N+3:3N+3]
    ρn = Qn[1:N+1]
    vn = Qn[N+2:2N+2]
    Un = Qn[2N+3:3N+3]
    p = (γ-1)*ρ.*U
    return [
      (ρ    - ρn    )/dt + upwind_D*(ρ.*v);
      (ρ.*v - ρn.*vn)/dt + upwind_D*(ρ.*v.*v) + centered_D*p + safe(0.5*f/D*ρ.*v.*v + ρ*g);
      (ρ.*U - ρn.*Un)/dt + upwind_D*(ρ.*v.*U) - safe(q) + p.*centered_D*v
    ]
  end

# Jacobian
  e(i) = [ k == i ? 1 : 0 for k in 1:3N+3 ]

  function Jac(Q, Qn)
    δ = 1e-4
    J = zeros(3N+3,3N+3)
    for j in 1:3N+3
      J[:,j] = (F(Q + δ*e(j), Qn) - F(Q, Qn))/δ
    end
    return J
  end

# Initialize fields
  ρ0 = 0.5
  v0 = 1.0
  T0 = 110 + 273.15
  U0 = Cv*T0
  
  q = [zi <= Lh ? Qval/(Lh*A_flow) : 0 for zi in z]

  ρ_initial = ρ0*ones(N+1)
  v_initial = v0*ones(N+1)
  U_initial = U0*ones(N+1)

  Qn = [ρ_initial; v_initial; U_initial]

  ρ_save = [ρ_initial]
  v_save = [v_initial]
  U_save = [U_initial]
  t_save = [0.0]

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
          ρ = Qk[1:N+1]
          v = Qk[N+2:2N+2]
          U = Qk[2N+3:3N+3]
          push!(ρ_save, ρ)
          push!(v_save, v)
          push!(U_save, U)
          push!(t_save, t[n])
          println()
          println("$(round(n/(nt+1)*100, digits=1))% ")
      end
      Qn = copy(Qk)
  end

# Calculate derived fields
  n_save = length(ρ_save)
  p_save = [(γ-1)*ρ_save[j].*U_save[j]/1e2    for j in 1:n_save]
  T_save = [U_save[j]/Cv .- 273.15 for j in 1:n_save]

  field_dict = Dict(
    "ρ" => ρ_save,
    "v" => v_save,
    "U" => U_save,

    "p" => p_save,
    "T" => T_save
  )

# Plotting

  select = ["p","v","T"]
  
  ## Dynamical limits
  # for n in 1:n_save
  #   plots = [plot(z,field_dict[name][n], title=name, formatter=:plain) for name in select]
  #   xlabel!("t = $(t_save[n])")
  #   p = plot(plots..., layout=(length(select), 1))
  #   display(p)
  # end

  ## Static limits
  field_min = Dict()
  field_max = Dict()
  
  for item in field_dict
    name, field_save = item
    matrix = hcat(field_save...)
    field_min[name] = minimum(matrix)
    field_max[name] = maximum(matrix)
  end
  
  for n in 1:n_save
    global p 
    plots = [plot(z,field_dict[name][n], ylims=(field_min[name], field_max[name]), title=name, formatter=:plain) for name in select]
    xlabel!("t = $(t_save[n])")
    p = plot(plots..., layout=(length(select), 1))
    display(p)
    sleep(0.1)
  end

  # simulation = "Backward Upwind Primitive NumJac"
  # p[:plot_title] = simulation
  # plot(p)
  # savefig(simulation)