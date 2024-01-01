using Plots
using SparseArrays
using LinearAlgebra

tf = 0.7
dt = 0.02
dz = 0.01
  
# Change this
  f = 0.5
  Qval = 100
  g = 0
  Δn = 2

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

  big_D = sparse(
    [upwind_D zeros(N+1,N+1) zeros(N+1,N+1);
     zeros(N+1,N+1) upwind_D zeros(N+1,N+1);
     zeros(N+1,N+1) zeros(N+1,N+1) upwind_D])

  function safe(vector)
    return [0; vector[2:end]]
  end

  function convection(Q)
    ρ = Q[1:N+1]
    G = Q[N+2:2N+2]
    e = Q[2N+3:3N+3]
    return [
      G;
      G.^2 ./ ρ + (γ-1)*e;
      e.*G./ρ
    ]
  end

  function source(Q)
    ρ = Q[1:N+1]  
    G = Q[N+2:2N+2]
    e = Q[2N+3:3N+3]
    v = G./ρ
    return [
      zeros(N+1);
      safe(0.5*f/D*(G.^2 ./ ρ) + ρ*g);
      safe(-q + (γ-1)*e.*centered_D*v)
    ]
  end

  F(Qk, Qn) = (Qk - Qn)/dt + big_D*convection(Qk) + source(Qk)

# Jacobian
  J1 = 1/dt*sparse(I,N+1,N+1)
  J2 = upwind_D
  J3 = zeros(N+1,N+1)
  J6 = (γ-1)*upwind_D

  function Jac(Q)
      ρ = Q[1:N+1]
      G = Q[N+2:2N+2]
      e = Q[2N+3:3N+3]

      # J4
      lower_diagonal_4 = [   1/dz*(G[i-1]/ρ[i-1])^2            for i in 2:N+1]
      diagonal_4       = [ -(1/dz + 0.5*f/D)*(G[i]/ρ[i])^2 + g for i in 1:N+1]
      J4 = spdiagm(
        -1 => lower_diagonal_4,
         0 => diagonal_4)
      J4[1,:] = zeros(N+1)

      # J5
      lower_diagonal_5 = [ -2/dz*G[i-1]/ρ[i-1]           for i in 2:N+1]
      diagonal_5       = [ 1/dt + (2/dz + f/D)*G[i]/ρ[i] for i in 1:N+1]

      J5 = spdiagm(
        -1 => lower_diagonal_5,
         0 => diagonal_5)
      J5[1,:] = [1/dt; zeros(N)]

      # J7
      lower_diagonal_7 = [  1/dz*G[i-1]/ρ[i-1]^2 * (e[i-1] + 0.5*(γ-1)*e[i]) for i in 2:N+1]
      diagonal_7       = [ -1/dz*e[i]*G[i]/ρ[i]^2                            for i in 1:N+1]
      upper_diagonal_7 = [ -0.5*(γ-1)/dz*G[i+1]/ρ[i+1]^2                     for i in 1:N]

      J7 = spdiagm(
        -1 => lower_diagonal_7,
         0 => diagonal_7,
         1 => upper_diagonal_7)
      J7[1,:]     = zeros(N+1)
      J7[N+1,N]   = 1/dz*G[N]/ρ[N]^2 * (e[N] + (γ-1)*e[N+1])
      J7[N+1,N+1] = -γ/dz * e[N+1]*G[N+1]/ρ[N+1]^2

      # J8
      lower_diagonal_8 = [ -1/dz/ρ[i-1] * (e[i-1] + 0.5*(γ-1)*e[i]) for i in 2:N+1]
      diagonal_8       = [ 1/dz*e[i]/ρ[i]                           for i in 1:N+1]
      upper_diagonal_8 = [ 0.5*(γ-1)/dz*e[i]/ρ[i+1]                 for i in 1:N]

      J8 = spdiagm(
        -1 => lower_diagonal_8,
         0 => diagonal_8,
         1 => upper_diagonal_8)

      J8[1,:]     = zeros(N+1)
      J8[N+1,N]   = -1/dz/ρ[N]*(e[N] + (γ-1)*e[N+1])
      J8[N+1,N+1] = γ/dz*e[N+1]/ρ[N+1]

      # J9
      lower_diagonal_9 = [ -1/dz*G[i-1]/ρ[i-1] for i in 2:N+1]
      diagonal_9       = [0; [1/dt + 1/dz*G[i]/ρ[i] + 0.5*(γ-1)/dz*(G[i+1]/ρ[i+1] - G[i-1]/ρ[i-1]) for i in 2:N]; 0]

      J9 = spdiagm(
        -1 => lower_diagonal_9,
         0 => diagonal_9)

      J9[1,:]     = [1/dt; zeros(N)]
      J9[N+1,N]   = -1/dz * G[N]/ρ[N]
      J9[N+1,N+1] = 1/dt + γ/dz*G[N+1]/ρ[N+1] - (γ-1)/dz * G[N]/ρ[N]


      return dropzeros([J1 J2 J3; J4 J5 J6; J7 J8 J9])
  end

# Initialize fields
  ρ0 = 0.5
  v0 = 1.0
  T0 = 110 + 273.15
  U0 = Cv*T0
  
  q = [zi <= Lh ? Qval/(Lh*A_flow) : 0 for zi in z]

  ρ_initial = ρ0*   ones(N+1)
  G_initial = ρ0*v0*ones(N+1)
  e_initial = ρ0*U0*ones(N+1)

  Qn = [ρ_initial; G_initial; e_initial]

  ρ_save = [ρ_initial]
  G_save = [G_initial]
  e_save = [e_initial]
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
          x = Jac(Qk)\-F(Qk,Qn)
          Qk = Qk + x
          res = norm(x)
          print(k)
      end
      println()
      if n % Δn == 0
          ρ = Qk[1:N+1]
          G = Qk[N+2:2N+2]
          e = Qk[2N+3:3N+3]
          push!(ρ_save, ρ)
          push!(G_save, G)
          push!(e_save, e)
          push!(t_save, t[n])
          println()
          println("$(round(n/(nt+1)*100, digits=1))% ")
      end
      Qn = copy(Qk)
  end

# Calculate primitive fields
  n_save = length(ρ_save)
  v_save = [G_save[j]./ρ_save[j]   for j in 1:n_save]
  U_save = [e_save[j]./ρ_save[j]   for j in 1:n_save]
  p_save = [(γ-1)*e_save[j]/1e2    for j in 1:n_save]
  T_save = [U_save[j]/Cv .- 273.15 for j in 1:n_save]

  field_dict = Dict(
    "ρ" => ρ_save,
    "G" => G_save,
    "e" => e_save,

    "p" => p_save,
    "v" => v_save,
    "T" => T_save,

    "U" => U_save
  )

# Plotting

  select = ["ρ","v","T","p"]
  
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
    p = plot(plots...)#, layout=(length(select), 1))
    display(p)
  end

  # simulation = "Backward Upwind Conserved AnJac"
  # p[:plot_title] = simulation
  # plot(p)
  # savefig(simulation)