using Plots
using SparseArrays
using LinearAlgebra

tf = 1.0
dt = 0.01
dz = 0.001
  
# Change this
  f = 0
  Qval = 0.1
  g = 0
  Δn = 1

# Grid
  L  = 0.2
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
  little_D = spdiagm(
    -1 => -1*ones(N),
     1 =>    ones(N)
  )
  little_D[1,:] = zeros(N+1)
  little_D[end,:] = [zeros(N-1); -2; 2]
  little_D = 0.5/dz*dropzeros(little_D)

  big_D = sparse(
    [little_D zeros(N+1,N+1) zeros(N+1,N+1);
     zeros(N+1,N+1) little_D zeros(N+1,N+1);
     zeros(N+1,N+1) zeros(N+1,N+1) little_D])

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
      safe(-q + (γ-1)*e.*little_D*v)
    ]
  end

  F(Qk, Qn) = (Qk - Qn)/dt + big_D*convection(Qk) + source(Qk)

# Jacobian

  J1 = 1/dt*sparse(I,N+1,N+1)

  J2 = spdiagm(
      -1 => -1*ones(N),
       1 =>    ones(N)
     )
  J2[1,:] = zeros(N+1)
  J2[end,end] = 2
  J2[end,end-1] = -2
  J2 = 0.5/dz*J2

  J3 = zeros(N+1,N+1)

  J6 = 0.5*(γ-1)/dz*spdiagm(
      -1 => -1*ones(N),
       1 =>    ones(N)
     )
  J6[1,:] = zeros(N+1)
  J6[end,end-1] = -(γ-1)/dz
  J6[end,end]   =  (γ-1)/dz

  function Jac(Q)
      ρ = Q[1:N+1]
      G = Q[N+2:2N+2]
      e = Q[2N+3:3N+3]

      # J4
      lower_diagonal_4 = [  0.5/dz *(G[j]/ρ[j])^2     for j in 1:N  ]
      diagonal_4       = [ -0.5*f/D*(G[j]/ρ[j])^2 + g for j in 1:N+1]
      upper_diagonal_4 = [ -0.5/dz *(G[j]/ρ[j])^2     for j in 2:N+1]
      J4 = spdiagm(              
      -1 => lower_diagonal_4,
       0 => diagonal_4,
       1 => upper_diagonal_4)
      J4[1,:]     = zeros(N+1)
      J4[N+1,N]   = 1/dz*(G[N]/ρ[N])^2
      J4[N+1,N+1] = -(G[N+1]/ρ[N+1])^2*(1/dz + 0.5*f/D) + g

      # J5
      lower_diagonal_5 = [     -1/dz*G[j]/ρ[j] for j in 1:N  ]
      diagonal_5       = [1/dt + f/D*G[j]/ρ[j] for j in 1:N+1]
      upper_diagonal_5 = [      1/dz*G[j]/ρ[j] for j in 2:N+1]
      J5 = spdiagm(              
      -1 => lower_diagonal_5,
       0 => diagonal_5,
       1 => upper_diagonal_5)
     J5[1,:]     = [1/dt; zeros(N)]
     J5[N+1,N]   = -2/dz*G[N]*ρ[N]
     J5[N+1,N+1] = 1/dt + G[N+1]/ρ[N+1]*(2/dz + f/D)

     # J7
     lower_diagonal_7 = [ 0.5/dz*G[j-1]/ρ[j-1]^2*(e[j-1] *(γ-1)*e[j]) for j in 2:N+1]
     upper_diagonal_7 = [-0.5/dz*G[j+1]/ρ[j+1]^2*(e[j+1] *(γ-1)*e[j]) for j in 1:N  ]
      J7 = spdiagm(              
      -1 => lower_diagonal_7,
       1 => upper_diagonal_7)
      J7[1,:] = zeros(N+1)
      J7[N+1,N  ] =  1/dz*G[N]/ρ[N]^2*(e[N] + (γ-1)*e[N+1])
      J7[N+1,N+1] = -γ/dz*e[N+1]*G[N+1]/ρ[N+1]^2

     # J8
     lower_diagonal_8 = [-(0.5/dz)/ρ[j-1]*(e[j-1] + (γ-1)*e[j]) for j in 2:N+1]
     upper_diagonal_8 = [ (0.5/dz)/ρ[j+1]*(e[j+1] + (γ-1)*e[j]) for j in 1:N  ]
      J8 = spdiagm(              
      -1 => lower_diagonal_8,
       1 => upper_diagonal_8)
     J8[1,:] = zeros(N+1)
     J8[N+1,N  ] = -1/dz/ρ[N]*(e[N] - e[N+1])
     J8[N+1,N+1] = γ*e[N+1]/ρ[N+1]

     # J9
     lower_diagonal_9 = [-0.5/dz*G[j]/ρ[j] for j in 1:N]
     diagonal_9       = [0; [1/dt + 0.5(γ-1)/dz*(G[j+1]/ρ[j+1] - G[j-1]/ρ[j-1]) for j in 2:N]; 0]
     upper_diagonal_9 = [ 0.5/dz*G[j]/ρ[j] for j in 2:N+1]
      J9 = spdiagm(              
      -1 => lower_diagonal_9,
       0 => diagonal_9,
       1 => upper_diagonal_9)
     J9[1,:] = [1/dt; zeros(N)]
     J9[N+1,N  ] = -1/dz*G[N]/ρ[N]
     J9[N+1,N+1] = 1/dt + 1/dz*G[N+1]/ρ[N+1] + (γ-1)/dz*(G[N+1]/ρ[N+1] - G[N]/ρ[N])

      return [J1 J2 J3; J4 J5 J6; J7 J8 J9]
      # return J6
  end

# Initialize fields
  ρ0 = 0.5
  v0 = 1.0
  T0 = 110 + 273.15
  # T0 = 1e
  U0 = Cv*T0
  
  # q = [zi < Lh ? Qval/(Lh*A_flow) : 0 for zi in z]
  q = [zi < Lh ? Qval : 0 for zi in z]

  ρ_initial = ρ0*   ones(N+1)
  G_initial = ρ0*v0*ones(N+1)
  e_initial = ρ0*U0*ones(N+1)

  # ρ_initial = ones(N+1)
  # G_initial = 5*ones(N+1)
  # e_initial = 10*ones(N+1)

  Qn = [ρ_initial; G_initial; e_initial]

  ρ_save = [ρ_initial]
  G_save = [G_initial]
  e_save = [e_initial]
  t_save = [0.0]

# Main loop
  tol = 1e-6
  Qk = copy(Qn)
  for n in 1:nt-1
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
          println("$(round(n/nt*100, digits=1))% ")
      end
      Qn = copy(Qk)
  end

# Calculate primitive fields
  n_save = length(ρ_save)
  v_save = [G_save[j]./ρ_save[j] for j in 1:n_save]
  U_save = [e_save[j]./ρ_save[j] for j in 1:n_save]
  p_save = [     (γ-1)*e_save[j] for j in 1:n_save]
  T_save = [U_save[j]/Cv .- 273.15 for j in 1:n_save]

  field_dict = Dict(
    "ρ" => ρ_save,
    "G" => G_save,
    "e" => e_save,

    "p" => p_save,
    "v" => v_save,
    "T" => T_save
  )

# Plotting

  select = ["G", "v", "T"]
  
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
    plots = [plot(z,field_dict[name][n], ylims=(field_min[name], field_max[name]), title=name, formatter=:plain) for name in select]
    xlabel!("t = $(t_save[n])")
    p = plot(plots..., layout=(length(select), 1))
    display(p)
  end