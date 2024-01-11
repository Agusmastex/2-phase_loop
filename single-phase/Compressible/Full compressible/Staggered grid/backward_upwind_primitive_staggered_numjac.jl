using Plots
using SparseArrays
using LinearAlgebra

tf = 1.0
dt = 0.01
dz = 0.1
  
# Change this
  f = 0
  Qval = 100
  g = 0
  Δn = 1

# Grid
  L  = 0.5
  z  = 0:dz:L
  t  = 0:dt:tf
  nt = length(t)
  N  = length(z) - 1
  
  z_cells = 0.5dz:dz:L
  z_edges = z[2:end]

# Constants
  γ = 1.4
  D = 0.05
  Lh = 0.5*L
  A_flow = 0.25*π*D^2
  Cv = 1.3e3

# Residual function
  upwind_D = spdiagm(
      -1 => -ones(N-1),
       0 =>  ones(N)
  )
  upwind_D[1,:] = zeros(N)
  upwind_D = 1/dz * dropzeros(upwind_D)
  
  forward_D = spdiagm(
      0 => -ones(N),
      1 =>  ones(N-1)
  )
  forward_D[1,:] = zeros(N)
  forward_D[end,:] = [zeros(N-2); -1; 1]
  forward_D = 1/dz * dropzeros(forward_D)

  average_D = spdiagm(
    0 => ones(N),
    1 => ones(N-1)
  )
  average_D[end,end] = 2
  average_D = 0.5 * dropzeros(average_D)

  function safe(vector)
    return [0; vector[2:end]]
  end

  function F(Qk, Qn)
    ρn = Qn[1:N]
    vn = Qn[N+1:2N]
    Un = Qn[2N+1:3N]

    ρ = Qk[1:N]
    v = Qk[N+1:2N]
    U = Qk[2N+1:3N]

    p = ρ.*U*(γ-1)

    return [
        (ρ - ρn)/dt + upwind_D*(ρ.*v);
        (v - vn)/dt + v.*upwind_D*(v) + (1 ./ (average_D*ρ)) .* forward_D*(p) + safe(0.5*f/D * (v.*v) .+ g);
        (ρ.*U - ρn.*Un)/dt + upwind_D*(ρ.*U.*v) - safe(q) + p.*upwind_D*(v)

    ]
  end


# Jacobian
  basis(i) = [ k == i ? 1 : 0 for k in 1:3N ]

  function Jac(Q, Qn)
    δ = 1e-4
    J = zeros(3N,3N)
    for j in 1:3N
      J[:,j] = (F(Q + δ*basis(j), Qn) - F(Q, Qn))/δ
    end
    return J
  end

# Initialize fields
  ρ0 = 0.5
  v0 = 1.0
  T0 = 110 + 273.15
  U0 = Cv*T0
  
  q = [zi <= Lh ? Qval/(Lh*A_flow) : 0 for zi in z_cells]

  ρ_initial = ρ0*ones(N)
  v_initial = v0*ones(N)
  U_initial = U0*ones(N)

  Qn = [ρ_initial; v_initial; U_initial]

  ρ_save = [ρ_initial]
  v_save = [v_initial]
  U_save = [U_initial]
  t_save = [0.0]

# Main loop
  tol = 1e-3
  Qk = copy(Qn)
  for n in 1:nt
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
          ρ = Qk[1:N]
          v = Qk[N+1:2N]
          U = Qk[2N+1:3N]
          push!(ρ_save, ρ)
          push!(v_save, v)
          push!(U_save, U)
          push!(t_save, t[n])
          println()
          println("$(round(n/nt*100, digits=1))% ")
      end
      Qn = copy(Qk)
  end

# Calculate derived fields
  n_save = length(ρ_save)
  p_save = [(γ-1)*ρ_save[j].*U_save[j]/1e3    for j in 1:n_save]
  T_save = [U_save[j]/Cv .- 273.15 for j in 1:n_save]

  field_dict = Dict(
    "ρ" => ρ_save,
    "v" => v_save,
    "U" => U_save,

    "p" => p_save,
    "T" => T_save
  )

# Plotting

for n in 1:n_save
   p1 = plot(z_cells, ρ_save[n], title = "ρ")
   p2 = plot(z_edges, v_save[n], title = "v")
   p3 = plot(z_cells, U_save[n], title = "U")
   p4 = plot(z_cells, p_save[n], title = "p")
   xlabel!("t = $(t_save[n])")
   p  = plot(p1,p2,p3,p4, xlims=(0,L), formatter=:plain)
  # display(p4)
   display(p)
end
  
  # simulation = "Backward Upwind Primitive Staggered NumJac"
  # p[:plot_title] = simulation
  # plot(p)
  # savefig(simulation)