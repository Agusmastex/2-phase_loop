using Plots
using SparseArrays

tf = 125.0
dt = 0.1
dz = 0.05
  
# Change this
  f = 0.5
  Qval = 1e4
  Δn = 10
  β  = 400e-6
  Vz_init = 0.0

# Grid
  L  = 0.5
  z  = 0:dz:L
  t  = 0:dt:tf
  nt = length(t) - 1
  N  = length(z) - 1

# Constants
  D  = 0.05
  Lh = 0.5*L
  ρ0 = 1000
  Cp = 4184
  T0 = 25
  A_flow = 0.25*π*D^2
  Q_vec = [zi <= Lh ? Qval/(ρ0*Cp*Lh*A_flow) : 0 for zi in z]
  Q_vec[1] = 0 # Enforce boundary condition

# Functions
  upwind_D = spdiagm(
    -1 => -ones(N),
     0 =>  [0; ones(N)])
  upwind_D = 1/dz*dropzeros(upwind_D)

  back3_D = spdiagm(
    -2 => ones(N-1),
    -1 => -4*ones(N),
     0 => 3*ones(N+1)
  )
  back3_D[1,:] = zeros(N+1)
  back3_D[2,:] = [-2; 2; zeros(N-1)]
  back3_D = 0.5/dz*dropzeros(back3_D)

  slope(Vz, T_avg) = β*(T_avg - T0) - 0.5*f/D*Vz*abs(Vz)
  average(x) = (sum(x[2:end-1]) + 0.5(x[1] + x[end]))/(length(x) - 1)


# Initialize fields
  T_init = 25.0*ones(N+1)

  Vz_save = [Vz_init]
  T_save  = [T_init]
  t_save  = [0.0]

# Solve

  # RK2 + Beam-Warming
  T_init = 25.0*ones(N+1)

  Vz_save = [Vz_init]
  T_save  = [T_init]
  t_save  = [0.0]

  T_old = T_init
  Vz_old = Vz_init
  
  for n in 1:nt
    global T_old, Vz_old
    Vz_half = Vz_old + 0.5*dt*slope(Vz_old,average(T_old))
    Vz = Vz_old + dt*slope(Vz_half, average(T_old))
    T  = T_old  + dt*(Q_vec - Vz_old*back3_D*T_old)

    if n % Δn == 0
      push!(Vz_save, Vz)
      push!(T_save, T)
      push!(t_save, t[n+1])
    end

    T_old = copy(T)
    Vz_old = copy(Vz)
  end
  Vz_save_BW_RK2 = copy(Vz_save)
  T_save_BW_RK2 = copy(T_save)

  # RK2 + Upwind
  T_init = 25.0*ones(N+1)

  Vz_save = [Vz_init]
  T_save  = [T_init]
  t_save  = [0.0]

  T_old = T_init
  Vz_old = Vz_init

  for n in 1:nt
    global T_old, Vz_old
    Vz_half = Vz_old + 0.5*dt*slope(Vz_old,average(T_old))
    Vz = Vz_old + dt*slope(Vz_half, average(T_old))
    T  = T_old  + dt*(Q_vec - Vz_old*upwind_D*T_old)

    if n % Δn == 0
      push!(Vz_save, Vz)
      push!(T_save, T)
      push!(t_save, t[n+1])
    end

    T_old = copy(T)
    Vz_old = copy(Vz)
  end
  Vz_save_upwind_RK2 = copy(Vz_save)
  T_save_upwind_RK2 = copy(T_save)

  # Forward Euler + Upwind
  T_init = 25.0*ones(N+1)

  Vz_save = [Vz_init]
  T_save  = [T_init]
  t_save  = [0.0]

  T_old = T_init
  Vz_old = Vz_init

  for n in 1:nt
    global T_old, Vz_old
    Vz = Vz_old + dt*slope(Vz_old,average(T_old))
    T  = T_old  + dt*(Q_vec - Vz_old*upwind_D*T_old)

    if n % Δn == 0
      push!(Vz_save, Vz)
      push!(T_save, T)
      push!(t_save, t[n+1])
    end

    T_old = copy(T)
    Vz_old = copy(Vz)
  end
  Vz_save_upwind_euler = copy(Vz_save)
  T_save_upwind_euler = copy(T_save)
  
# Steady-state

Vz_steady = 2*(β*Qval/(ρ0*Cp*π*D*f)*(1 - 0.5*Lh/L))^(1/3)
T_steady_max = T0 + Qval/(ρ0*Cp*A_flow*Vz_steady)

# Plotting
  n_save = length(Vz_save)
  for j in 1:n_save
    p1 = plot(z,[
      T_save_upwind_euler[j],
      T_save_upwind_RK2[j],
      # T_save_BW_RK2[j]
    ])
    hline!([T_steady_max])
    p2 = plot(t_save[1:j],[
      Vz_save_upwind_euler[1:j],
      Vz_save_upwind_RK2[1:j],
      # Vz_save_BW_RK2[1:j]
    ])
    hline!([Vz_steady])
    xlabel!("t = $(t_save[j])")
    p = plot(p1,p2, layout=(2,1))
    display(p)
  end

