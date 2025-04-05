using Plots
using CoolProp
using SparseArrays

# Geometry
  D = 9.1e-3
  H = 2.3275
  W = 2.02
  L = 2(H + W)
  Lh = 0.575
  Lc = 0.2
  s0_heater = W + 0.05
  s0_cooler = 2W + H + 0.05

# Grid
  N = 100
  ds = L/N
  s = 0:ds:L
  s = copy(s[1:end-1])

# Constants
  Qh = 3e3
  A  = 0.25*π*D^2
  e = 1e-6

  g = zeros(N)
  g[W .< s .< W + H] .= -9.8
  g[2W + H .< s .< 2W + 2H] .= 9.8

# Matrices
  upwind_D = spdiagm(
      -1 => -ones(N-1),
       0 =>  ones(N)
  )
  upwind_D = 1/ds * upwind_D
  
  average_M = spdiagm(
      -1 => ones(N-1),
       0 => ones(N)
  )
  average_M = 0.5*average_M

  updiag(vector) = 1/ds*spdiagm(
    -1 => -vector[1:end-1],
     0 =>  vector
  )

# Functions 

  function calculate_rho(P,h)
      ρ = [PropsSI("D", "P", P[i], "H", h[i], "Water") for i in 1:N]
      return ρ
   end
  
  function calculate_v(ρ)
    global v0
      A = updiag(ρ)
      A[1,1] = 1
      b = [v0; zeros(N-1)]
      v = A\b
      return v
  end
  
  function calculate_P(ρ,v,f)
      A = 1 ./ (average_M*ρ) .* upwind_D
      b = g - 0.5*f/D .* v.^2 - (average_M*v) .* (upwind_D*v)
      A[1,:] = [1; zeros(N-1)]
      b[1] = P0
      P = A\b
      return P
  end
  
  function calculate_h(ρ,v,P,Qc)
    global h0, S
      A = updiag(ρ.*v)
      b = (average_M*v) .* (upwind_D*P) + S(Qc)
      A[1,:] = [1; zeros(N-1)]
      b[1] = h0
      h = A\b
      return h
  end

  function calculate_f(Re)
    global e
    f1 = (-2.457*log.((7 ./ Re).^0.9 .+ 0.27*e/D)).^16
    f2 = (37530 ./ Re).^12
    f  = 8*((8 ./ Re).^12 + (f1 + f2).^(-1.5)).^(1/12)
    return f
  end

  function calculate_Re(ρ,v,P,h)
    μ = [PropsSI("V", "P", P[i], "H", h[i], "Water") for i in 1:N]
    Re = ρ.*v*D./μ
    return Re
  end

  function S(Qc)
    S_vec = zeros(N)
    S_vec[s0_heater .< s .<  s0_heater + Lh] .= Qh/(A*Lh)
    S_vec[s0_cooler .< s .< s0_cooler + Lc] .= -Qc/(A*Lc)
    return S_vec
  end

# Initialize

  Qc = Qh
  v0 = 0.3
  P0 = 57.98e5

  T_sat = PropsSI("T", "P", P0, "Q", 0, "Water")
  T_subcooling = 12
  T0 = T_sat - T_subcooling
  h0 = PropsSI("H", "P", P0, "T", T0, "Water")

  ρ = 1000*ones(N)
  v = calculate_v(ρ)
  f = 0.02*ones(N)
  P = calculate_P(ρ,v,f)
  h = calculate_h(ρ,v,P,Qc)

# Main loop
 
  crit = 500
  nit = 0
  nitmax = 500
  while abs(P[1] - P[end]) > crit && nit < nitmax
    global ρ, v, P, h, v0, h0, nit, Qc
    ρ = calculate_rho(P,h)
    Re = calculate_Re(ρ,v,P,h)
    f = calculate_f(Re)
    v0 = v0*(P[end]/P[1])^4
    Qc = Qc*(h[end]/h[1])^8
    h0 = h[end]
    v = calculate_v(ρ)
    P = calculate_P(ρ,v,f)
    h = calculate_h(ρ,v,P,Qc)
    nit = nit + 1
  end
  println("Iterations: $nit")
  println("P[1] - P[end] = $(round(abs(P[1] - P[end]),digits=2)) [Pa]")
  println("Qc = $Qc [W]")
  println("m = $(ρ[1]*v[1]*A) [kg/s]")

# Derived fields
 
  function calculate_T(P,h)
    T = [PropsSI("T", "P", P[i], "H", h[i], "Water") for i in 1:N]
    T = T .- 273.15
    return T
  end
  
  function calculate_Q(P,h)
    Q = zeros(N)
    for i in 1:N
      Q[i] = PropsSI("Q", "P", P[i], "H", h[i], "Water")
      if Q[i] == -1
        Q[i] = 0
      end
    end
    return Q
  end

  function calculate_alpha(P,Q)
    α = zeros(N)
    x = Q
    for i in 1:N
      ρg = PropsSI("D", "P",P[i],"Q",1,"Water")
      ρl = PropsSI("D", "P",P[i],"Q",0,"Water")
      α[i] = (x[i]/ρg)/(x[i]/ρg + (1-x[i])/ρl)
    end
    return α
  end
  
  T = calculate_T(P,h)
  Q = calculate_Q(P,h)
  α = calculate_alpha(P,Q)
  P = P/1e5
  h = h/1e3
  
  field_dict = Dict(
    "ρ" => ρ,
    "v" => v,
    "P" => P,
    "h" => h,
    "T" => T,
    "Q" => Q,
    "α" => α
  )
  
# Plot

select = ["ρ", "v", "P", "Q"]
select = ["v", "α", "P", "h"]

plots = []
for name in select
  field = field_dict[name]
  a_plot = plot(s, field, title=name)
  push!(plots, a_plot)
end
plot(plots...)