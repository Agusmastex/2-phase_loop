using Plots
using CoolProp
using SparseArrays
using BandedMatrices

@time begin
# Geometry
  D = 9.1e-3
  H = 2.3275
  W = 2.02
  L = 2(H + W)
  Lh = 0.575
  Lc = 0.2
  s0_heater = W + 0.05
  s0_cooler = 2W + H + 0.05

# Auxiliary functions 

  function plot_this(list, names)
    plots = []
    for (item, name) in zip(list, names)
      p = plot(s,item, title=name)
      push!(plots,p)
    end
    p = plot(plots...)
    return p
  end

# Grid
  N = 20
  ds = L/N
  s = 0:ds:(L-ds)

# Constants
  Qh = 3e3
  A  = 0.25*π*D^2
  e = 1e-6

  g = zeros(N)
  g[W .< s .< W + H] .= -9.8
  g[2W + H .< s .< 2W + 2H] .= 9.8

# This will later be variable

  function S(Qc)
    S_vec = zeros(N)
    S_vec[s0_heater .< s .<  s0_heater + Lh] .= Qh/(A*Lh)
    S_vec[s0_cooler .< s .< s0_cooler + Lc] .= -Qc/(A*Lc)
    return S_vec
  end

# Matrices
  upwind_D = spdiagm(
      -1 => -ones(N-1),
       0 =>  ones(N)
  )
  #   upwind_D[1,N] = -1
  upwind_D = 1/ds * upwind_D
  
  average_M = spdiagm(
      -1 => ones(N-1),
       0 => ones(N)
  )
  average_M[1,N] = 1
  average_M = 0.5*average_M

  updiag(vector) = 1/ds*spdiagm(
    -1 => -vector[1:end-1],
     0 =>  vector
  )

# Functions 

  function calculate_rho!(ρ,P,h)
      ρ .= [PropsSI("D", "P", P[i], "H", h[i], "Water") for i in range(1,N)]
  end
  
  Av = updiag(ones(N))
  bv = zeros(N)
  v = zeros(N)
  function calculate_v!(v,Av,ρ,bv,v0)
    Av[Band(0)] .= ρ
    Av[Band(-1)] .= -ρ[1:end-1]
      Av[1,1] = 1.
      bv[1] = v0
      v .= Av\bv
  end
  
  Ap = copy(upwind_D)
  cache = zeros(N)
  P = zeros(N)
  bp = zeros(N)
  function calculate_P!(P,Ap,ρ,v,f,bp,cache)
    cache .= 1 ./ (average_M*ρ)
      Ap[Band(0)] .=  cache .* upwind_D[Band(0)]
      Ap[Band(-1)] .=  cache[2:end] .* upwind_D[Band(-1)]
      bp .= g - 0.5*f/D .* v.^2 - (average_M*v) .* (upwind_D*v)
      Ap[1,1] = 1.
      bp[1] = P0
      P .= Ap\bp
  end
  
  function calculate_h(ρ,v,P,Qc,h0)
    global S
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
    μ = zeros(N)
    for i in 1:N
      μ[i] = PropsSI("V", "P", P[i], "H", h[i], "Water")
    end
    Re = ρ.*v*D./μ
    return Re
  end

# Initialize

  Qc = Qh
  v0 = 0.5
  P0 = 57.98e5

  T_sat = PropsSI("T", "P", P0, "Q", 0, "Water")
  T_subcooling = 12
  T0 = T_sat - T_subcooling
  h0 = PropsSI("H", "P", P0, "T", T0, "Water")

  ρ = 1000*ones(N)
  calculate_v!(v,Av,ρ,bv,v0)
  f = 0.02*ones(N)
  calculate_P!(P,Ap,ρ,v,f,bp,cache)
  h = calculate_h(ρ,v,P,Qc,h0)
  # display(plot_this([ρ,v,P,h], ["ρ","v","P","h"]))


# Main loop
 
  crit = 500
  nit = 0
  nitmax = 500
  function advance(nit, v0, h0, Qc)
    while abs(P[1] - P[end]) > crit && nit < nitmax
      global h
      calculate_rho!(ρ,P,h)
      Re = calculate_Re(ρ,v,P,h)
      f = calculate_f(Re)
      v0 = v0*(P[end]/P[1])^4
      Qc = Qc*(h[end]/h[1])^8
      h0 = h[end]
      calculate_v!(v,Av,ρ,bv,v0)
      calculate_P!(P,Ap,ρ,v,f,bp,cache)
      h = calculate_h(ρ,v,P,Qc,h0)
      #display(plot_this([ρ,v,P,h], ["ρ","v","P","h"]))
      # println("Qc = $Qc")
      nit = nit + 1
    end
    return nit, v0, h0, Qc
  end
  total_it, v0_final, h0_final, Qc_final = advance(nit, v0, h0, Qc)
  println("Iterations: $total_it")
  println("P[1] - P[end] = $(round(abs(P[1] - P[end]),digits=2)) [Pa]")
  println("Qc = $Qc_final [W]")
  println("m = $(ρ[1]*v[1]*A) [kg/s]")

end
  
  
# Derived fields
 
  function calculate_T(P,h)
    T = zeros(N)
    for i in 1:N
      T[i] = PropsSI("T", "P", P[i], "H", h[i], "Water")
    end
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
      ρg = PropsSI("D","P",P[i],"Q",1,"Water")
      ρl = PropsSI("D","P",P[i],"Q",0,"Water")
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
  
  println()

# Plot

# s = circshift(s, 1)

# select = ["ρ", "v", "P", "Q"]
select = ["v", "α", "P", "h"]
select = ["ρ", "v", "h", "P"]

plots = []
for name in select
  field = field_dict[name]
  a_plot = plot(s, field, title=name)
  # a_plot = plot(s, circshift(field,1), title=name)
  # xaxis = minimum(field)
  # plot!([s0_heater, s0_heater + Lh], [xaxis, xaxis], marker=:circle)
  # plot!([s0_cooler, s0_cooler + Lc], [xaxis, xaxis], marker=:circle)
  push!(plots, a_plot)
end
plot(plots...)