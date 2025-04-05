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
  N = 50
  ds = L/N
  s = 0:ds:L
  s = copy(s[1:end-1])

# Constants
  A  = 0.25*π*D^2
  e = 1e-6
  P0 = 57.98e5
  T_sat = PropsSI("T", "P", P0, "Q", 0, "Water")
  T_subcooling = 11.12
  T0 = T_sat - T_subcooling

  g = zeros(N)
  for i in 1:N
    if W < s[i] < W + H
        g[i] = -9.8
    elseif 2W + H < s[i] < 2W + 2H
        g[i] = 9.8
    end
  end


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
      ρ = zeros(N)
      for i in 1:N
        ρ[i] = PropsSI("D", "P", P[i], "H", h[i], "Water")
      end
      return ρ
   end
  
  function calculate_v(ρ,v0)
      A = updiag(ρ)
      A[1,1] = 1
      b = [v0; zeros(N-1)]
      v = A\b
      return v
  end
  
  function calculate_P(ρ,v,f)
    global P0
      A = 1 ./ (average_M*ρ) .* upwind_D
      b = g - 0.5*f/D .* v.^2 - (average_M*v) .* (upwind_D*v)
      A[1,:] = [1; zeros(N-1)]
      b[1] = P0
      P = A\b
      return P
  end
  
  function calculate_h(ρ,v,P,h0,Qh,Qc)
      A = updiag(ρ.*v)
      b = (average_M*v) .* (upwind_D*P) + S(Qh,Qc)
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

  function S(Qh,Qc)
    S_vec = zeros(N)
    for i in 1:N
      if  s0_heater < s[i] <  s0_heater + Lh
          S_vec[i] =  Qh/(A*Lh)
      elseif s0_cooler < s[i] < s0_cooler + Lc
          S_vec[i] = -Qc/(A*Lc)
      end
    end
    return S_vec
  end

function solve(Qh)
  # Initialize
  
    Qc = Qh
    v0 = 0.5
    h0 = PropsSI("H", "P", P0, "T", T0, "Water")
  
    ρ = 1000*ones(N)
    v = calculate_v(ρ,v0)
    f = 0.02*ones(N)
    P = calculate_P(ρ,v,f)
    h = calculate_h(ρ,v,P,h0,Qh,Qc)
  
  # Main loop
   
    crit = 500
    nit = 0
    nitmax = 500
    while abs(P[1] - P[end]) > crit && nit < nitmax
      # global ρ, v, P, h
      ρ = calculate_rho(P,h)
      Re = calculate_Re(ρ,v,P,h)
      f = calculate_f(Re)
      v0 = v0*(P[end]/P[1])^4
      Qc = Qc*(h[end]/h[1])^8
      h0 = h[end]
      v = calculate_v(ρ,v0)
      P = calculate_P(ρ,v,f)
      h = calculate_h(ρ,v,P,h0,Qh,Qc)
      # display(plot_this([ρ,v,P,h], ["ρ","v","P","h"]))
      # println("Qc = $Qc")
      nit = nit + 1
    end

    m = ρ[1]*v[1]*A
    println("Iterations: $nit")
    println("P[1] - P[end] = $(round(abs(P[1] - P[end]),digits=2)) [Pa]")
    println("Qc = $Qc [W]")
    println("m = $m [kg/s]")

  diagnostic_dict = Dict(
    "nit" => nit,
    "Qc"  => Qc,
    "m"   => m,
  )

  # return [ρ,v,P,h], diagnostic_dict
  return Qc
  end