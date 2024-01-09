using Plots

tf = 0.00000001
dt = 0.000000000001 
dz = 0.0001

# Change this
  Qval = 10000
  f₀ = 0
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
  A = 0.25*π*D^2
  Lh = 0.5*L
  Cv = 1.3e3

# Flux and source functions
  function f(Q)
      ρ, m, e = Q
      return [m,
              m^2/ρ + e*(γ-1),
              e*m/ρ]
  end
  
  function ψ(Ql,Qc,Qr,i)
      ρi, mi, ei = Qc
      ρip1 = Qr[1]
      ρim1 = Ql[1]
      mip1 = Qr[2]
      mim1 = Ql[2]
      ∂V∂z = (mip1/ρip1 - mim1/ρim1)/(2*dz)
      return [0,
              -0.5*f₀/D*mi*abs(mi)/ρi - ρi*g,
              q[i] - ei*(γ-1)*∂V∂z]
  end
  
  function ψend(Qleft, Qend)
      ρ,m,e = Qend
      ρl, ml, el = Qleft
      ∂V∂z = (m/ρ - ml/ρl)/dz
      return [0,
              -0.5*f₀/D*m*abs(m)/ρ - ρ*g,
               -e*(γ-1)*∂V∂z]
  end

# Initialize fields
  ρ0 = 0.5
  v0 = 1
  T0 = 110 + 273.15
  U0 = Cv*T0
  
  q = [zi < Lh ? Qval/(A*Lh) : 0 for zi in z]
  # q = [zi < Lh ? Qval : 0 for zi in z]
  
  ρ = ρ0   *ones(N+1)
  m = ρ0*v0*ones(N+1)
  e = ρ0*U0*ones(N+1)
  
  Q = [ρ m e]
  
  ρ_save = [ρ]
  m_save = [m]
  e_save = [e]
  t_save = [0.0]
  

# Main loop
  Q_new = copy(Q)
  @time for n in 1:nt
      global Q, Qnew
      for i in 2:N
          Q_new[i,:]  = 0.5*(Q[i-1,:] + Q[i+1,:]) - 0.5*dt/dz*(f(Q[i+1,:]) - f(Q[i-1,:])) + dt*ψ(Q[i-1,:], Q[i,:], Q[i+1,:], i)
      end
          Q_new[end,:] = Q[end,:] - dt/dz*(f(Q[end,:]) - f(Q[end-1,:])) + dt*ψend(Q[end-1,:], Q[end,:])
  
      Q[1,1] = ρ0
      Q[1,2] = ρ0*v0
      Q[1,3] = ρ0*U0
      Q = Q_new
      if n % Δn == 0 
          println("$(round(n/nt*100, digits=1)) %")
          push!(ρ_save, Q[:,1])
          push!(m_save, Q[:,2])
          push!(e_save, Q[:,3])
          push!(t_save, t[n])
      end
  end

# Calculate primitive fields
  n_save = length(ρ_save)
  v_save = [m_save[n]./ρ_save[n]   for n in 1:n_save]
  p_save = [e_save[n]*(γ-1)/1e5    for n in 1:n_save]
  U_save = [e_save[n]./ρ_save[n]   for n in 1:n_save]
  T_save = [U_save[n]/Cv .- 273.15 for n in 1:n_save]

  field_save = Dict(
    "ρ" => ρ_save,
    "G" => m_save,
    "e" => e_save,

    "p" => p_save,
    "v" => v_save,
    "T" => T_save
  )

# Plotting

select = ["p", "v", "T"]

for n in 1:n_save
  plots = [plot(z,field_save[name][n], title=name, formatter=:plain) for name in select]
  xlabel!("t = $(t_save[n])")
  p = plot(plots..., layout=(length(select), 1))
  display(p)
end