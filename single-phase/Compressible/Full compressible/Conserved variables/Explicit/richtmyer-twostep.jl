using Plots

tf = 0.005
dt = 0.00001
t  = 0:dt:tf
nt = length(t) - 1

L  = 1.0
dz = 0.01
z  = 0:dz:L
N  = length(z) - 1

γ = 1.4
D = 0.01
M = 18/1000
# g = 9.8
g = 0

function f(Q)
    ρ, m, e = Q
    return [m,
            m^2/ρ + e/M*(γ-1),
            e*m/ρ]
end


function ψ(Ql,Qc,Qr,i)
    ρi, mi, ei = Qc
    ρip1 = Qr[1]
    ρim1 = Ql[1]
    mip1 = Qr[2]
    mim1 = Ql[2]
    if i != N+1
        ∂V∂z = (mip1/ρip1 - mim1/ρim1)/(2*dz)
    else
        ∂V∂z = (mi/ρi - mim1/ρim1)/dz
    end
    return [0,
            -0.5*f₀/D*mi*abs(mi)/ρi - ρi*g,
            q[i] - 0*ei/M*(γ-1)*∂V∂z]
end


A = π*D^2/4
Lh = 0.5*L

Qval = 1
f₀ = 0.0

q = [zi < Lh ? Qval/(A*Lh) : 0 for zi in z]


T0 = 110 + 273.15
Cv = 1.5

ρ0 = 0.5
v0 = 1
U0 = Cv*T0

ρ = ρ0*ones(N+1)
m = ρ0*v0*ones(N+1)
e = ρ0*U0*ones(N+1)

Q = [ρ m e]

ρ_save = [ρ]
m_save = [m]
e_save = [e]

Q_new = copy(Q)

for n in 1:nt
    global Q, Qnew
    for i in 2:N
        Q_halfleft  = 0.5*(Q[i,:] + Q[i-1,:]) - 0.5*dt/dz*(f(Q[i,  :]) - f(Q[i-1,:])) #+ 0.5*dt*ψ(Q[i-1,:],Q[i,:],Q[i+1,:],i)
        Q_halfright = 0.5*(Q[i,:] + Q[i+1,:]) - 0.5*dt/dz*(f(Q[i+1,:]) - f(Q[i  ,:])) #+ 0.5*dt*ψ(Q[i-1,:],Q[i,:],Q[i+1,:],i)
        F_left  = f(Q_halfleft)
        F_right = f(Q_halfright)
        Q_new[i,:] = Q[i,:] - dt/dz*(F_right - F_left) + dt*ψ(Q[i-1,:], Q[i,:], Q[i+1,:], i)
    end
    Q[1,1] = ρ0
    Q[1,2] = ρ0*v0
    Q[1,3] = ρ0*U0
    Q = Q_new
    push!(ρ_save, Q[:,1])
    push!(m_save, Q[:,2])
    push!(e_save, Q[:,3])
end

for n in 1:nt
    p1 = plot(z, ρ_save[n])
    p2 = plot(z, m_save[n])
    p3 = plot(z, e_save[n],
    xlabel="t = $(t[n])")
    p = plot(p1,p2,p3, layout=(3,1))
    display(p)
end