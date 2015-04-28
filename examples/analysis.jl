using quantumoptics
using multimode
using Optim

# Particles
const E0 = 0.2
const Nparticles = 1
const Nlevels = 12
const potential = multimode.BoxPotential
const particles = multimode.Bosons(Nparticles, Nlevels, E0, potential)

# Cavity mode
const index1 = 3
const cavitymodes1 = 12
const delta1 = -2.
const eta1 = 4.0
const kappa1 = 1.
const U0_1 = 0
const mode1 = multimode.CavityMode(index1, cavitymodes1, delta1, eta1, kappa1, U0_1)

# Multimode system
const scaling = 1./3
const system = multimode.MultimodeSystem(scaling, particles, [mode1])

# Bases
basis_particles = multimode.basis(particles)
basis_mode1 = multimode.basis(mode1)

# Initial state
const ψ₀ = basis_ket(basis_particles, 1) ⊗ basis_ket(basis_mode1, 1)
const ρ₀ = ψ₀ ⊗ dagger(ψ₀)

function fmin(x)
    N = length(x)//2
    y = x[1:N] + 1im*x[N+1:end]
    y /= norm(y,2)
    return multimode.multimode_steadystate.Hfunctional(system, Ket(basis_particles, y))
end
y0 = basis_ket(basis_particles, 1).data
N = length(y0)
x0 = [real(y0),imag(y0)]
x = Optim.optimize(fmin, x0).minimum
y = x[1:N] + 1im*x[N+1:end]
y /= norm(y,2)
x1_opt = Ket(basis_particles, y)
ρg1_opt = x1_opt ⊗ dagger(x1_opt)


function fmin2(x)
    N = length(x)//2
    y = x[1:N] + 1im*x[N+1:end]
    y /= norm(y,2)
    y = Ket(basis_particles, y)
    y_ = y - (dagger(x1_opt)*y)*x1_opt
    y_ /= sqrt(dagger(y_)*y_)
    return multimode.multimode_steadystate.Hfunctional(system, y_)
end
y0 = basis_ket(basis_particles, 1).data
N = length(y0)
x0 = [real(y0),imag(y0)]
x = Optim.optimize(fmin2, x0).minimum
y = x[1:N] + 1im*x[N+1:end]
y /= norm(y,2)
y = Ket(basis_particles, y)
y_ = y - (dagger(x1_opt)*y)*x1_opt
n = sqrt(dagger(y_)*y_)
y_ /= n
println("Norm: ", dagger(y_)*y_)
x2_opt = y_
ρg2_opt = x2_opt ⊗ dagger(x2_opt)


function fmin3(x)
    N = length(x)//2
    y = x[1:N] + 1im*x[N+1:end]
    y /= norm(y,2)
    y = Ket(basis_particles, y)
    y_ = y - (dagger(x1_opt)*y)*x1_opt - (dagger(x2_opt)*y)*x2_opt
    y_ /= sqrt(dagger(y_)*y_)
    return multimode.multimode_steadystate.Hfunctional(system, y_)
end
y0 = basis_ket(basis_particles, 1).data
N = length(y0)
x0 = [real(y0),imag(y0)]
x = Optim.optimize(fmin3, x0).minimum
y = x[1:N] + 1im*x[N+1:end]
y /= norm(y,2)
y = Ket(basis_particles, y)
y_ = y - (dagger(x1_opt)*y)*x1_opt - (dagger(x2_opt)*y)*x2_opt
n = sqrt(dagger(y_)*y_)
y_ /= n
println("Norm: ", dagger(y_)*y_)
x3_opt = y_
ρg3_opt = x3_opt ⊗ dagger(x3_opt)


T = [0.,100.0]

tout, ρt = multimode.timeevolution_master(system, T, ρ₀)

ρp = ptrace(ρt[end], [2])
ρf = ptrace(ρt[end], [1])
ρf /= trace(ρf)

U, S, V = svd(ρt[end].data)
println(S[1:5])


x1 = Ket(ρ₀.basis_l, U[:,1])
x2 = Ket(ρ₀.basis_l, U[:,2])
x3 = Ket(ρ₀.basis_l, U[:,3])
x4 = Ket(ρ₀.basis_l, U[:,4])
x5 = Ket(ρ₀.basis_l, U[:,5])

ρg1 = x1 ⊗ dagger(x1)
ρg2 = x2 ⊗ dagger(x2)
ρg3 = x3 ⊗ dagger(x3)
ρg4 = x4 ⊗ dagger(x4)
ρg5 = x5 ⊗ dagger(x5)

ρf1 = ptrace(ρg1, [1])
ρf2 = ptrace(ρg2, [1])
ρf3 = ptrace(ρg3, [1])
ρf4 = ptrace(ρg4, [1])
ρf5 = ptrace(ρg5, [1])

function f1(x)
    state = coherent_state(basis_mode1, complex(x[1], x[2]))
    return tracedistance(state ⊗ dagger(state), ρf1)
end
x = Optim.optimize(f1, Float64[0,0]).minimum
α1 = x[1] + 1im*x[2]

function f2(x)
    state = coherent_state(basis_mode1, complex(x[1], x[2]))
    return tracedistance(state ⊗ dagger(state), ρf2)
end
x = Optim.optimize(f2, Float64[0,0]).minimum
α2 = x[1] + 1im*x[2]

function f3(x)
    state = coherent_state(basis_mode1, complex(x[1], x[2]))
    return tracedistance(state ⊗ dagger(state), ρf3)
end
x = Optim.optimize(f3, Float64[0,0]).minimum
α3 = x[1] + 1im*x[2]

function f4(x)
    state = coherent_state(basis_mode1, complex(x[1], x[2]))
    return tracedistance(state ⊗ dagger(state), ρf4)
end
x = Optim.optimize(f4, Float64[0,0]).minimum
α4 = x[1] + 1im*x[2]

function f5(x)
    state = coherent_state(basis_mode1, complex(x[1], x[2]))
    return tracedistance(state ⊗ dagger(state), ρf5)
end
x = Optim.optimize(f5, Float64[0,0]).minimum
α5 = x[1] + 1im*x[2]



println("α₁ = ", α1)
println("α₂ = ", α2)
println("α₃ = ", α3)
println("α₄ = ", α4)
println("α₅ = ", α5)

Halpha1 = multimode.multimode_steadystate.Halpha(system, α1)
Halpha2 = multimode.multimode_steadystate.Halpha(system, α2)
Halpha3 = multimode.multimode_steadystate.Halpha(system, α3)

B_ = multimode.quantum.B(1, system)
D, V = eig(Hermitian(Halpha1.data))
x1_ = Ket(basis_particles, vec(V[:,1]))
x2_ = Ket(basis_particles, vec(V[:,2]))
x3_ = Ket(basis_particles, vec(V[:,3]))
x4_ = Ket(basis_particles, vec(V[:,4]))
x5_ = Ket(basis_particles, vec(V[:,5]))

println("α₁ = ", eta1*expect(B_, x1_)/(delta1+1im*kappa1))
println("α₂ = ", eta1*expect(B_, x2_)/(delta1+1im*kappa1))
println("α₃ = ", eta1*expect(B_, x3_)/(delta1+1im*kappa1))
println("α₄ = ", eta1*expect(B_, x4_)/(delta1+1im*kappa1))
println("α₅ = ", eta1*expect(B_, x5_)/(delta1+1im*kappa1))

println("α₁opt = ", eta1*expect(B_, x1_opt)/(delta1+1im*kappa1))
println("α₂opt = ", eta1*expect(B_, x2_opt)/(delta1+1im*kappa1))
println("α₃opt = ", eta1*expect(B_, x3_opt)/(delta1+1im*kappa1))

ρg1_ = x1_ ⊗ dagger(x1_)
ρg2_ = x2_ ⊗ dagger(x2_)
ρg3_ = x3_ ⊗ dagger(x3_)
ρg4_ = x4_ ⊗ dagger(x4_)
ρg5_ = x5_ ⊗ dagger(x5_)

x = linspace(-1,1,100)

n1 = multimode.particledensity(particles, x, ptrace(ρg1, [2]))
n2 = multimode.particledensity(particles, x, ptrace(ρg2, [2]))
n3 = multimode.particledensity(particles, x, ptrace(ρg3, [2]))
n4 = multimode.particledensity(particles, x, ptrace(ρg4, [2]))
n5 = multimode.particledensity(particles, x, ptrace(ρg5, [2]))

n1_ = multimode.particledensity(particles, x, ρg1_)
n2_ = multimode.particledensity(particles, x, ρg2_)
n3_ = multimode.particledensity(particles, x, ρg3_)
n4_ = multimode.particledensity(particles, x, ρg4_)
n5_ = multimode.particledensity(particles, x, ρg5_)

n1_opt = multimode.particledensity(particles, x, ρg1_opt)
n2_opt = multimode.particledensity(particles, x, ρg2_opt)
n3_opt = multimode.particledensity(particles, x, ρg3_opt)

using PyCall
@pyimport matplotlib.pyplot as plt
plt.figure()
plt.plot(x, n1, "k--")
plt.plot(x, n1_, "b")
plt.plot(x, n1_opt, "b--")
plt.plot(x, n2, "k--")
plt.plot(x, n2_opt, "g--")
plt.plot(x, n3_opt, "r--")
plt.plot(x, n2_, "g")
# plt.plot(x, n3, "k--")
plt.plot(x, n3_, "r")
# plt.plot(x, n4, "k--")
# plt.plot(x, n4_, "c")
# plt.plot(x, n5, "k--")
# plt.plot(x, n5_, "orange")
plt.show()
@assert false

tout, ρt = multimode.timeevolution_master(system, T, ρg1)
U, S, V = svd(ρt[end].data)

println(S[1:5])



ρg = S[1]*ρg1 + S[2]*ρg2

x = linspace(-1,1,100)
n1 = multimode.particledensity(particles, x, ptrace(ρg1, [2]))
n1 = multimode.particledensity(particles, x, ptrace(ρg1, [2]))
n2 = multimode.particledensity(particles, x, ptrace(ρg2, [2]))
n3 = multimode.particledensity(particles, x, ptrace(ρg3, [2]))
n4 = multimode.particledensity(particles, x, ptrace(ρg4, [2]))
n5 = multimode.particledensity(particles, x, ptrace(ρg5, [2]))

using PyCall
@pyimport matplotlib.pyplot as plt

q = qfunc(ptrace(ρg1, [1]), [-4:0.1:4], [-4:0.1:4])
plt.imshow(q)

plt.figure()
q = qfunc(ptrace(ρg2, [1]), [-4:0.1:4], [-4:0.1:4])
plt.imshow(q)

plt.figure()
q = qfunc(ptrace(ρg3, [1]), [-4:0.1:4], [-4:0.1:4])
plt.imshow(q)

plt.figure()
q = qfunc(ptrace(ρg4, [1]), [-4:0.1:4], [-4:0.1:4])
plt.imshow(q)

plt.figure()
q = qfunc(ptrace(ρg5, [1]), [-4:0.1:4], [-4:0.1:4])
plt.imshow(q)



# plt.figure()
# q = qfunc(ρf-rho_f_try, [-4:0.1:4], [-4:0.1:4])
# plt.imshow(q)
# plt.colorbar()

plt.figure()
plt.plot(x, n1)
plt.plot(x, n2)
plt.plot(x, n3)
plt.plot(x, n4)
plt.plot(x, n5)
plt.show()

# Hp = multimode.Hamiltonian(particle)
# B = multimode.B(1, system)

# ξ = 7.

# Hplus = Hp + ξ*B
# Hminus = Hp - ξ*B

# Dplus, Vplus = eig(Hplus.data)
# Dminus, Vminus = eig(Hminus.data)

# vplus1 = Ket(particle.basis, Vplus[:,1])
# vminus1 = Ket(particle.basis, Vminus[:,1])
# vplus2 = Ket(particle.basis, Vplus[:,2])
# vminus2 = Ket(particle.basis, Vminus[:,2])
# vplus3 = Ket(particle.basis, Vplus[:,3])
# vminus3 = Ket(particle.basis, Vminus[:,3])
# vplus4 = Ket(particle.basis, Vplus[:,4])
# vminus4 = Ket(particle.basis, Vminus[:,4])

# #vplus = sqrt(0.8)*vplus1 + sqrt(0.15)*vplus2 + sqrt(0.05)*vplus4
# #vminus = sqrt(0.8)*vminus1 + sqrt(0.15)*vminus2 + sqrt(0.05)*vminus4
# #vplus = sqrt(0.98)*vplus1 + sqrt(0.02)*vplus2
# #vminus = sqrt(0.98)*vminus1 + sqrt(0.02)*vminus2
# vplus = vplus1
# vminus = vminus1

# n_xplus = particledensity(particle, [-1:0.01:1], vplus)
# n_xminus = particledensity(particle, [-1:0.01:1], vminus)

# exp_B = dagger(vplus)*B*vplus
# η = sqrt((delta1^2 + kappa1^2)/(delta1*exp_B)*ξ)
# α = η*exp_B/(delta1 + 1im*kappa1)

# println("<B> = ", exp_B)
# println("η = ", η)
# println("α = ", α)
# println("n = ", abs2(α))

# # for i=4:5
# #     vplus = Ket(particle.basis, Vplus[:,i])
# #     vminus = Ket(particle.basis, Vminus[:,i])
# #     n_xplus = particledensity(particle, [-1:0.01:1], vplus)
# #     n_xminus = particledensity(particle, [-1:0.01:1], vminus)
# #     plt.plot((n_xplus + n_xminus)/2)
# # end

# # plt.show()
# # @assert false
# plt.plot(n_xplus)
# plt.plot(n_xminus)
# plt.figure()
# plt.plot((n_xplus + n_xminus)/2)
# #plt.plot(n_xplus)

# mode1.eta = η
# ρplus₀ = tensor(vplus, dagger(vplus))
# ρminus₀ = tensor(vminus, dagger(vminus))
# Ψfplus₀ = coherent_state(mode1.basis, α)
# Ψfminus₀ = coherent_state(mode1.basis, -α)
# #Ψf₀ = coherent_state(mode1.basis, 0)
# #println(expect(destroy(mode1.basis), Ψf₀))
# ρfplus₀ = tensor(Ψfplus₀, dagger(Ψfplus₀))
# ρfminus₀ = tensor(Ψfminus₀, dagger(Ψfminus₀))
# ρ₀ = 0.5*tensor(ρplus₀, ρfplus₀) + 0.5*tensor(ρminus₀, ρfminus₀)

# #println(expect(destroy(mode1.basis), ρf₀))
# #@assert false
# #ρ₀ = tensor(basis_ket(system.basis, 1), basis_bra(system.basis, 1))
# T, ρ_t = multimode.evolve_master(system, [0.:1:5.], ρ₀)
# exp_n = expect(number(mode1.basis), ptrace(ρ_t[end], [1]))
# println("<a^t a> = ", exp_n)
# n2 = particledensity(particle, [-1:0.01:1], ptrace(ρ_t[2], [2]))
# n3 = particledensity(particle, [-1:0.01:1], ptrace(ρ_t[3], [2]))
# n4 = particledensity(particle, [-1:0.01:1], ptrace(ρ_t[4], [2]))
# n = particledensity(particle, [-1:0.01:1], ptrace(ρ_t[end], [2]))
# q = qfunc(ptrace(ρ_t[end], [1]), [-3:0.1:3], [-3:0.1:3])

# using PyCall
# @pyimport matplotlib.pyplot as plt
# plt.plot(n2)
# plt.plot(n3)
# plt.plot(n4)
# plt.plot(n)
# plt.figure()
# plt.imshow(q)
# plt.show()
# #println([calc_exp_B(ρ) for ρ in ρ_t])
# #println(real(diag(ptrace(ρ_t[1], [2]).data)))
# #println(real(diag(ptrace(ρ_t[end], [2]).data)))


# #println(norm((Hplus*v1-D[1]*v1).data))

