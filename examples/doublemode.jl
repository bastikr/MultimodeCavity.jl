using quantumoptics
using multimode
using Optim

# Particles
const E0 = 0.2
const Nparticles = 1
const Nlevels = 6
const potential = multimode.BoxPotential
const particles = multimode.Bosons(Nparticles, Nlevels, E0, potential)

# Cavity mode
const index1 = 3
const cavitymodes1 = 10
const delta1 = -1.
const eta1 = 2.0
const kappa1 = 1.
const U0_1 = 0
const mode1 = multimode.CavityMode(index1, cavitymodes1, delta1, eta1, kappa1, U0_1)

# Cavity mode
const index2 = 6
const cavitymodes2 = 12
const delta2 = -1.
const eta2 = 3.0
const kappa2 = 1.
const U0_2 = 0
const mode2 = multimode.CavityMode(index2, cavitymodes2, delta2, eta2, kappa2, U0_2)

# Multimode system
const scaling = 1./3
const system = multimode.MultimodeSystem(scaling, particles, [mode1, mode2])

# Bases
basis_particles = multimode.basis(particles)
basis_mode1 = multimode.basis(mode1)
basis_mode2 = multimode.basis(mode2)

# Initial state
const ψ₀ = basis_ket(basis_particles, 1) ⊗ basis_ket(basis_mode1, 1) ⊗ basis_ket(basis_mode2, 1)
const ρ₀ = ψ₀ ⊗ dagger(ψ₀)

T = [0.:1.0]

tout, ρt = multimode.timeevolution_master(system, T, ρ₀, reltol=1e-7)
ρp = ptrace(ρt[end], [2,3])
ρf1 = ptrace(ρt[end], [1,3])
ρf2 = ptrace(ρt[end], [1,2])

function f1(x)
    state = coherent_state(basis_mode1, complex(x[1], x[2]))
    return tracedistance(state ⊗ dagger(state), ρf1)
end

function f2(x)
    state_p = coherent_state(basis_mode2, complex(x[1], x[2]))
    state_m = coherent_state(basis_mode2, complex(-x[1], -x[2]))
    return tracedistance(0.5*state_m ⊗ dagger(state_m) + 0.5*state_p ⊗ dagger(state_p), ρf2)
end

println(Optim.optimize(f1, Float64[0,0]))
println(Optim.optimize(f2, Float64[0,0]))

q1 = qfunc(ρf1, [-4:0.1:4], [-4:0.1:4])
q2 = qfunc(ρf2, [-4:0.1:4], [-4:0.1:4])

x = linspace(-1,1,100)
n = multimode.particledensity(particles, x, ρp)

using PyCall
@pyimport matplotlib.pyplot as plt
plt.imshow(q1)
plt.figure()
plt.imshow(q2)
plt.figure()
plt.plot(x, n)
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

