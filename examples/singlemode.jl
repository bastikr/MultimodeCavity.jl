using quantumoptics
using multimode
using Optim

# Particles
const E0 = 0.2
const Nparticles = 1
const Nlevels = 15
const potential = multimode.quantum.BoxPotential
const particles = multimode.quantum.Bosons(Nparticles, Nlevels, E0, potential)

# Cavity mode
const index1 = 2
const cavitymodes1 = 20
const delta1 = -1.
const eta1 = 4.0
const kappa1 = 1.
const U0_1 = 0
const mode1 = multimode.quantum.CavityMode(index1, cavitymodes1, delta1, eta1, kappa1, U0_1)

# Multimode system
const scaling = 1.
const system = multimode.quantum.MultimodeSystem(scaling, particles, [mode1])

# Bases
basis_total = multimode.quantum.basis(system)
basis_particles = multimode.quantum.basis(particles)
basis_mode1 = multimode.quantum.basis(mode1)

# Initial state
const ψ₀ = basis_ket(basis_particles, 1) ⊗ basis_ket(basis_mode1, 1)
const ρ₀ = ψ₀ ⊗ dagger(ψ₀)

T = [0.,10.0]

tout, ρt = multimode.quantum.timeevolution_master(system, T, ρ₀, reltol=1e-7)
ρp = ptrace(ρt[end], [2])
ρf = ptrace(ρt[end], [1])
ρf /= trace(ρf)

using PyCall
@pyimport matplotlib.pyplot as plt

H = multimode.quantum.Hamiltonian(system)
U, S, V = svd(ρt[end].data)
x = Float64[-1:0.01:1]

for i=1:9
    ϕ = Ket(basis_total, vec(V[:,i]))
    q = qfunc(ptrace(ϕ⊗dagger(ϕ), [1]), [-4:0.1:4], [-4:0.1:4])
    plt.figure(1)
    plt.subplot(3,3,i)
    plt.imshow(q)

    n = multimode.quantum.particledensity(particles, x, ptrace(ϕ⊗dagger(ϕ), [2]))
    plt.figure(2)
    plt.subplot(3,3,i)
    plt.plot(x, n)
end

H = multimode.quantum.Hamiltonian(particles)
U, S, V = svd(ρp.data)

for i=1:9
    ϕ = Ket(basis_particles, vec(V[:,i]))
    n = multimode.quantum.particledensity(particles, x, ϕ⊗dagger(ϕ))
    plt.figure(3)
    plt.subplot(3,3,i)
    plt.plot(x, n)
end
plt.show()
