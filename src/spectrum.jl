using multimode
using QuantumOptics

# Particles
const E0 = 0.2
const Nparticles = 1
const Nlevels = 14
const potential = multimode.BoxPotential
const particles = multimode.Bosons(Nparticles, Nlevels, E0, potential)

# Cavity mode
const index1 = 12
const cavitymodes1 = 14
const delta1 = -4.
# const eta1 = 2.0
const kappa1 = 1.
const U0_1 = -2.


# Multimode system
const scaling = 1./4


function spectrum(system, T)
    Ψ₀ = basis_ket(multimode.quantum.basis(system), 1)
    ρ₀ = Ψ₀⊗dagger(Ψ₀)
    Tss = [0,100]
    tout, ρ_t = multimode.quantum.timeevolution_master(system, Tss, ρ₀)
    ρ_ss = ρ_t[end]

    a = embed(multimode.quantum.basis(system), [2], [destroy(multimode.quantum.basis(system.modes[1]))])
    aρ_ss = a*ρ_ss

    x = Complex128[]
    fout(t, rho) = push!(x, expect(dagger(a), rho))

    multimode.quantum.timeevolution_master(system, T, aρ_ss; fout=fout)
    N = length(x)
    x = [reverse(x), x]
    S = abs2(fft(x))
    S = [S[N+1:end], S[1:N]]
    S /= maximum(S)
    return S
end

T = [0:0.1:250]
#etas = [0.5:0.1:3.5]
etas = [0.1:0.5:3.5]
#etas = [2.0:0.1:2.0]
S = zeros(Float64, length(etas), 2*length(T))
for (i, eta1) in enumerate(etas)
    println("η=", eta1)
    mode1 = multimode.CavityMode(index1, cavitymodes1, delta1, eta1, kappa1, U0_1)
    system = multimode.MultimodeSystem(scaling, particles, [mode1])
    S[i,:] = spectrum(system, T)
end

writecsv("spectrum_mode12_b.csv", S)

# using PyCall
# @pyimport matplotlib.pyplot as plt
# plt.plot(S)
# plt.show()
