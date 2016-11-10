using multimode
using QuantumOptics

E0 = 0.2
s = 0.25
Nparticles = 1
Nparticlemodes = 16

index1 = 12
cavitymodes1 = 14
delta1 = -4.
eta1 = 2.0
kappa1 = 1.
U0_1 = 0.

function QuantumOptics.qfunc(system)
    Ψ₀ = basis_ket(system.basis, 1)
    ρ₀ = Ψ₀⊗dagger(Ψ₀)
    Tss = [0,100]

    T, ρ_t = multimode.evolve_master(system, Tss, ρ₀)
    ρ_ss = ρ_t[end]

    rho_f = ptrace(ρ_ss, [1])
    rho_p = ptrace(ρ_ss, [2])

    x = linspace(-3,3,200)
    return qfunc(rho_f, x, x)
end


function n(system)
    Ψ₀ = basis_ket(system.basis, 1)
    ρ₀ = Ψ₀⊗dagger(Ψ₀)
    Tss = [0,20]

    T, ρ_t = multimode.evolve_master(system, Tss, ρ₀)
    ρ_ss = ρ_t[end]

    rho_f = ptrace(ρ_ss, [1])
    rho_p = ptrace(ρ_ss, [2])

    a = destroy(system.modes[1].basis)

    return expect(dagger(a)*a, rho_f)
end

#particle = Bosons(Nparticles, Nparticlemodes, E0, BoxPotential)
#eta1 = 3.0
#mode1 = CavityMode(index1, cavitymodes1, delta1, eta1, kappa1, U0_1)
#system = MultimodeSystem(s, particle, [mode1])
#Ψ₀ = basis_ket(system.basis, 1)
#println(n(system))
#println(multimode.evolve_mcwf(system, [0,1], Ψ₀))


# particle = Bosons(Nparticles, Nparticlemodes, E0, BoxPotential)
# eta1 = 3.0
# mode1 = CavityMode(index1, cavitymodes1, delta1, eta1, kappa1, U0_1)
# system = MultimodeSystem(s, particle, [mode1])
# M = qfunc(system)
# writecsv("qfunc_mode12_eta=3.0.csv", M)

# eta1 = 2.
#etas = [0:0.1:3.5]
deltas = [-3:0.1:3]
exp_n = Float64[]
particle = Bosons(Nparticles, Nparticlemodes, E0, BoxPotential)
for (i, delta1) in enumerate(deltas)
     println(delta1)
     mode1 = CavityMode(index1, cavitymodes1, delta1, eta1, kappa1, U0_1)
     system = MultimodeSystem(s, particle, [mode1])
     push!(exp_n, abs(n(system)))
end
#
# writecsv("expn_mode12.csv", exp_n)

using PyCall
@pyimport matplotlib.pyplot as plt

plt.plot(deltas, exp_n)
plt.show()


