using multimode
using quantumoptics

E0 = 0.2
s = 0.25
Nparticles = 1
Nparticlemodes = 14

index1 = 12
cavitymodes1 = 14
delta1 = -4.
eta1 = 2.0
kappa1 = 1.
U0_1 = -2.


function spectrum(system, T)
    H = multimode.Hamiltonian(system)
    J = multimode.JumpOperators(system)
    Jdagger = map(dagger, J)
    #Hnh = H - 0.5im*Jdagger[1]*J[1] - 0.5im*Jdagger[2]*J[2]
    Hnh = H - 0.5im*Jdagger[1]*J[1]
    Hnh_dagger = dagger(Hnh)
    J_sparse = map(operators_sparse.SparseOperator, J)
    Jdagger_sparse = map(operators_sparse.SparseOperator, Jdagger)
    Hnh_sparse = operators_sparse.SparseOperator(Hnh)
    Hnh_dagger_sparse = operators_sparse.SparseOperator(Hnh_dagger)

    Ψ₀ = basis_ket(system.basis, 1)
    ρ₀ = Ψ₀⊗dagger(Ψ₀)
    Tss = [0,100]
    tout, ρ_t = timeevolution.master_nh(Tss, ρ₀, Hnh_sparse, J_sparse; Jdagger=Jdagger_sparse, Hdagger=Hnh_dagger_sparse)
    ρ_ss = ρ_t[end]

    a = embed(system.basis, [2], [destroy(system.modes[1].basis)])

    aρ_ss = a*ρ_ss

    x = Complex128[]

    fout(t, rho) = push!(x, expect(dagger(a), rho))
    timeevolution.master_nh(T, aρ_ss, Hnh_sparse, J_sparse; Jdagger=Jdagger_sparse, Hdagger=Hnh_dagger_sparse, fout=fout)
    N = length(x)
    x = [reverse(x), x]
    S = abs2(fft(x))
    S = [S[N+1:end], S[1:N]]
    S /= maximum(S)
    return S
end

T = [0:0.5:500]
#etas = [0.5:0.1:3.5]
etas = [0.1:0.1:3.5]
#etas = [2.0:0.1:2.0]
S = zeros(Float64, length(etas), 2*length(T))
particle = Bosons(Nparticles, Nparticlemodes, E0, BoxPotential)
for (i, eta1) in enumerate(etas)
    println(eta1)
    mode1 = CavityMode(index1, cavitymodes1, delta1, eta1, kappa1, U0_1)
    system = MultimodeSystem(s, particle, [mode1])
    S[i,:] = spectrum(system, T)
end

writecsv("spectrum_mode12.csv", S)

# using PyCall
# @pyimport matplotlib.pyplot as plt
# plt.plot(S)
# plt.show()
