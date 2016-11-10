using QuantumOptics
using multimode
using PyCall
@pyimport matplotlib.pyplot as plt

E0 = 0.2
s = 0.25
Nparticles = 1
Nparticlemodes = 22

index1 = 12
cavitymodes1 = 11
delta1 = -1.
eta1 = 5.0
kappa1 = 1.
U0_1 = 0

index2 = 20
cavitymodes2 = 11
delta2 = -1.
eta2 = 5.0
kappa2 = 1.
U0_2 = 0

particle = Bosons(Nparticles, Nparticlemodes, E0, BoxPotential)
mode1 = CavityMode(index1, cavitymodes1, delta1, eta1, kappa1, U0_1)
mode2 = CavityMode(index2, cavitymodes2, delta2, eta2, kappa2, U0_2)
system = MultimodeSystem(s, particle, [mode1, mode2])

Hp = multimode.Hamiltonian(particle)
B1 = multimode.B(1, system)
B2 = multimode.B(2, system)

ξ1 = -5.
ξ2 = 6.

Hplus = Hp + ξ1*B1 + ξ2*B2
Hminus = Hp - ξ1*B1 - ξ2*B2

Dplus, Vplus = eig(Hplus.data)
Dminus, Vminus = eig(Hminus.data)

println(Dplus)
println(Dminus)

vplus1 = Ket(particle.basis, Vplus[:,1])
vminus1 = Ket(particle.basis, Vminus[:,1])
vplus2 = Ket(particle.basis, Vplus[:,2])
vminus2 = Ket(particle.basis, Vminus[:,2])
vplus3 = Ket(particle.basis, Vplus[:,3])
vminus3 = Ket(particle.basis, Vminus[:,3])
vplus4 = Ket(particle.basis, Vplus[:,4])
vminus4 = Ket(particle.basis, Vminus[:,4])

#vplus = sqrt(0.8)*vplus1 + sqrt(0.15)*vplus2 + sqrt(0.05)*vplus4
#vminus = sqrt(0.8)*vminus1 + sqrt(0.15)*vminus2 + sqrt(0.05)*vminus4
#vplus = sqrt(0.98)*vplus1 + sqrt(0.02)*vplus2
#vminus = sqrt(0.98)*vminus1 + sqrt(0.02)*vminus2
vplus = vplus1
vminus = vminus1

n_xplus = particledensity(particle, [-1:0.01:1], vplus)
n_xminus = particledensity(particle, [-1:0.01:1], vminus)

exp_B1 = dagger(vplus)*B1*vplus
exp_B2 = dagger(vplus)*B2*vplus
η1 = sqrt((delta1^2 + kappa1^2)/(delta1*exp_B1)*ξ1)
η2 = sqrt((delta2^2 + kappa2^2)/(delta2*exp_B2)*ξ2)
α1 = η1*exp_B1/(delta1 + 1im*kappa1)
α2 = η2*exp_B2/(delta2 + 1im*kappa2)

println("<B1> = ", exp_B1)
println("<B2> = ", exp_B2)
println("η1 = ", η1)
println("η2 = ", η2)
println("α1 = ", α1)
println("α2 = ", α2)
println("n1 = ", abs2(α1))
println("n2 = ", abs2(α2))

# for i=4:5
#     vplus = Ket(particle.basis, Vplus[:,i])
#     vminus = Ket(particle.basis, Vminus[:,i])
#     n_xplus = particledensity(particle, [-1:0.01:1], vplus)
#     n_xminus = particledensity(particle, [-1:0.01:1], vminus)
#     plt.plot((n_xplus + n_xminus)/2)
# end

plt.plot(n_xplus)
#plt.plot(n_xminus)
#plt.figure()
#plt.plot((n_xplus + n_xminus)/2)

#plt.show()
#@assert false

mode1.eta = η1
mode2.eta = η2

ρplus₀ = tensor(vplus, dagger(vplus))
ρminus₀ = tensor(vminus, dagger(vminus))
Ψf1plus₀ = coherent_state(mode1.basis, α1)
Ψf1minus₀ = coherent_state(mode1.basis, -α1)
Ψf2plus₀ = coherent_state(mode2.basis, α2)
Ψf2minus₀ = coherent_state(mode2.basis, -α2)
#Ψf₀ = coherent_state(mode1.basis, 0)
#println(expect(destroy(mode1.basis), Ψf₀))
ρf1plus₀ = tensor(Ψf1plus₀, dagger(Ψf1plus₀))
ρf1minus₀ = tensor(Ψf1minus₀, dagger(Ψf1minus₀))
ρf2plus₀ = tensor(Ψf2plus₀, dagger(Ψf2plus₀))
ρf2minus₀ = tensor(Ψf2minus₀, dagger(Ψf2minus₀))
#ρ₀ = 0.5*(ρplus₀ ⊗ ρf1plus₀ ⊗ ρf2plus₀) + 0.5*(ρminus₀ ⊗ ρf1minus₀ ⊗ ρf2minus₀)
ρ₀ = ρplus₀ ⊗ ρf1plus₀ ⊗ ρf2plus₀

#println(expect(destroy(mode1.basis), ρf₀))
#@assert false
#ρ₀ = tensor(basis_ket(system.basis, 1), basis_bra(system.basis, 1))
T, ρ_t = multimode.evolve_master(system, [0.:1:5.], ρ₀)
exp_n1 = expect(number(mode1.basis), ptrace(ρ_t[end], [1,3]))
exp_n2 = expect(number(mode2.basis), ptrace(ρ_t[end], [1,2]))
println("<n1> = ", exp_n1)
println("<n2> = ", exp_n2)


n2 = particledensity(particle, [-1:0.01:1], ptrace(ρ_t[2], [2,3]))
n3 = particledensity(particle, [-1:0.01:1], ptrace(ρ_t[3], [2,3]))
n4 = particledensity(particle, [-1:0.01:1], ptrace(ρ_t[4], [2,3]))
n = particledensity(particle, [-1:0.01:1], ptrace(ρ_t[end], [2,3]))

q1 = qfunc(ptrace(ρ_t[end], [1,3]), [-3:0.1:3], [-3:0.1:3])
q2 = qfunc(ptrace(ρ_t[end], [1,2]), [-3:0.1:3], [-3:0.1:3])

plt.plot(n2)
plt.plot(n3)
plt.plot(n4)
plt.plot(n)
plt.figure()
plt.imshow(q1)
plt.figure()
plt.imshow(q2)
plt.show()
#println([calc_exp_B(ρ) for ρ in ρ_t])
#println(real(diag(ptrace(ρ_t[1], [2]).data)))
#println(real(diag(ptrace(ρ_t[end], [2]).data)))


#println(norm((Hplus*v1-D[1]*v1).data))