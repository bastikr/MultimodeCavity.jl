using MultimodeCavity
const mm = MultimodeCavity

Nparticles = 2
Nlevels = 5
E0 = 0.1
P = mm.potentials.BoxPotential

index = 3
Nmodes = 10
delta = -1
eta = 1
kappa = 2
U0 = 1

s = 0.4

particles_bosons = mm.quantum.Bosons(Nparticles, Nlevels, E0, P)
particles_fermions = mm.quantum.Fermions(Nparticles, Nlevels, E0, P)
mode1 = mm.quantum.CavityMode(index, Nmodes, delta, eta, kappa, U0)
mode2 = mm.quantum.CavityMode(index, Nmodes, delta, eta, kappa, U0)
system = mm.quantum.MultimodeSystem(s, particles_bosons, [mode1, mode2])
system = mm.quantum.MultimodeSystem(s, particles_fermions, [mode1, mode2])

b = mm.quantum.basis(particles_bosons)
b = mm.quantum.basis(particles_fermions)
b = mm.quantum.basis(mode1)
b = mm.quantum.basis(system)

An = mm.quantum.An(1, system)
An = mm.quantum.An(2, system)

Bn = mm.quantum.Bn(1, system)
Bn = mm.quantum.Bn(2, system)

H = mm.quantum.Hamiltonian(particles_fermions)
H = mm.quantum.Hamiltonian(particles_bosons)
H = mm.quantum.Hamiltonian(mode1)
H = mm.quantum.Hamiltonian(system)

J = mm.quantum.Jump_operators(system)
