using MultimodeCavity
const mm = MultimodeCavity

N = 100
mass = 0.5

k = 1.
delta = -1
eta = 1
kappa = 2
U0 = 1

particles = mm.classical.Particles(N, mass)
mode1 = mm.classical.CavityMode(k, delta, eta, kappa, U0)
mode2 = mm.classical.CavityMode(k, delta, eta, kappa, U0)
system = mm.classical.MultimodeSystem(particles, [mode1, mode2])
