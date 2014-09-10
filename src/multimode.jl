module multimodecavity

export MultimodeSystem,
    multiparticlebasis, BosonicBasis, FermionicBasis,
    potentials, BoxPotential

using quantumoptics

include("multiparticlebasis.jl")
include("potential_box.jl")

using .multiparticlebasis
using .potentials


type Particle
    Nparticles::Int
    Nlevels::Int
    potential::Potential
    basis::NParticleBasis
end

Particle(potential::Potential, basis::NParticleBasis) = Particle(basis.Nparticles, basis.Nlevels, potential, basis)
Bosons(Nparticles::Int, Nlevels::Int, potential::Potential) = Particle(Nparticles, Nlevels, potential, BosonicBasis(Nparticles, Nlevels))
Fermions(Nparticles::Int, Nlevels::Int, potential::Potential) = Particle(Nparticles, Nlevels, potential, FermionicBasis(Nparticles, Nlevels))


type CavityMode
    index::Int
    Nmodes::Int
    delta::Float64
    eta::Float64
    kappa::Float64
    U0::Float64
    basis::FockBasis
end

CavityMode(index, Nmodes, delta, eta, kappa, U0) = CavityMode(index, Nmodes, delta, eta, kappa, U0, FockBasis(Nmodes))


type MultimodeSystem
    particles::Particle
    modes::Vector{CavityMode}
    scaling::Float64
end

A(n::Int, scaling::Float64, p::Particle) = multiparticlebasis.weighted_cdaggeri_cj(p.basis, (i,j)->p.potential.A(scaling, n, i, j))
A(n::Int, system::MultimodeSystem) = A(n, system.scaling, system.particles)

B(n::Int, scaling::Float64, p::Particle) = multiparticlebasis.weighted_cdaggeri_cj(p.basis, (i,j)->p.potential.B(scaling, n, i, j))
B(n::Int, system::MultimodeSystem) = B(n, system.scaling, system.particles)

end