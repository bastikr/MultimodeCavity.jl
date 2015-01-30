module system

export Particles, Bosons, Fermions, CavityMode, MultimodeSystem

using ..potentials

abstract Particles

type Bosons <: Particles
    Nparticles::Int
    Nlevels::Int
    E0::Float64
    potential::Potential
end

type Fermions <: Particles
    Nparticles::Int
    Nlevels::Int
    E0::Float64
    potential::Potential
end

type CavityMode
    index::Int
    Nmodes::Int
    delta::Float64
    eta::Float64
    kappa::Float64
    U0::Float64
end

type MultimodeSystem
    scaling::Float64
    particles::Particles
    modes::Vector{CavityMode}
end

end # module
