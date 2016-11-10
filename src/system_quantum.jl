module system_quantum

export Bosons, Fermions, CavityMode, MultimodeSystem,
        An, Bn, basis, particledensity, Hamiltonian, Jump_operators


using QuantumOptics
using ..potentials
# using ..multiparticlebasis


# Define system

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


# Define basis of Hilbert space

basis(s::Bosons) = BosonicNParticleBasis(s.Nparticles, s.Nlevels)
basis(s::Fermions) = FermionicNParticleBasis(s.Nparticles, s.Nlevels)
basis(s::CavityMode) = FockBasis(s.Nmodes)
basis(s::MultimodeSystem) = CompositeBasis(basis(s.particles), [basis(mode) for mode=s.modes]...)


# Define interaction between cavity modes and particles

function An(scaling, n, p)
    M = [p.potential.A(scaling, mode_index, i, j) for i=1:p.Nlevels, j=1:p.Nlevels]
    op = DenseOperator(GenericBasis(p.Nlevels), complex(M))
    nparticleoperator_1(basis(p), op)
end

function Bn(scaling, n, p)
    M = [p.potential.B(scaling, n, i, j) for i=1:p.Nlevels, j=1:p.Nlevels]
    op = DenseOperator(GenericBasis(p.Nlevels), complex(M))
    nparticleoperator_1(basis(p), op)
end

An(submode::Int, system::MultimodeSystem) = An(system.scaling, system.modes[submode].index, system.particles)
Bn(submode::Int, system::MultimodeSystem) = Bn(system.scaling, system.modes[submode].index, system.particles)


# Define Hamiltonian
function Hamiltonian(p::Particles)
    M = [p.potential.E(p.E0, j) for j=1:p.Nlevels]
    op = DenseOperator(GenericBasis(p.Nlevels), spdiagm(complex(M)))
    nparticleoperator_1(basis(p), op)
end

Hamiltonian(mode::CavityMode) = (create(basis(mode))*destroy(basis(mode))*(-mode.delta))

function Hamiltonian(system::MultimodeSystem)
    b = basis(system)
    H = embed(b, 1, Hamiltonian(system.particles))

    for i=1:length(system.modes)
        mode = system.modes[i]
        modebasis = basis(mode)
        An_ = An(system.scaling, mode.index, system.particles)
        Bn_ = Bn(system.scaling, mode.index, system.particles)
        H += embed(b, [1,i+1], [An_,number(modebasis)*mode.U0])
        H += embed(b, [1,i+1], [Bn_,(destroy(modebasis)+create(modebasis))*mode.eta])
        H += embed(b, i+1, Hamiltonian(mode))
    end
    return H
end

function Jump_operators(system::MultimodeSystem)
    b = basis(system)
    J = Operator[]
    for i=1:length(system.modes)
        mode = system.modes[i]
        a = destroy(basis(mode))
        j = embed(b, i+1, sqrt(2*mode.kappa)*a)
        push!(J, j)
    end
    return J
end

# function multiparticlebasis.particledensity(system::Particles, x::Vector{Float64}, rho::Operator)
#     @assert basis(system)==rho.basis_l
#     @assert basis(system)==rho.basis_r
#     basis_functions = Vector{Float64}[]
#     for i=1:system.Nlevels
#         push!(basis_functions, system.potential.basis_function(i, x))
#     end
#     return particledensity(basis_functions, rho)
# end


end # module
