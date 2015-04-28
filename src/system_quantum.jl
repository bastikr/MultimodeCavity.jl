module system_quantum

export Bosons, Fermions, CavityMode, MultimodeSystem,
        A, B, basis, particledensity, Hamiltonian, Jump_operators


using quantumoptics
using ..potentials
using ..multiparticlebasis


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

basis(s::Bosons) = BosonicBasis(s.Nparticles, s.Nlevels)
basis(s::Fermions) = FermionicBasis(s.Nparticles, s.Nlevels)
basis(s::CavityMode) = FockBasis(s.Nmodes)
basis(s::MultimodeSystem) = CompositeBasis(basis(s.particles), [basis(mode) for mode=s.modes]...)


# Define interaction between cavity modes and particles

const _cicj = multiparticlebasis.weighted_cdaggeri_cj

A(scaling::Float64, mode_index::Int, p::Particles) = _cicj(basis(p), (i,j)->p.potential.A(scaling, mode_index, i, j))
A(submode::Int, system::MultimodeSystem) = A(system.scaling, system.modes[submode].index, system.particles)

B(scaling::Float64, mode_index::Int, p::Particles) = _cicj(basis(p), (i,j)->p.potential.B(scaling, mode_index, i, j))
B(submode::Int, system::MultimodeSystem) = B(system.scaling, system.modes[submode].index, system.particles)


# Define Hamiltonian

function Hamiltonian(p::Particles)
    b = basis(p)
    N = length(b.occupations)
    energies = zeros(Complex128, N)
    for i=1:N
        occ = b.occupations[i]
        for j=1:length(occ)
            energies[i] += occ[j]*p.potential.E(p.E0, j)
        end
    end
    return Operator(b, diagm(energies))
end

Hamiltonian(mode::CavityMode) = (create(basis(mode))*destroy(basis(mode))*(-mode.delta))

function Hamiltonian(system::MultimodeSystem)
    b = basis(system)
    H = embed(b, 1, Hamiltonian(system.particles))

    for i=1:length(system.modes)
        mode = system.modes[i]
        modebasis = basis(mode)
        An = A(system.scaling, mode.index, system.particles)
        Bn = B(system.scaling, mode.index, system.particles)
        H += embed(b, [1,i+1], [An,number(modebasis)*mode.U0])
        H += embed(b, [1,i+1], [Bn,(destroy(modebasis)+create(modebasis))*mode.eta])
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

function multiparticlebasis.particledensity(system::Particles, x::Vector{Float64}, rho::AbstractOperator)
    @assert basis(system)==rho.basis_l
    @assert basis(system)==rho.basis_r
    basis_functions = Vector{Float64}[]
    for i=1:system.Nlevels
        push!(basis_functions, system.potential.basis_function(i, x))
    end
    return particledensity(basis_functions, rho)
end


end # module
