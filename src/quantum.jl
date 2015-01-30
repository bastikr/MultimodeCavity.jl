module quantum

export A, B, basis, Hamiltonian, Jump_operators, timeevolution_master, timeevolution_mcwf

using quantumoptics
using ..system
using ..multiparticlebasis

basis(s::Bosons) = BosonicBasis(s.Nparticles, s.Nlevels)
basis(s::Fermions) = FermionicBasis(s.Nparticles, s.Nlevels)
basis(s::CavityMode) = FockBasis(s.Nmodes)
basis(s::MultimodeSystem) = CompositeBasis(basis(s.particles), [basis(mode) for mode=s.modes]...)

const _cicj = multiparticlebasis.weighted_cdaggeri_cj

A(scaling::Float64, mode_index::Int, p::Particles) = _cicj(basis(p), (i,j)->p.potential.A(scaling, mode_index, i, j))
A(submode::Int, system::MultimodeSystem) = A(system.scaling, system.modes[submode].index, system.particles)

B(scaling::Float64, mode_index::Int, p::Particles) = _cicj(basis(p), (i,j)->p.potential.B(scaling, mode_index, i, j))
B(submode::Int, system::MultimodeSystem) = B(system.scaling, system.modes[submode].index, system.particles)

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

function particledensity(system::Particles, x::Vector{Float64}, rho)
    basis_functions = Vector{Float64}[]
    for i=1:sys.Nlevels
        push!(basis_functions, sys.potential.basis_function(i, x))
    end
    return particledensity(basis_functions, rho)
end

# function particle_density(system, rho_p, resolution=200)
#     X = linspace(-1, 1, resolution)
#     Nlevels = system.particles.Nlevels
#     basis_functions = [system.particles.potential.basis_function(i, X) for i=1:Nlevels]
#     nx = zeros(Float64, resolution)
#     for i=1:Nlevels, j=1:Nlevels
#         c_ij = cdaggeri_cj(system.particles.basis, i, j)
#         vi = basis_functions[i]
#         vj = basis_functions[j]
#         nx += real(sum(diag((c_ij*rho_p).data)))*conj(vi).*vj
#     end
#     return nx
# end

function timeevolution_master(system::MultimodeSystem, T, ρ₀; kwargs...)
    H = Hamiltonian(system)
    Hsparse = SparseOperator(H)
    J = Jump_operators(system)
    Jsparse = map(SparseOperator, J)
    return timeevolution.master(T, ρ₀, Hsparse, Jsparse; kwargs...)
end

function timeevolution_mcwf(system::MultimodeSystem, T, Ψ₀::Ket; kwargs...)
    H = Hamiltonian(system)
    Hsparse = SparseOperator(H)
    J = Jump_operators(system)
    Jsparse = map(SparseOperator, J)
    return timeevolution.mcwf(T, Ψ₀, Hsparse, Jsparse; kwargs...)
end

end # module
