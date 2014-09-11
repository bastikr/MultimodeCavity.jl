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
    E0::Float64
    potential::Potential
    basis::NParticleBasis
end

Particle(V::Potential, basis::NParticleBasis, E0::Float64) = Particle(basis.Nparticles, basis.Nlevels, E0, V, basis)
Bosons(Np::Int, Nl::Int, E0::Float64, V::Potential) = Particle(Np, Nl, E0, V, BosonicBasis(Np, Nl))
Fermions(Np::Int, Nl::Int, E0::Float64, V::Potential) = Particle(Np, Nl, E0, V, FermionicBasis(Np, Nl))


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
    scaling::Float64
    particles::Particle
    modes::Vector{CavityMode}
    basis::CompositeBasis
    MultimodeSystem(scaling, particles, modes) = new(scaling, particles, modes, bases.compose(particles.basis, [m.basis for m=modes]...))
end

const _cicj = multiparticlebasis.weighted_cdaggeri_cj

A(scaling::Float64, mode_index::Int, p::Particle) = _cicj(p.basis, (i,j)->p.potential.A(scaling, mode_index, i, j))
A(submode::Int, system::MultimodeSystem) = A(system.scaling, system.modes[submode].index, system.particles)

B(scaling::Float64, mode_index::Int, p::Particle) = _cicj(p.basis, (i,j)->p.potential.B(scaling, mode_index, i, j))
B(submode::Int, system::MultimodeSystem) = B(system.scaling, system.modes[submode].index, system.particles)


function embed(basis::CompositeBasis, indices::Vector{Int}, operators::Vector)
    op_total = (1 in indices ? operators[findfirst(indices, 1)] : identity(basis.bases[1]))
    for i=2:length(basis.bases)
        op = (i in indices ? operators[findfirst(indices, i)] : identity(basis.bases[i]))
        op_total = tensor(op_total, op)
    end
    return op_total
end

function Hamiltonian(p::Particle)
    N = length(p.basis.occupations)
    energies = zeros(Complex128, N)
    for i=1:N
        occ = p.basis.occupations[i]
        for j=1:length(occ)
            energies[i] += occ[j]*p.potential.E(p.E0, j)
        end
    end
    return Operator(p.basis, diagm(energies))
end

Hamiltonian(mode::CavityMode) = (create(mode.basis)*destroy(mode.basis)*(-mode.delta))



function Hamiltonian(system::MultimodeSystem)
    Hp = Hamiltonian(system.particles)
    H = embed(system.basis, [1], [Hp])

    for i=1:length(system.modes)
        mode = system.modes[i]
        An = A(system.scaling, mode.index, system.particles)
        Bn = B(system.scaling, mode.index, system.particles)
        n = number(mode.basis)
        a = destroy(mode.basis)
        aH = create(mode.basis)
        H += embed(system.basis, [1,i+1], [An,n*mode.U0])
        H += embed(system.basis, [1,i+1], [Bn,(a+aH)*mode.eta])
        H += embed(system.basis, [i+1], [Hamiltonian(mode)])
    end
    return H
end

function JumpOperators(system::MultimodeSystem)
    J = Operator[]
    for i=1:length(system.modes)
        mode = system.modes[i]
        a = destroy(mode.basis)
        j = embed(system.basis, [i+1], [sqrt(mode.kappa)*a])
        push!(J, j)
    end
    return J
end



E0 = 0.2
s = 0.25
Nparticles = 1
Nparticlemodes = 20

index0 = 12
cavitymodes0 = 20
delta0 = -4.
eta0 = 4.5
kappa0 = 1.
U0 = -2.

particle = Bosons(Nparticles, Nparticlemodes, E0, BoxPotential)
mode1 = CavityMode(index0, cavitymodes0, delta0, eta0, kappa0, U0)
system = MultimodeSystem(s, particle, [mode1])

H = Hamiltonian(system)

J = JumpOperators(system)
Jdagger = map(dagger, J)
J_sparse = map(operators_sparse.SparseOperator, J)
Jdagger_sparse = map(operators_sparse.SparseOperator, Jdagger)

#Hnh = H - 0.5im*sum([Jdagger[i]*J[i] for i=1:length(J)])
Hnh = H - 0.5im*Jdagger[1]*J[1]
Hnh_dagger = dagger(Hnh)
Hnh_sparse = operators_sparse.SparseOperator(Hnh)
Hnh_dagger_sparse = operators_sparse.SparseOperator(Hnh_dagger)



Ψ₀ = basis_ket(system.basis, 1)
ρ₀ = Ψ₀⊗dagger(Ψ₀)
T = [0,1]

# tout, ρ = timeevolution.master(T, ρ₀, H, J)
# @time tout, ρ = timeevolution.master(T, ρ₀, H, J)
# tout, ρnh = timeevolution.master_nh(T, ρ₀, Hnh, J; Jdagger=Jdagger, Hdagger=Hnh_dagger)
# @time tout, ρnh = timeevolution.master_nh(T, ρ₀, Hnh, J; Jdagger=Jdagger, Hdagger=Hnh_dagger)
tout, ρnh_sparse = timeevolution.master_nh(T, ρ₀, Hnh_sparse, J_sparse; Jdagger=Jdagger_sparse, Hdagger=Hnh_dagger_sparse)
@time tout, ρnh_sparse = timeevolution.master_nh(T, ρ₀, Hnh_sparse, J_sparse; Jdagger=Jdagger_sparse, Hdagger=Hnh_dagger_sparse)


println("<H> = ", expect(H, ρ[end]))
println("<H> = ", expect(H, ρnh_sparse[end]))
println("<H> = ", expect(H, ρnh_sparse[end]))


end