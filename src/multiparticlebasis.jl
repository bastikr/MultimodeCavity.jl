module multiparticlebasis

export NParticleBasis, BosonicNParticleBasis, FermionicNParticleBasis, weighted_cdaggeri_cj, cdaggeri_cj, particledensity

using QuantumOptics


function distribute_bosons(particlenumber::Int, singleparticledimension::Int, index::Int=1, occupations::Vector{Int}=zeros(Int,singleparticledimension), results::Vector{Vector{Int}}=Vector{Int}[])
    if index==singleparticledimension
        occupations[index] = particlenumber
        push!(results, copy(occupations))
    else
        for n=particlenumber:-1:0
            occupations[index] = n
            distribute_bosons(particlenumber-n, singleparticledimension, index+1, occupations, results)
        end
    end
    return results
end

function distribute_fermions(particlenumber::Int, singleparticledimension::Int, index::Int=1, occupations::Vector{Int}=zeros(Int,singleparticledimension), results::Vector{Vector{Int}}=Vector{Int}[])
    if (singleparticledimension-index)+1<particlenumber
        return results
    end
    if index==singleparticledimension
        occupations[index] = particlenumber
        push!(results, copy(occupations))
    else
        for n=min(1,particlenumber):-1:0
            occupations[index] = n
            distribute_fermions(particlenumber-n, singleparticledimension, index+1, occupations, results)
        end
    end
    return results
end


abstract NParticleBasis <: Basis

type BosonicNParticleBasis <: NParticleBasis
    shape::Vector{Int}
    particlenumber::Int
    particlebasis::Basis
    occupations::Vector{Vector{Int}}
    function BosonicNParticleBasis(particlenumber, particlebasis, occupations)
        if particlenumber < 0
            throw(ArgumentError("Can't have less than zero particles."))
        end
        for occupation=occupations
            if length(particlebasis) != length(occupation)
                throw(ArgumentError("Dimension of single particle basis has to be equal to the dimension of the N-particle basis vector."))
            end
            if any(occupation.<=0)
                throw(ArgumentError("Occupation numbers smaller than zero not possible."))
            end
            if sum(occupation) != particlenumber
                throw(ArgumentError("Total occupation has to be equal to the particle number."))
            end
        end
        new([length(particlebasis)], particlenumber, particlebasis, occupations)
    end
end

type FermionicNParticleBasis <: NParticleBasis
    shape::Vector{Int}
    particlenumber::Int
    particlebasis::Basis
    occupations::Vector{Vector{Int}}
    function FermionicNParticleBasis(particlenumber::Int, particlebasis::Basis, occupations::Vector{Vector{Int}})
        if particlenumber < 0
            throw(ArgumentError("Can't have less than zero particles."))
        end
        for occupation=occupations
            if length(particlebasis) != length(occupation)
                throw(ArgumentError("Dimension of single particle basis has to be equal to the dimension of the N-particle basis vector."))
            end
            if any(occupation.<=0)
                throw(ArgumentError("Occupation numbers smaller than zero not possible."))
            end
            if sum(occupation) != particlenumber
                throw(ArgumentError("Total occupation has to be equal to the particle number."))
            end
            if any(occupation.>1)
                throw(ArgumentError("Occupation numbers greater than zero not possible for Fermions."))
            end
        end
        new([length(particlebasis)], particlenumber, particlebasis, occupations)
    end
end

BosonicNParticleBasis(particlenumber::Int, particlebasisdimension::Int) = BosonicNParticleBasis(particlenumber, GenericBasis([particlebasisdimension]), distribute_bosons(particlenumber, particlebasisdimension))
FermionicNParticleBasis(particlenumber::Int, particlebasisdimension::Int) = FermionicNParticleBasis(particlenumber, GenericBasis([particlebasisdimension]), distribute_fermions(particlenumber, particlebasisdimension))

BosonicNParticleBasis(particlenumber::Int, particlebasis::Basis) = BosonicNParticleBasis(particlenumber, particlebasis, distribute_bosons(particlenumber, length(particlebasis)))
FermionicNParticleBasis(particlenumber::Int, particlebasis::Basis) = FermionicNParticleBasis(particlenumber, particlebasis, distribute_fermions(particlenumber, length(particlebasis)))


=={T<:NParticleBasis}(b1::T, b2::T) = (b1.particlenumber==b2.particlenumber && b1.particlebasis==b2.particlebasis)

function multiparticleoperator(multiparticlebasis::NParticleBasis, op::Operator, rank::Int)
    subbasis = [multiparticlebasis.particlebasis for i=1:N]
    particlebasisdimension = length(multiparticlebasis.particlebasis)
    for indices=product([particlebasisdimension for i=1:N]...)

    end
end

function cdaggeri_cj(m::Int,n::Int,occ_m::Vector{Int},occ_n::Vector{Int})
    idx_create = 0
    idx_destroy = 0
    for i=1:length(occ_m)
        if occ_m[i]==occ_n[i]
            continue
        end
        delta = occ_n[i] - occ_m[i]
        if delta==-1
            if idx_destroy != 0
                return 0, 0
            end
            idx_destroy = i
            N_destroy = occ_n[i]
        elseif delta==1
            if idx_create != 0
                return 0, 0
            end
            idx_create = i
            N_create = occ_n[i]+1
        else
            return 0, 0
        end
    end
    return idx_create, idx_destroy
end

function weighted_cdaggeri_cj(basis::NParticleBasis, f::Function)
    result = Operator(basis)
    x = result.data
    N = basis.shape[1]
    for m=1:N, n=1:N
        if m==n
            for i=1:basis.singleparticledimension
                x[m,m] += f(i,i)*basis.occupations[m][i]
            end
        end
        idx_create, idx_destroy = cdaggeri_cj(m, n, basis.occupations[m], basis.occupations[n])
        if idx_create!=0 && idx_destroy!=0
            Ni = basis.occupations[n][idx_create]
            Nj = basis.occupations[n][idx_destroy]+1
            x[m,n] += f(idx_create,idx_destroy)*sqrt(Ni*Nj)
        end
    end
    return result
end


function cdaggeri_cj(basis::NParticleBasis, i::Int, j::Int)
    result = Operator(basis)
    N = length(basis.occupations)
    for n=1:N, m=1:N
        occ_n = deepcopy(basis.occupations[n])
        occ_m = deepcopy(basis.occupations[m])
        Ni = occ_n[i]
        Nj = occ_m[j]
        if Ni==0 || Nj==0
            continue
        end
        occ_n[i] -= 1
        occ_m[j] -= 1
        if occ_n == occ_m
            result.data[n,m] = sqrt(Ni*Nj)
        end
    end
    return result
end

function particledensity(basis_functions::Vector, rho)
    N_positionpoints = length(basis_functions[1])
    n = zeros(Float64, N_positionpoints)
    for i=1:N_positionpoints
        f(idx_create::Int, idx_destroy::Int) = basis_functions[idx_create][i]*basis_functions[idx_destroy][i]
        n_op = weighted_cdaggeri_cj(basis(rho), f)
        n[i] = real(expect(n_op, rho))
    end
    return n
end

end