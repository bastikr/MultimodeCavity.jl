module multiparticlebasis

export NParticleBasis, BosonicBasis, FermionicBasis, weighted_cdaggeri_cj

using quantumoptics

function distribute_bosons(Nparticles::Int, Nlevels::Int, index::Int=1, occupations::Vector{Int}=zeros(Int,Nlevels), results::Vector{Vector{Int}}=Vector{Int}[])
    if index==Nlevels
        occupations[index] = Nparticles
        push!(results, copy(occupations))
    else
        for n=Nparticles:-1:0
            occupations[index] = n
            distribute_bosons(Nparticles-n, Nlevels, index+1, occupations, results)
        end
    end
    return results
end

function distribute_fermions(Nparticles::Int, Nlevels::Int, index::Int=1, occupations::Vector{Int}=zeros(Int,Nlevels), results::Vector{Vector{Int}}=Vector{Int}[])
    if (Nlevels-index)+1<Nparticles
        return results
    end
    if index==Nlevels
        occupations[index] = Nparticles
        push!(results, copy(occupations))
    else
        for n=min(1,Nparticles):-1:0
            occupations[index] = n
            distribute_fermions(Nparticles-n, Nlevels, index+1, occupations, results)
        end
    end
    return results
end

abstract NParticleBasis <: Basis

type BosonicBasis <: NParticleBasis
    shape::Vector{Int}
    Nparticles::Int
    Nlevels::Int
    occupations::Vector{Vector{Int}}
    BosonicBasis(Nparticles, Nlevels, occupations) = new([length(occupations)], Nparticles, Nlevels, occupations)    
end

type FermionicBasis <: NParticleBasis
    shape::Vector{Int}
    Nparticles::Int
    Nlevels::Int
    occupations::Vector{Vector{Int}}
    FermionicBasis(Nparticles, Nlevels, occupations) = new([length(occupations)], Nparticles, Nlevels, occupations)
end

BosonicBasis(Nparticles::Int, Nlevels::Int) = BosonicBasis(Nparticles, Nlevels, distribute_bosons(Nparticles, Nlevels))
FermionicBasis(Nparticles::Int, Nlevels::Int) = FermionicBasis(Nparticles, Nlevels, distribute_fermions(Nparticles, Nlevels))

=={T<:NParticleBasis}(b1::T, b2::T) = (b1.Nparticles==b2.Nparticles && b1.Nlevels==b2.Nlevels)

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
            for i=1:basis.Nlevels
                x[m,m] += f(i,i)*basis.occupations[m][i]
            end
        end
        idx_create, idx_destroy = cdaggeri_cj(m, n, basis.occupations[m], basis.occupations[n])
        if idx_create!=0 && idx_destroy!=0
            Ni = basis.occupations[n][idx_create]
            Nj = basis.occupations[n][idx_destroy]+1
            #println("Ni=",Ni,";Nj=",Nj,": f(",idx_create,",",idx_destroy,") = ",f(idx_create,idx_destroy))
            x[m,n] += f(idx_create,idx_destroy)*sqrt(Ni*Nj)
        end
    end
    return result
end



end