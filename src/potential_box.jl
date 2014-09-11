module potentials

export Potential, BoxPotential

type Potential
    A::Function
    B::Function
    E::Function
end

fcos(s, n::Int, i::Int, j::Int) = cospi((i+j)/2)*sinc((i+j+2*n*s)/2)
fsin(s, n::Int, i::Int, j::Int) = sinpi((i+j+n)/2)*sinc((i+j+n*s)/2)

function A(s, n::Int, i::Int, j::Int)
    @assert n>=0
    @assert i>0
    @assert j>0
    # i += 1
    # j += 1
    # n += 1
    k = (i==j ? 0.5 : 0.)
    k + (-1)^n/4 * (fcos(s,n,i,j) + fcos(s,-n,i,j) - fcos(s,n,i,-j) - fcos(s,-n,i,-j))
end

function B(s, n::Int, i::Int, j::Int)
    @assert n>=0
    @assert i>0
    @assert j>0
    # i += 1
    # j += 1
    # n += 1
    return 0.5 * (-fsin(s,n,i,j) + fsin(s,n,-i,j) + fsin(s,n,i,-j) + fsin(s,-n,i,+j))
end

E(E0, i::Int) = E0*i^2

const BoxPotential = Potential(A, B, E)

end