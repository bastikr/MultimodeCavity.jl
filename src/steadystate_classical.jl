module steadystate_classical

export OpticalPotential, evaluate_opticalpotential

using ..system_classical

type OpticalPotential
    N::Int
    k::Vector{Float64}
    ξ::Vector{Float64}
    function OpticalPotential(k::Vector{Float64}, ξ::Vector{Float64})
        @assert length(k)==length(ξ)
        x = new()
        x.N = length(k)
        x.k = k
        x.ξ = ξ
        return x
    end
end

OpticalPotential(s::MultimodeSystem, state::ClassicalState) = OpticalPotential(Float64[mode.k for mode=s.modes], Float64[2*s.modes[n].eta*αn(state,n)[1] for n=1:length(s.modes)])


function evaluate_opticalpotential(U::OpticalPotential, derivative::Int, x)
    f = 0.
    trigfunc = iseven(derivative) ? sin : cos
    sign = iseven(div(derivative, 2)) ? 1. : -1.
    for n=1:U.N
        f += sign*U.ξ[n]*U.k[n]^derivative*trigfunc(U.k[n]*x)
    end
    return f
end

function localextremum(U::OpticalPotential, x0::Float64, eps::Float64=1e-5, dxmax::Float64=1., maxiterations=100)
    function f(x::Float64)
        a = evaluate_opticalpotential(U, 1, x)
        b = evaluate_opticalpotential(U, 2, x)
        abs(b)<0.1 ? NaN : -a/b
    end
    fx = 1
    i = 0
    while abs(fx)>eps
        fx = f(x0)
        i += 1
        if i>maxiterations || isnan(fx)
            return NaN
        end
        x0 = x0 + fx
    end
    return x0
end

function extrema(U::OpticalPotential)
    N = 10*prod(U.k)/(2*pi)^U.N
    eps = 1e-3/N
    solutions = Float64[]
    for x0=linspace(0.,1.,N)
        x = localextremum(U, x0, eps)
        if (!isnan(x)) && (0.<=x<1.) && (!any([abs(y-x)<eps for y=solutions]))
            push!(solutions, x)
        end
    end
    return solutions
end

minima(U::OpticalPotential) = filter(x->evaluate_opticalpotential(U,2,x)>0, extrema(U))
maxima(U::OpticalPotential) = filter(x->evaluate_opticalpotential(U,2,x)<0, extrema(U))
global_minimum(U::OpticalPotential) = (P = extrema(U); length(P)==0 ? NaN : P[indmin(evaluate_opticalpotential(U,0,P))])
global_maximum(U::OpticalPotential) = (P = extrema(U); length(P)==0 ? NaN : P[indmax(evaluate_opticalpotential(U,0,P))])

end #module
