module system_classical

export Particles, CavityMode, MultimodeSystem,
        ClassicalState, x, v, αn, splitstate,
        orderparameter, orderparameters

using ArrayViews


# Define system

type Particles
    N::Int
    mass::Float64
end

type CavityMode
    k::Float64
    delta::Float64
    eta::Float64
    kappa::Float64
    U0::Float64
end

type MultimodeSystem
    particles::Particles
    modes::Vector{CavityMode}
end


# Define state of system

type ClassicalState
    particlenumber::Int
    cavitymodenumber::Int
    data::Vector{Float64}
    function ClassicalState(particlenumber::Int, cavitymodenumber::Int, data::Vector{Float64})
        @assert particlenumber > 0
        @assert cavitymodenumber > 0
        @assert length(data)==2*(particlenumber+cavitymodenumber)
        x = new()
        x.particlenumber = particlenumber
        x.cavitymodenumber = cavitymodenumber
        x.data = data
        return x
    end
end

ClassicalState(particlenumber::Int, cavitymodenumber::Int) = ClassicalState(particlenumber::Int, cavitymodenumber::Int, zeros(Float64, 2*(particlenumber + cavitymodenumber)))
ClassicalState(s::MultimodeSystem, data::Vector{Float64}) = ClassicalState(s.particles.N, length(s.modes), data)
ClassicalState(s::MultimodeSystem) = ClassicalState(s.particles.N, length(s.modes))

x(state::ClassicalState) = view(state.data, 1:state.particlenumber)
v(state::ClassicalState) = view(state.data, state.particlenumber+1:2*state.particlenumber)
αn(state::ClassicalState, n::Int) = view(state.data, 2*state.particlenumber+2*(n-1)+1:2*state.particlenumber+2*n)

splitstate(N::Int, M::Int, data::Vector{Float64}) = view(data, 1:N), view(data, N+1:2*N), view(data, 2*N+1:2*(N+M))

orderparameter(k::Float64, state::ClassicalState) = 1./state.particlenumber*sum(sin(k*x(state)))
orderparameters(K::Vector{Float64}, state::ClassicalState) = Float64[orderparameter(k, state) for k=K]
orderparameters(s::MultimodeSystem, state::ClassicalState) = orderparameters(Float64[mode.k for mode=s.modes], state)


end # module
