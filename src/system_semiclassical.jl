module system_semiclassical

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

end # module
