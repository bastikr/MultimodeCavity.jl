module multimode

export potentials, BoxPotential,
        system_quantum,
        system_classical

include("potential_box.jl")
# include("multiparticlebasis.jl")
include("system_quantum.jl")
include("system_classical.jl")

include("timeevolution_quantum.jl")
include("timeevolution_classical.jl")
include("steadystate_classical.jl")


using .potentials
using .system_quantum
using .system_classical
using .timeevolution_quantum
using .timeevolution_classical
using .steadystate_classical

end # module
