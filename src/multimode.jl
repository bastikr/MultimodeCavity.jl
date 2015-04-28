module multimode

module quantum
include("potential_box.jl")
include("multiparticlebasis.jl")
include("system_quantum.jl")
include("timeevolution_quantum.jl")
using .potentials
using .system_quantum
using .timeevolution_quantum
end # quantum

module classical
include("system_classical.jl")
include("timeevolution_classical.jl")
include("steadystate_classical.jl")
using .system_classical
using .timeevolution_classical
using .steadystate_classical
end # classical

end # module
