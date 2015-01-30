module multimode

include("potential_box.jl")
include("system.jl")
include("multiparticlebasis.jl")
include("quantum.jl")

using .potentials
using .system
using .quantum

end # module
