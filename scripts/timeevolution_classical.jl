#!/usr/bin/env julia
using ArgParse

s = ArgParseSettings(
        description = "Calculate classical self-organization model."
)

@add_arg_table s begin
    "--Nparticles"
        arg_type = Int
        required = true
    "--mass"
        arg_type = Float64
        required = true
    "--indices"
        required = true
    "--deltas"
        required = true
    "--etas"
        required = true
    "--kappas"
        required = true
    "--U0s"
        required = true
    "--T"
        required = true
    "--seed"
        arg_type = Int
        required = true
    "--o"
        required = true
end

parameters = parse_args(s)

# Output
const odir = parameters["o"]
@assert isdir(odir)

# Integration
const seed = parameters["seed"]
const T = float(eval(parse(parameters["T"])))
@assert T[1]<T[end]

# Particles
const Nparticles = parameters["Nparticles"]
const particlemass = parameters["mass"]

# Modes
const indices = float(eval(parse(parameters["indices"])))
const Nmodes = length(indices)

const deltas = float(eval(parse(parameters["deltas"])))
const etas = float(eval(parse(parameters["etas"])))
const kappas = float(eval(parse(parameters["kappas"])))
const U0s = float(eval(parse(parameters["U0s"])))

@assert length(deltas)==Nmodes
@assert length(etas)==Nmodes
@assert length(kappas)==Nmodes
@assert length(U0s)==Nmodes

using multimode

const particles = multimode.classical.Particles(Nparticles, particlemass)
const modes = multimode.classical.CavityMode[]

for i=1:Nmodes
    push!(modes, multimode.classical.CavityMode(
                    float(indices[i]*pi),
                    float(deltas[i]),
                    float(etas[i]),
                    float(kappas[i]),
                    float(U0s[i])))
end

const system = multimode.classical.MultimodeSystem(particles, modes)

keyparameters = Dict(
    "N"=>parameters["Nparticles"],
    # "m"=>parameters["mass"],
    "index"=>parameters["indices"],
    "eta"=>parameters["etas"],
    # "delta"=>parameters["deltas"],
    # "kappa"=>parameters["kappas"],
    # "U0"=>parameters["U0s"],
    # "T"=>T[end],
    )

function dict2filename(d::Dict; extension=".dat")
    pairs = ["$(key)=$(item)" for (key, item) in d]
    return join(sort(pairs), ";")*extension
end
outname = "timeevolution_"*dict2filename(keyparameters)
outpath = joinpath(odir, outname)


# # Modes
# const index1 = 4
# const k1 = float(index1*pi)
# const delta1 = -10.
# const eta1 = 50/sqrt(Nparticles)
# const kappa1 = 10.0
# const U0_1 = 0.0/sqrt(Nparticles)
# const mode1 = multimode.classical.CavityMode(k1, delta1, eta1, kappa1, U0_1)

# const index2 = 6
# const k2 = float(index2*pi)
# const delta2 = -10.
# #const eta2 = 36.07135/sqrt(Nparticles)
# const eta2 = 34/sqrt(Nparticles)
# const kappa2 = 10.0
# const U0_2 = 0.0/sqrt(Nparticles)
# const mode2 = multimode.classical.CavityMode(k2, delta2, eta2, kappa2, U0_2)


# const system = multimode.classical.MultimodeSystem(particles, [mode1, mode2])
# const targetpath = "data/timeevolution_modes_$(index1)_$(index2).dat"

# Initial state
state0 = multimode.classical.ClassicalState(Nparticles, length(system.modes))
x0 = multimode.classical.x(state0)
srand(seed)
for j=1:Nparticles
    x0[j] = rand()
end

# Time evolution
tout, states = multimode.classical.timeevolution(T, system, state0)

# Write results to file
f = open(outpath, "w")
write(f, "# System: ")
write(f, string(system))
write(f, "\n# t state")
for i=1:length(T)
    write(f, "\n$(T[i]);")
    write(f, string(states[i]))
end
close(f)
