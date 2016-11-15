module timeevolution_quantum

export timeevolution_master, timeevolution_mcwf

using QuantumOptics
using ..quantum


function timeevolution_master(system::MultimodeSystem, T, ρ₀; kwargs...)
    H = Hamiltonian(system)
    Hsparse = SparseOperator(H)
    J = Jump_operators(system)
    Jsparse = map(SparseOperator, J)
    return timeevolution.master(T, ρ₀, Hsparse, Jsparse; kwargs...)
end

function timeevolution_mcwf(system::MultimodeSystem, T, Ψ₀::Ket; kwargs...)
    H = Hamiltonian(system)
    Hsparse = SparseOperator(H)
    J = Jump_operators(system)
    Jsparse = map(SparseOperator, J)
    return timeevolution.mcwf(T, Ψ₀, Hsparse, Jsparse; kwargs...)
end

end # module