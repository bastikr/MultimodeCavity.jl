module multimode_steadystate

using quantumoptics
using ..system
using ..quantum

function Halpha(system::MultimodeSystem, alpha::Number)
    H = Hamiltonian(system.particles)
    for n=1:length(system.modes)
        An = A(n, system)
        Bn = B(n, system)
        mode = system.modes[n]
        #H += 2*mode.eta*real(alpha)*Bn
        H += mode.eta^2*Bn^2*mode.delta/(mode.delta^2 + mode.kappa^2)
    end
    return H
end

function Hfunctional(system::MultimodeSystem, psi::Ket)
    E = expect(Hamiltonian(system.particles), psi)
    for n=1:length(system.modes)
        An = A(n, system)
        Bn = expect(B(n, system), psi)
        mode = system.modes[n]
        #H += 2*mode.eta*real(alpha)*Bn
        E += mode.eta^2*Bn^2*mode.delta/(mode.delta^2 + mode.kappa^2)
    end
    println("E = ", real(E))
    return real(E)
end

function Hξ(system::MultimodeSystem, ξ::Float64)
    H = Hamiltonian(system.particles)
    for n=1:length(system.modes)
        Bn = B(n, system)
        H += ξ*Bn
    end
    return H
end

# function steady_particle(system::MultimodeSystem, alpha<:Number...)
#     Hplus = Hamiltonian(system.particles)
#     Hminus = Hamiltonian(system.particles)
#     for n=1:length(system.modes)
#         An = A(n, system)
#         Bn = B(n, system)
#         mode = system.modes[n]
#         Hplus += mode.U0*abs2(alpha[n]*An) + 2*mode.eta*real(alpha[n])*Bn
#         Hminus += mode.U0*abs2(alpha[n]*An) - 2*mode.eta*real(alpha[n])*Bn
#     end

# end

# function steady_modes(system::MultimodeSystem, xplus::Ket)
#     alpha = Complex128[]
#     for n=1:length(system.modes)
#         a = dagger(xplus)*A(n, system)*xplus
#         b = dagger(xplus)*B(n, system)*xplus
#         push(alpha, b0)
#     end
# end

end #multimode_steadystate