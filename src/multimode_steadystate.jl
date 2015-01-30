module multimode_steadystate

using quantumoptics
using multimode

function steady_particle(system::MultimodeSystem, alpha<:Number...)
    Hplus = Hamiltonian(system.particles)
    Hminus = Hamiltonian(system.particles)
    for n=1:length(system.modes)
        An = A(n, system)
        Bn = B(n, system)
        mode = system.modes[n]
        Hplus += mode.U0*abs2(alpha[n]*An) + 2*mode.eta*real(alpha[n])*Bn
        Hminus += mode.U0*abs2(alpha[n]*An) - 2*mode.eta*real(alpha[n])*Bn
    end

end

function steady_modes(system::MultimodeSystem, xplus::Ket)
    alpha = Complex128[]
    for n=1:length(system.modes)
        a = dagger(xplus)*A(n, system)*xplus
        b = dagger(xplus)*B(n, system)*xplus
        push(alpha, b0)
    end
end

end #multimode_steadystate