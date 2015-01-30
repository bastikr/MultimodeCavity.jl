using multimode

function steady_particle(system::MultimodeSystem, alphas::Vector)
    H = Hamiltonian(system.particles)
    for n=1:length(system.modes)
        An = A(n, system)
        Bn = B(n, system)
        H += system.U0*abs2(alphas[n])
    end
end