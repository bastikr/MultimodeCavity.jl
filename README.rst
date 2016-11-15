MultimodeCavity.jl
==================

**MultimodeCavity.jl** is a numerical framework written in `Julia <http://julialang.org/>`_ used to simulate particles coupled to a multi-mode cavity.

Example
-------

.. code-block:: julia

    using CollectiveSpins

    using QuantumOptics
    using MultimodeCavity
    const mm = MultimodeCavity
    using Optim

    # Particles
    const E0 = 0.2
    const Nparticles = 1
    const Nlevels = 10
    const potential = mm.quantum.BoxPotential
    const particles = mm.quantum.Bosons(Nparticles, Nlevels, E0, potential)

    # Cavity mode
    const index1 = 2
    const cavitymodes1 = 12
    const delta1 = -1.
    const eta1 = 4.0
    const kappa1 = 1.
    const U0_1 = 0
    const mode1 = mm.quantum.CavityMode(index1, cavitymodes1, delta1, eta1, kappa1, U0_1)

    # Multimode system
    const scaling = 1.
    const system = mm.quantum.MultimodeSystem(scaling, particles, [mode1])

    # Bases
    basis_total = mm.quantum.basis(system)
    basis_particles = mm.quantum.basis(particles)
    basis_mode1 = mm.quantum.basis(mode1)

    # Initial state
    const ψ₀ = basis_ket(basis_particles, 1) ⊗ basis_ket(basis_mode1, 1)
    const ρ₀ = ψ₀ ⊗ dagger(ψ₀)

    T = [0.,10.0]

    tout, ρt = mm.timeevolution_master(system, T, ρ₀, reltol=1e-7)


Documentation
-------------

The documentation written with `Sphinx <http://www.sphinx-doc.org/>`_ using the `Sphinx-Julia <https://github.com/bastikr/sphinx-julia>`_ plugin is available at

    https://bastikr.github.io/MultimodeCavity.jl/
