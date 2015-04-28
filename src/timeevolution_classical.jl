module timeevolution_classical

export timeevolution, timeevolution_stochastic

using quantumoptics
using ..system_classical


function timeevolution(T, S::MultimodeSystem, state0::ClassicalState; fout=nothing)
    T = float(T)
    Nparticles = state0.particlenumber
    Ncavitymodes = state0.cavitymodenumber
    mass = S.particles.mass
    function f(t, y::Vector{Float64}, dy::Vector{Float64})
        x, v, α = splitstate(Nparticles, Ncavitymodes, y)
        dx, dv, dα = splitstate(Nparticles, Ncavitymodes, dy)
        for j=1:Nparticles
            dx[j] = v[j]
            dv[j] = 0
        end
        for n=1:Ncavitymodes
            kn = S.modes[n].k
            αn = α[2*n-1] + 1im*α[2*n]
            ηn = S.modes[n].eta
            U0n = S.modes[n].U0
            Δn = S.modes[n].delta
            κn = S.modes[n].kappa
            ηn_eff = 2*ηn*real(αn)
            U0n_eff = 2*U0n*abs2(αn)
            θn = 0.
            Ωn = 0.
            for j=1:Nparticles
                dv[j] -= 1./mass*kn*cos(kn*x[j])*(U0n_eff*sin(kn*x[j]) + ηn_eff)
                θn += sin(kn*x[j])
                Ωn += sin(kn*x[j])^2
            end
            dαn = (1im*(Δn - U0n*Ωn) - κn)*αn - 1im*ηn*θn
            dα[2*n-1] = real(dαn)
            dα[2*n] = imag(dαn)
        end
    end
    if fout==nothing
        t_out = Float64[]
        state_out = ClassicalState[]
        function fout_(t, y::Vector{Float64})
            push!(t_out, t)
            push!(state_out, ClassicalState(Nparticles, Ncavitymodes, deepcopy(y)))
        end
        quantumoptics.ode_dopri.ode(f, T, state0.data; fout=fout_)
        return t_out, state_out
    else
        return quantumoptics.ode_dopri.ode(f, T, state0.data; fout=(t,y)->fout(t, ClassicalState(Nparticlemodes, Ncavitymodes, y)))
    end
end


function step(f, κ::Vector{Float64}, t0::Float64, dt::Float64, state0::Vector{Float64})
    T, states = quantumoptics.ode_dopri.ode(f, [t0, t0+dt], state0, display_initialvalue=false)
    state1 = states[1]
    return state1
end

function ode_stochastic(f, κ::Vector{Float64}, dt::Float64, T::Vector{Float64}, state0::Vector{Float64}, fout)
    t = T[1]
    tfinal = T[end]
    for tfixed = T[2:end]
        while t<tfixed
            tnext = min(t+dt, tfixed)
            state0 = step(f, κ, t, tnext-t, state0)
            t = tnext
        end
        fout(t, state0)
    end
end

function timeevolution_stochastic(T, S::MultimodeSystem, state0::ClassicalState; fout=nothing, dt=1e-3)
    T = float(T)
    Nparticles = state0.particlenumber
    Ncavitymodes = state0.cavitymodenumber
    mass = S.particles.mass
    function f(t, y::Vector{Float64}, dy::Vector{Float64})
        x, v, α = splitstate(Nparticles, Ncavitymodes, y)
        dx, dv, dα = splitstate(Nparticles, Ncavitymodes, dy)
        for j=1:Nparticles
            dx[j] = v[j]
            dv[j] = 0
        end
        for n=1:Ncavitymodes
            kn = S.modes[n].k
            αn = α[2*n-1] + 1im*α[2*n]
            ηn = S.modes[n].eta
            U0n = S.modes[n].U0
            Δn = S.modes[n].delta
            κn = S.modes[n].kappa
            ηn_eff = 2*ηn*real(αn)
            U0n_eff = 2*U0n*abs2(αn)
            θn = 0.
            Ωn = 0.
            for j=1:Nparticles
                dv[j] -= 1./mass*kn*cos(kn*x[j])*(U0n_eff*sin(kn*x[j]) + ηn_eff)
                θn += sin(kn*x[j])
                Ωn += sin(kn*x[j])^2
            end
            dαn = (1im*(Δn - U0n*Ωn) - κn)*αn - 1im*ηn*θn
            dα[2*n-1] = real(dαn)
            dα[2*n] = imag(dαn)
        end
    end
    κ = Float64[S.modes[n].kappa for n=1:Ncavitymodes]
    if fout==nothing
        t_out = Float64[]
        state_out = ClassicalState[]
        function fout_(t, y::Vector{Float64})
            push!(t_out, t)
            push!(state_out, ClassicalState(Nparticles, Ncavitymodes, deepcopy(y)))
        end
        ode_stochastic(f, κ, dt, T, state0.data, fout_)
        return t_out, state_out
    else
        return ode_stochastic(f, κ, dt, T, state0.data, fout=(t,y)->fout(t, MFState(Nparticlemodes, Ncavitymodes, y)))
    end
end

end # module
