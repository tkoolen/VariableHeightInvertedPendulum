type SimulationResults
    model
    state0::Vector{Float64}
    ts::Vector{Float64}
    xs::Vector{Float64}
    zs::Vector{Float64}
    xds::Vector{Float64}
    zds::Vector{Float64}
    us::Vector{Float64}
    leg_forces::Vector{Vector{Float64}}
    normalized_leg_force_intensities::Vector{Float64}

    function SimulationResults(model, state0, ts)
        ts, states = ode45((t, state) -> odefun(model, state), state0, ts, points=:specified)
        xs = [getx(model, state)::Float64 for state in states]
        zs = [getz(model, state)::Float64 for state in states]
        xds = [getxd(model, state)::Float64 for state in states]
        zds = [getzd(model, state)::Float64 for state in states]
        us = [leg_u(model, state) for state in states]
        leg_forces = [u * [getx(model, state); getz(model, state)] for (u, state) in zip(us, states)]
        normalized_leg_force_intensities = [dot(f, [x; z] / (model.g * norm([x; z]))) for (x, z, f) in zip(xs, zs, leg_forces)]

        new(model, state0, ts, xs, zs, xds, zds, us, leg_forces, normalized_leg_force_intensities)
    end
end
