type LIPM
    name::ASCIIString
    g::Float64
    z0::Float64
    ω::Float64
    LIPM(name, g, z0) = new(name, g, z0, getomega(g, z0))
end

getx(model::LIPM, state::Vector{Float64}) = state[1]
getxd(model::LIPM, state::Vector{Float64}) = state[2]
getz(model::LIPM, state::Vector{Float64}) = model.z0
getzd(model::LIPM, state::Vector{Float64}) = 0.
getztraj(model::LIPM, x::Float64) = model.z0
leg_u(model::LIPM, state::Vector{Float64}) = model.g / model.z0
normalized_leg_force(model::LIPM, state::Vector{Float64}) = leg_u(model, state) * [getx(model, state); getz(model, state)]

function odefun(model::LIPM, state::Vector)
    x = getx(model, state)
    xd = getxd(model, state)
    xdd = model.ω^2 * x
    [xd; xdd]
end

function orbital_energy(model::LIPM, state::Vector)
    x = getx(model, state)
    xd = getxd(model, state)
    1/2 * xd^2 - model.g / (2 * model.z0) * x^2
end

function velocity_given_orbital_energy(model::LIPM, x, Eo)
    xd = copysign(sqrt(2 * Eo + model.g / model.z0 * x^2), -x)
end
