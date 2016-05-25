type LIPM
    g::Float64
    z0::Float64
    ω::Float64
    LIPM(g, z0) = new(g, z0, getomega(g, z0))
end

getx(lipm::LIPM, state::Vector{Float64}) = state[1]
getxd(lipm::LIPM, state::Vector{Float64}) = state[2]
getz(lipm::LIPM, state::Vector{Float64}) = lipm.z0
getzd(lipm::LIPM, state::Vector{Float64}) = 0.
getztraj(lipm::LIPM, x::Float64) = lipm.z0
normalized_leg_force(lipm::LIPM, state::Vector{Float64}) = lipm.g * [getx(lipm, state) / lipm.z0; 1]

function odefun(lipm::LIPM, state::Vector)
    x = getx(lipm, state)
    xd = getxd(lipm, state)
    xdd = lipm.ω^2 * x
    [xd; xdd]
end

function orbital_energy(lipm::LIPM, state::Vector)
    x = getx(lipm, state)
    xd = getxd(lipm, state)
    1/2 * xd^2 - lipm.g / (2 * lipm.z0) * x^2
end

function velocity_given_orbital_energy(lipm::LIPM, x, Eo)
    xd = copysign(sqrt(2 * Eo + lipm.g / lipm.z0 * x^2), -x)
end
