type PolynomialVariableHeightModel
    g::Float64
    f::Poly{Float64}
    fp::Poly{Float64}
    fpp::Poly{Float64}
    h::Poly{Float64}

    function PolynomialVariableHeightModel(g::Float64, f::Poly{Float64})
        fp = polyder(f)
        fpp = polyder(fp)
        x = Poly([0., 1.])
        h = f - fp * x
        new(g, f, fp, fpp, h)
    end
end

getx(model::PolynomialVariableHeightModel, state::Vector{Float64}) = state[1]
getxd(model::PolynomialVariableHeightModel, state::Vector{Float64}) = state[2]
getz(model::PolynomialVariableHeightModel, state::Vector{Float64}) = polyval(model.f, getx(model, state))
getzd(model::PolynomialVariableHeightModel, state::Vector{Float64}) = polyval(polyder(model.f), getx(model, state)) * getxd(model, state)
getztraj(model::PolynomialVariableHeightModel, x::Float64) = polyval(model.f, x)

function normalized_leg_force(model::PolynomialVariableHeightModel, state::Vector{Float64})
    g = model.g
    x = getx(model, state)
    z = getz(model, state)
    fp = polyval(model.fp, x)
    fpp = polyval(model.fpp, x)
    xd = getxd(model, state)
    force = (g + fpp * xd^2) / (z - fp * x) * [x; z]
end

function odefun(model::PolynomialVariableHeightModel, state::Vector)
    xd = getxd(model, state)
    xdd = normalized_leg_force(model, state)[1]
    [xd; xdd]
end

function orbital_energy(model::PolynomialVariableHeightModel, state::Vector)
    g = model.g
    x = getx(model, state)
    xd = getxd(model, state)
    h = polyval(model.h, x)
    # Eo = 1./2. * xd^2 * h^2 - g * polyval(polyint(h * Poly([0., 1.])), x)
    f = polyval(model.f, x)
    Eo = 1/2 * xd^2 * h^2 + g * x^2 * f - 3 * g * polyval(polyint(model.f * Poly([0., 1.])), x)
    # println(Eo - Eo2)
end

function velocity_given_orbital_energy(model::PolynomialVariableHeightModel, x, Eo)
    g = model.g
    h = polyval(model.h, x)
    f = polyval(model.f, x)
    xd = sqrt((2 * Eo - 2 * g * x^2 * f + 6 * g * polyval(polyint(model.f * Poly([0., 1.])), x)) / h^2)
    # xd = copysign(sqrt(2 * (Eo + g * polyval(polyint(h * Poly([0., 1.])), x)) / h^2), -x)
end
