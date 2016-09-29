type PolynomialVariableHeightModel
    name::String
    g::Float64
    f::Poly{Float64}
    fp::Poly{Float64}
    fpp::Poly{Float64}
    h::Poly{Float64}

    function PolynomialVariableHeightModel(name, g::Float64, f::Poly{Float64})
        fp = polyder(f)
        fpp = polyder(fp)
        x = Poly([0., 1.])
        h = f - fp * x
        new(name, g, f, fp, fpp, h)
    end
end

getx(model::PolynomialVariableHeightModel, state::Vector{Float64}) = state[1]
getxd(model::PolynomialVariableHeightModel, state::Vector{Float64}) = state[2]
getz(model::PolynomialVariableHeightModel, state::Vector{Float64}) = polyval(model.f, getx(model, state))
getzd(model::PolynomialVariableHeightModel, state::Vector{Float64}) = polyval(polyder(model.f), getx(model, state)) * getxd(model, state)
getztraj(model::PolynomialVariableHeightModel, state0, x::Float64) = polyval(model.f, x)
getzf(model::PolynomialVariableHeightModel) = polyval(model.f, 0.)

function leg_u(model::PolynomialVariableHeightModel, state::Vector{Float64})
    g = model.g
    x = getx(model, state)
    z = getz(model, state)
    fp = polyval(model.fp, x)
    fpp = polyval(model.fpp, x)
    xd = getxd(model, state)
    u = (g + fpp * xd^2) / (z - fp * x)
end

function odefun(model::PolynomialVariableHeightModel, state::Vector)
    xd = getxd(model, state)
    xdd = leg_u(model, state) * getx(model, state)
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

function zero_orbital_energy_trajectory{T}(x0::T, z0::T, xd0::T, zd0::T, g::T, zf::T)
    Axf = ForwardDiff.jacobian(coefs -> [polyval(Poly(coefs), 0)], ones(T, 4))
    Ax0 = ForwardDiff.jacobian(coefs -> [polyval(Poly(coefs), x0)], ones(T, 4))
    Adzdx0 = ForwardDiff.jacobian(coefs -> [polyval(polyder(Poly(coefs)), x0)], ones(T, 4))
    Aorbital = zeros(T, 1, 4)
    for i = 0 : 3
        Aorbital[i + 1] = 2 * g * (1 - i) / (i + 2) * x0^(i + 2); # orbital energy
    end
    A = [Axf; Ax0; Adzdx0; Aorbital]

    b = [zf; z0; zd0 / xd0; (xd0 * z0 - zd0 * x0)^2]

    coefs = A \ b
    f = Poly(coefs)
end

function is_force_always_nonnegative_condition2(g, x0, z0, xd0, zd0, zf)
    # extracted out so that it can be used in contourf
    # should nonnegative to be in valid region
    a = xd0 ./ x0
    b = zd0 - a .* z0
    ((20.*b+7.*a.^(-1).*g)-(3.^(1/2).*(g.*(3.*a.^(-2).*g+40.*zf)).^(1/2)))
end

function is_force_always_nonnegative(g, x0, z0, xd0, zd0, zf)
    opposite_velocity = xd0 ./ x0
    return opposite_velocity < 0 & is_force_always_nonnegative_condition2(g, x0, z0, xd0, zd0, zf) >= 0
end
