getomega(g, z) = √(g / z)

height_trajectory(model, state0, x)

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


function zero_orbital_energy_cubic{T}(x0::T, z0::T, xd0::T, zd0::T, g::T, zf::T)
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
    a0 = xd0 ./ x0
    b0 = zd0 - a .* z0
    sqrt(9 * g^2 + 120 * a0.^2 * g * zf) - (7 * g + 20 * a0 .* b0)
end

function is_force_always_nonnegative(g, x0, z0, xd0, zd0, zf)
    a0 = xd0 ./ x0
    return a0 < 0 & is_force_always_nonnegative_condition2(g, x0, z0, xd0, zd0, zf) >= 0
end

function cubic_control_law(g, zf, x, z, xd, zd)
    a = xd / x
    b = zd - a * z
    -7 * a^2 + (3 * zf * a^3 - g * a) / b - (10 * a^3 * b) / g
end

# function leg_u(model::PolynomialVariableHeightModel, state::Vector{Float64})
#     g = model.g
#     x = getx(model, state)
#     z = getz(model, state)
#     fp = polyval(model.fp, x)
#     fpp = polyval(model.fpp, x)
#     xd = getxd(model, state)
#     u = (g + fpp * xd^2) / (z - fp * x)
# end
