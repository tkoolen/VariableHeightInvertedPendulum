type ControlledVariableHeightInvertedPendulum{C}
    name::String
    g::Float64
    zf::Float64
    controller::C
end

getx(model::ControlledVariableHeightInvertedPendulum, state) = state[1]
getz(model::ControlledVariableHeightInvertedPendulum, state) = state[2]
getxd(model::ControlledVariableHeightInvertedPendulum, state) = state[3]
getzd(model::ControlledVariableHeightInvertedPendulum, state) = state[4]
getzf(model::ControlledVariableHeightInvertedPendulum) = model.zf
unpack(model::ControlledVariableHeightInvertedPendulum, state) = (getx(model, state), getz(model, state), getxd(model, state), getzd(model, state))

function cubic_orbital_energy_controller(model::ControlledVariableHeightInvertedPendulum, state)
    x, z, xd, zd = unpack(model, state)
    g = model.g
    zf = model.zf
    a = xd / x
    b = zd - a * z
    u = -7 * a^2 + (3 * zf * a^3 - g * a) / b - (10 * a^3 * b) / g
end

function cubic_clipped_controller(model::ControlledVariableHeightInvertedPendulum, state)
    u = max(cubic_orbital_energy_controller(model, state), 0.)
end

function odefun(model::ControlledVariableHeightInvertedPendulum, state)
    u = model.controller(model, state)
    x, z, xd, zd = unpack(model, state)
    xdd = getx(model, state) * u
    zdd = getz(model, state) * u - model.g
    [xd; zd; xdd; zdd]
end

# function getztraj(model::CubicControlledVariableHeightInvertedPendulum, state0, x)
#     x0, z0, xd0, zd0 = unpack(model, state0)
#     g = model.g
#     zf = model.zf
#     z = (1/2).*g.^(-1).*x0.^(-5).*xd0.^(-1).*(10.*x.*(x+(-1).*x0).^2.*
#     xd0.*(xd0.*z0+(-1).*x0.*zd0).^2+g.*x0.^2.*((-1).*x.*(5.*x.^2+(-8)
#     .*x.*x0+x0.^2).*xd0.*z0+x.*(5.*x+(-3).*x0).*(x+(-1).*x0).*x0.*zd0+
#     (x+(-1).*x0).^2.*((-5).*x+2.*x0).*xd0.*zf))
# end

leg_u(model::ControlledVariableHeightInvertedPendulum, state::Vector{Float64}) = model.controller(model, state)
