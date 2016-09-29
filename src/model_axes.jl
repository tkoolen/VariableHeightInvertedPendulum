type ModelAxes
    ax
    velocity_scale
    force_scale
    leg
    com
    fgrav
    fgravtext
    fleg
    flegtext
    qdot
    qdottext
    ztraj_sim
    show_gravity
    show_leg_force
    show_qdot

    function ModelAxes(ax, results::SimulationResults, xrange, font_size, restrict_ztraj, show_gravity, show_leg_force)
        ret = new(ax)

        com_radius = 0.05
        xrange_model = [xrange[1] - 1.5 * com_radius; xrange[2] + 1.5 * com_radius]
        model = results.model
        x0 = getx(model, results.state0)
        z0 = getz(model, results.state0)
        zf = getzf(model)
        zref = max(z0, zf)
        zrange_model = [-0.1 * zref; 1.3 * zref]
        # ret.force_scale = (diff(zrange_model)[1] / 4) / max(map(f -> f[2], leg_forces)...)
        ret.force_scale = (diff(zrange_model)[1] / 4) / model.g
        ret.velocity_scale = (diff(xrange_model)[1] / 3)

        plt[:sca](ax)

        # z trajectory
        xs_linspace = restrict_ztraj ? linspace(x0, 0, 100) : linspace(xrange_model[1], xrange_model[2], 100)
        # ztraj = plot(xs_linspace, [getztraj(model, state0, x)::Float64 for x in xs_linspace], "r--", zorder = 7, label = "trajectory")

        # foot
        plot([0], [0], "ko", zorder = 4)

        xlabel(L"x" * " [m]")
        ylabel(L"z" * " [m]")
        plt[:axis]("equal")
        xlim(xrange_model)
        ylim(zrange_model)
        xticks([-0.4; -0.2; 0.; 0.2; 0.4])

        # dynamic objects
        # plt[:sca](ax)
        ret.leg = plot([], [], color = "black")[1]
        ret.com = CoMSymbol(ax, com_radius, 2.; zorder = 4)
        function add_arrow_patch(ax, color)
            arrow = patches.FancyArrowPatch([0.; 0.], [0.; 0.], arrowstyle = "simple", fc = color, ec = color, mutation_scale = 20, zorder = 5)
            ax[:add_patch](arrow)
        end

        if show_gravity
            ret.fgrav = add_arrow_patch(ax, force_color)
            ret.fgravtext = text(0., 0., L"m \mathbf{g}", fontsize = font_size)
            ret.fgravtext[:set_clip_on](true)
        end
        ret.show_gravity = show_gravity

        if show_leg_force
            ret.fleg = add_arrow_patch(ax, force_color)
            ret.flegtext = text(0., 0., L"\mathbf{f}_{\mathrm{gr}}", fontsize = font_size)
            ret.flegtext[:set_clip_on](true)
        end
        ret.show_leg_force = show_leg_force

        ret.qdot = add_arrow_patch(ax, velocity_color)
        ret.qdottext = text(0., 0., L"\dot{\mathbf{q}}_0", fontsize = font_size)
        ret.qdottext[:set_clip_on](true)
        ret.show_qdot = true

        ret.ztraj_sim = plot([], [], color = ztraj_color, zorder = 6, label = "simulation")[1]

        # legend(handles = [ztraj; ztraj_sim], frameon = true, handlelength = 3, handletextpad = 0, borderaxespad = 0.5, fontsize = 12)
        ret
    end
end

function update(model_ax::ModelAxes, results::SimulationResults, state_ind_range::UnitRange{Int64}, current_state_ind::Int64)
    model = results.model
    x = results.xs[current_state_ind]
    xd = results.xds[current_state_ind]
    z = results.zs[current_state_ind]
    zd = results.zds[current_state_ind]

    model_ax.leg[:set_data]([0; x], [0; z])
    set_center(model_ax.com, [x; z])

    if model_ax.show_gravity
        fgrav_text_offset = [-0.1; 0.]
        fgrav_tip = [x; z] + model_ax.force_scale * [0; -model.g]
        model_ax.fgrav[:set_positions]([x; z], fgrav_tip)
        model_ax.fgravtext[:set_position](fgrav_tip + fgrav_text_offset)
    end

    if model_ax.show_leg_force
        fleg_text_offset = [0.03; 0.]
        fleg_tip = model_ax.force_scale * results.leg_forces[current_state_ind]
        model_ax.fleg[:set_positions]([0; 0], fleg_tip)
        model_ax.flegtext[:set_position](fleg_tip + fleg_text_offset)
    end

    if model_ax.show_qdot
        qdot_text_offset = [0.; 0.03]
        qdot_tip = [x; z] + model_ax.velocity_scale * [xd; zd]
        model_ax.qdot[:set_positions]([x; z], qdot_tip)
        model_ax.qdottext[:set_position](qdot_tip + qdot_text_offset)
    else
        model_ax.qdot[:set_visible](false)
        model_ax.qdottext[:set_visible](false)
    end

    model_ax.ztraj_sim[:set_data](results.xs[state_ind_range], results.zs[state_ind_range])
    nothing
end
