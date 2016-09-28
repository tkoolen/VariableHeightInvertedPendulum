function sim_movie(model, state0, filename::ASCIIString;
    staticplot = false, fps = 30., realtimerate = 0.25, dpi = 100, simtime = 1., stilltime = 0.5,
    restrict_ztraj = false, show_region = false, model_only = false, show_orbital_energy = false, font_size = 15, fig_size = (8., 6.))

    ztraj_color = "#364b9f" # blue
    force_color = "#f37620" # orange
    velocity_color = "#be2025" # red

    ioff()
    fig = figure(figsize = fig_size)
    plt[:rc]("font", size = font_size, family = "serif")

    if model_only
        ss_ax = nothing
        force_ax = nothing
        model_ax = gca()
    else
        ss_ax = plt[:subplot2grid]((2,2), (0,0))
        force_ax = plt[:subplot2grid]((2,2), (1,0))
        model_ax = plt[:subplot2grid]((2,2), (0, 1), rowspan=2)
    end

    x0 = getx(model, state0)
    xd0 = getxd(model, state0)
    z0 = getz(model, state0)
    zd0 = getzd(model, state0)
    dzdx0 = zd0 / xd0
    zf = getztraj(model, state0, 0.)
    ω = getomega(model.g, z0 - dzdx0 * x0)

    xrange = 1.2 * [min(x0, 0.); max(x0, 0.)]
    xdrange = sort(-ω * xrange)

    numstillframes = stilltime * fps / realtimerate
    numsimframes = simtime * fps / realtimerate
    ts = linspace(0, simtime, numsimframes)

    # title
    # realtimerate_text = (@sprintf "%0.2f" realtimerate) * "x"
    # tit = suptitle(realtimerate_text)

    ts, states = ode45((t, state) -> odefun(model, state), state0, ts, points=:specified)
    xs = [getx(model, state)::Float64 for state in states]
    zs = [getz(model, state)::Float64 for state in states]
    xds = [getxd(model, state)::Float64 for state in states]
    zds = [getzd(model, state)::Float64 for state in states]
    us = [leg_u(model, state) for state in states]
    leg_forces = [u * [getx(model, state); getz(model, state)] for (u, state) in zip(us, states)]
    normalized_leg_force_intensities = [dot(f, [x; z] / (model.g * norm([x; z]))) for (x, z, f) in zip(xs, zs, leg_forces)]

    # state space
    if ss_ax != nothing
        let
            plt[:sca](ss_ax)
            xlim(xrange)
            ylim(xdrange)
            xticks([-0.3; -0.2; -0.1; 0.])

            # axis lines
            axhline(linewidth = 0.5, color = "black")
            axvline(linewidth = 0.5, color = "black")

            # eigenvectors
            # plot(xrange,  ω * xrange, "k--", linewidth = 1.0, zorder = 1)
            icpline = plot(xrange, -ω * xrange, "k--", linewidth = 1.0, zorder = 1, label = L"x + \sqrt{\frac{z_\mathrm{f}}{g}} \dot{x} = 0")
            legend(handles = collect(icpline), frameon = true, handlelength = 3, handletextpad = 0, borderaxespad = 0.5, fontsize = 14, loc = "lower left")

            # ICP line annotation
            # xarrowtip = xrange[2] - 0.1 * diff(xrange)
            # arrowtip = [xarrowtip; -ω * xarrowtip]
            # offset = [-diff(xrange) / 2; 0]
            # annotate(L"x + \sqrt{\frac{z}{g}} \dot{x} = 0", xy = arrowtip, xytext = arrowtip + offset, xycoords="data", arrowprops = Dict("facecolor"=>"black"))

            # axis labels
            xlabel(L"x" * " [m]")
            ylabel(L"\dot x" * " [m/s]")
        end

        # dynamic objects
        plt[:sca](ss_ax)
        traj = plot([], [], color = "red")[1]
        point = plot([], [], "ro")[1]
        validcontour = nothing
    end

    # force
    if force_ax != nothing
        let
            plt[:sca](force_ax)
            xlim(xrange)
            ylim(min(0, normalized_leg_force_intensities...), 1.1 * max(normalized_leg_force_intensities...))
            xticks([-0.3; -0.2; -0.1; 0.])

            # axis line
            axhline(linewidth = 0.5, color = "black")

            # axis labels
            xlabel(L"x" * " [m]")
            # ylabel(L"\mathbf{f}_{\mathrm{leg}} \cdot \frac{\mathbf{q}}{|| \mathbf{q} ||}")
            # ylabeltext = staticplot ? L"\frac{\mathbf{f}_{\mathrm{gr}}}{mg} \cdot \frac{\mathbf{q}}{|| \mathbf{q} ||}" * " [-]" : "normalized leg force [-]"
            ylabel("normalized leg force [-]")
        end

        # dynamic objects
        plt[:sca](force_ax)
        utraj = plot([], [], color = force_color)[1]
        upoint = plot([], [], "o", color = force_color)[1]
    end

    # model
    if model_ax != nothing
        com_radius = 0.05
        xrange_model = [xrange[1] - 1.5 * com_radius; xrange[2] + 1.5 * com_radius]
        zrange_model = [-0.1 * z0; 1.3 * z0]
        # force_scale = (diff(zrange_model)[1] / 4) / max(map(f -> f[2], leg_forces)...)
        force_scale = (diff(zrange_model)[1] / 4) / model.g
        velocity_scale = (diff(xrange_model)[1] / 3)

        plt[:sca](model_ax)

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
        # plt[:sca](model_ax)
        leg = plot([], [], color = "black")[1]
        com = CoMSymbol(model_ax, com_radius, 2.; zorder = 4)
        function add_arrow_patch(ax, color)
            arrow = patches.FancyArrowPatch([0.; 0.], [0.; 0.], arrowstyle = "simple", fc = color, ec = color, mutation_scale = 20, zorder = 5)
            ax[:add_patch](arrow)
        end

        if !staticplot
            fgrav = add_arrow_patch(model_ax, force_color)
            fgravtext = text(0., 0., L"m \mathbf{g}", fontsize = font_size)
            fgravtext[:set_clip_on](true)

            fleg = add_arrow_patch(model_ax, force_color)
            flegtext = text(0., 0., L"\mathbf{f}_{\mathrm{gr}}", fontsize = font_size)
            flegtext[:set_clip_on](true)
        end

        qdot = add_arrow_patch(model_ax, velocity_color)
        qdottext = text(0., 0., L"\dot{\mathbf{q}}_0", fontsize = font_size)
        qdottext[:set_clip_on](true)

        ztraj_sim = plot([], [], color = ztraj_color, zorder = 6, label = "simulation")[1]

        # legend(handles = [ztraj; ztraj_sim], frameon = true, handlelength = 3, handletextpad = 0, borderaxespad = 0.5, fontsize = 12)
    end

    plt[:tight_layout]()
    plt[:subplots_adjust](top = 0.92)

    if staticplot
        trange = length(ts) : length(ts)
    else
        FFMpegWriter = get(anim.writers, "ffmpeg")
        writer = FFMpegWriter(fps = fps, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
        writer[:setup](fig, filename, dpi)
        trange = 1 : length(ts)
    end

    for i in trange
        state_ind = staticplot ? 1 : i
        x = xs[state_ind]
        z = zs[state_ind]
        xd = xds[state_ind]
        zd = zds[state_ind]

        if show_orbital_energy
            Eo = orbital_energy(model, [x; xd])
            tit[:set_text](realtimerate_text * ", Orbital energy: " * (@sprintf "%0.2f" Eo)) #* ", " * L"z = " * "$z, " * L"\dot{z} = " * "$zd")
        end

        # state space plot
        if ss_ax != nothing
            traj[:set_data](xs[1 : i], xds[1 : i])
            point[:set_data](xs[state_ind], xds[state_ind])

            # region of validity
            if show_region
                plt[:sca](ss_ax)
                if validcontour != nothing
                    for collection in validcontour[:collections]
                        collection[:remove]()
                    end
                end
                xlinspace = linspace(xrange..., 100)
                xdlinspace = linspace(xdrange..., 100)
                x_grid = [i::Float64 for i in xlinspace, j in xdlinspace]
                xd_grid = [j::Float64 for i in xlinspace, j in xdlinspace]
                opposite_velocity_mask = max(-sign((x_grid .* xd_grid)), 0)
                valid = is_force_always_nonnegative_condition2(model.g, x_grid, z, xd_grid, zd, zf) .* opposite_velocity_mask
                validcontour = contourf(x_grid, xd_grid, valid, levels = [0, Inf], colors = "g", alpha = 0.6, zorder = 1)
            end
        end

        # force plot
        if force_ax != nothing
            utraj[:set_data](xs[1 : i], normalized_leg_force_intensities[1 : i])
            upoint[:set_data](xs[state_ind], normalized_leg_force_intensities[state_ind])
        end

        # model plot
        if model_ax != nothing
            leg[:set_data]([0; x], [0; z])
            set_center(com, [x; z])

            if !staticplot
                fgrav_text_offset = [-0.1; 0.]
                fleg_text_offset = [0.03; 0.]

                fgrav_tip = [x; z] + force_scale * [0; -model.g]
                fgrav[:set_positions]([x; z], fgrav_tip)
                fgravtext[:set_position](fgrav_tip + fgrav_text_offset)

                fleg_tip = force_scale * leg_forces[state_ind]
                fleg[:set_positions]([0; 0], fleg_tip)
                flegtext[:set_position](fleg_tip + fleg_text_offset)
            end

            if staticplot || (i == 1 && numstillframes > 0)
                qdot_text_offset = [0.; 0.03]
                qdot_tip = [x; z] + velocity_scale * [xd; zd]
                qdot[:set_positions]([x; z], qdot_tip)
                qdottext[:set_position](qdot_tip + qdot_text_offset)
            else
                qdot[:set_visible](false)
                qdottext[:set_visible](false)
            end

            ztraj_sim[:set_data](xs[1 : i], zs[1 : i])
        end

        if !staticplot
            numframes = i == 1 ? numstillframes : 1
            for j in 1 : numframes
                writer[:grab_frame]()
            end
        end
    end

    if staticplot
        # ion()
        # figure()
        # show()
    else
        writer[:saving]()
        writer[:finish]()
    end

    plt[:rcdefaults]()

    fig
end
