type SimulationFigure
    fig
    axes
    model_ax
    show_orbital_energy

    function SimulationFigure(results::SimulationResults, fig_size, font_size, model_only, restrict_ztraj, show_region, show_icp_line, show_orbital_energy, show_gravity, show_leg_force)
        fig = figure(figsize = fig_size)
        plt[:rc]("font", size = font_size, family = "serif")

        # don't use type 3 fonts
        # http://phyletica.org/matplotlib-fonts/
        plt[:rc]("pdf", fonttype = 42)
        plt[:rc]("ps", fonttype = 42)

        axes = []

        model = results.model
        state0 = results.state0

        x0 = getx(model, state0)
        xd0 = getxd(model, state0)
        z0 = getz(model, state0)
        zd0 = getzd(model, state0)
        dzdx0 = zd0 / xd0
        ω = getomega(model.g, z0 - dzdx0 * x0)

        x0 = getx(model, state0)
        xrange = 1.2 * [min(x0, 0.); max(x0, 0.)]
        xdrange = sort(-ω * xrange)

        axes = []
        if model_only
            model_ax = ModelAxes(gca(), results, xrange, font_size, restrict_ztraj, show_gravity, show_leg_force)
        else
            push!(axes, StateSpaceAxes(plt[:subplot2grid]((2,2), (0,0)), xrange, xdrange, ω, show_region, show_icp_line))
            push!(axes, ForceAxes(plt[:subplot2grid]((2,2), (1,0)), xrange, results.normalized_leg_force_intensities))
            model_ax = ModelAxes(plt[:subplot2grid]((2,2), (0, 1), rowspan=2), results, xrange, font_size, restrict_ztraj, show_gravity, show_leg_force)
        end
        push!(axes, model_ax)

        plt[:tight_layout]()
        plt[:subplots_adjust](top = 0.92)

        new(fig, axes, model_ax, show_orbital_energy)
    end
end

function update(sim_fig::SimulationFigure, results::SimulationResults, state_ind_range::UnitRange{Int64}, current_state_ind::Int64)
    if sim_fig.show_orbital_energy
        x = results.xs[state_ind]
        xd = results.xds[state_ind]
        Eo = orbital_energy(model, [x; xd])
        tit[:set_text](realtimerate_text * ", Orbital energy: " * (@sprintf "%0.2f" Eo)) #* ", " * L"z = " * "$z, " * L"\dot{z} = " * "$zd")
    end

    for axes in sim_fig.axes
        update(axes, results, state_ind_range, current_state_ind)
    end

    nothing
end
