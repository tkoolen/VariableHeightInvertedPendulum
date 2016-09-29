type StateSpaceAxes
    ax
    xrange::Vector{Float64}
    xdrange::Vector{Float64}
    show_region
    traj
    point
    validcontour

    function StateSpaceAxes(ax, xrange, xdrange, ω, show_region, show_icp_line)
        # static objects
        let
            plt[:sca](ax)
            xlim(xrange)
            ylim(xdrange)
            xticks([-0.3; -0.2; -0.1; 0.])

            # axis lines
            axhline(linewidth = 0.5, color = "black")
            axvline(linewidth = 0.5, color = "black")

            # eigenvectors
            # plot(xrange,  ω * xrange, "k--", linewidth = 1.0, zorder = 1)
            if show_icp_line
                icpline = plot(xrange, -ω * xrange, "k--", linewidth = 1.0, zorder = 1, label = L"x + \sqrt{\frac{z_\mathrm{f}}{g}} \dot{x} = 0")
                legend(handles = collect(icpline), frameon = true, handlelength = 3, handletextpad = 0, borderaxespad = 0.5, fontsize = 14, loc = "lower left")
            end

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
        plt[:sca](ax)
        traj = plot([], [], color = "red")[1]
        point = plot([], [], "ro")[1]
        validcontour = nothing

        new(ax, xrange, xdrange, show_region, traj, point, validcontour)
    end
end

function update(ss_ax::StateSpaceAxes, results::SimulationResults, state_ind_range::UnitRange{Int64}, current_state_ind::Int64)
    ss_ax.traj[:set_data](results.xs[state_ind_range], results.xds[state_ind_range])
    ss_ax.point[:set_data](results.xs[current_state_ind], results.xds[current_state_ind])

    # region of validity
    if ss_ax.show_region
        model = results.model
        z = results.zs[current_state_ind]
        zd = results.zds[current_state_ind]
        zf = getzf(model)

        plt[:sca](ss_ax.ax)
        if ss_ax.validcontour != nothing
            for collection in ss_ax.validcontour[:collections]
                collection[:remove]()
            end
        end
        xlinspace = linspace(ss_ax.xrange..., 100)
        xdlinspace = linspace(ss_ax.xdrange..., 100)
        x_grid = [i::Float64 for i in xlinspace, j in xdlinspace]
        xd_grid = [j::Float64 for i in xlinspace, j in xdlinspace]
        opposite_velocity_mask = max(-sign((x_grid .* xd_grid)), 0)
        valid = is_force_always_nonnegative_condition2(model.g, x_grid, z, xd_grid, zd, zf) .* opposite_velocity_mask
        ss_ax.validcontour = contourf(x_grid, xd_grid, valid, levels = [0, Inf], colors = "g", alpha = 0.6, zorder = 1)
    end
    nothing
end
