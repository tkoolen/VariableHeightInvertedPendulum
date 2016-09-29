type ForceAxes
    ax
    utraj
    upoint

    function ForceAxes(ax, xrange, normalized_leg_force_intensities)
        # static objects
        let
            plt[:sca](ax)
            xlim(xrange)
            ylim(min(0, normalized_leg_force_intensities...), 1.1 * max(normalized_leg_force_intensities...))
            xticks(collect(tickmark_locations(xrange, 0.1)))

            # axis line
            axhline(linewidth = 0.5, color = "black")

            # axis labels
            xlabel(L"x" * " [m]")
            # ylabel(L"\mathbf{f}_{\mathrm{leg}} \cdot \frac{\mathbf{q}}{|| \mathbf{q} ||}")
            # ylabeltext = staticplot ? L"\frac{\mathbf{f}_{\mathrm{gr}}}{mg} \cdot \frac{\mathbf{q}}{|| \mathbf{q} ||}" * " [-]" : "normalized leg force [-]"
            ylabel("normalized leg force [-]")
        end

        # dynamic objects
        plt[:sca](ax)
        utraj = plot([], [], color = force_color)[1]
        upoint = plot([], [], "o", color = force_color)[1]

        new(ax, utraj, upoint)
    end
end

function update(force_ax::ForceAxes, results::SimulationResults, state_ind_range::UnitRange{Int64}, current_state_ind::Int64)
    force_ax.utraj[:set_data](results.xs[state_ind_range], results.normalized_leg_force_intensities[state_ind_range])
    force_ax.upoint[:set_data](results.xs[current_state_ind], results.normalized_leg_force_intensities[current_state_ind])
    nothing
end
