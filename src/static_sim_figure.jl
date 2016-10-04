function sim_figure(scenario::Scenario;
    simtime = 1., restrict_ztraj = false, show_region = false, show_icp_line = true, model_only = false, show_orbital_energy = false, font_size = 15, fig_size = (8., 5.5),
    show_gravity = false, show_leg_force = false, export_dir = "../figures", do_export = false, snapshots = 1)

    model = scenario.model
    state0 = scenario.initial_conditions

    ts = 0 : 0.005 : simtime
    results = SimulationResults(model, state0, ts)

    ioff()
    fig = SimulationFigure(results, fig_size, font_size, model_only, restrict_ztraj, show_region, show_icp_line, show_orbital_energy, show_gravity, show_leg_force)

    for i in 0 : snapshots - 1
      time_index = Int64(floor((length(ts) - 1) * i / (snapshots - 1))) + 1
      if i > 0
        fig.model_ax.show_qdot = false
      end
      update(fig, results, 1 : length(ts), time_index)

      if do_export
          !isdir(export_dir) && mkdir(export_dir)
          index_part = snapshots > 1 ? "_$i" : ""
          filename = "$export_dir/$(scenario.model.name)_$(scenario.name)$index_part.pdf"
          savefig(filename)
      end
    end

    plt[:rcdefaults]()
    nothing
end
