function sim_figure(scenario::Scenario;
    simtime = 1., restrict_ztraj = false, show_region = false, show_icp_line = true, model_only = false, show_orbital_energy = false, font_size = 15, fig_size = (8., 6.),
    show_gravity = false, show_leg_force = false, export_dir = "../figures", do_export = false)

    model = scenario.model
    state0 = scenario.initial_conditions

    ts = 0 : 0.005 : simtime
    results = SimulationResults(model, state0, ts)

    ioff()
    fig = SimulationFigure(results, fig_size, font_size, model_only, restrict_ztraj, show_region, show_icp_line, show_orbital_energy, show_gravity, show_leg_force)

    update(fig, results, 1 : length(ts), 1)

    if do_export
        !isdir(export_dir) && mkdir(export_dir)
        filename = "$export_dir/$(scenario.model.name)_$(scenario.name).pdf"
        savefig(filename)
    end

    plt[:rcdefaults]()
    nothing
end
