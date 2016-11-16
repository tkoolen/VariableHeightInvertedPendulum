using IJulia
# import IJulia: InlineDisplay

function sim_movie(scenario::Scenario;
    fps = 30., realtimerate = 0.25, dpi = 100, simtime = 1., stilltime = 0.5,
    restrict_ztraj = false, show_region = false, show_region2 = false, show_icp_line = false, model_only = false,
    show_orbital_energy = false, font_size = 15, fig_size = (8., 6.), show_gravity = true, show_leg_force = true, show_ballistic = false,
    export_dir = "../movies", do_export = false)

    model = scenario.model
    state0 = scenario.initial_conditions

    !isdir(export_dir) && mkdir(export_dir)
    filename = "$export_dir/$(scenario.model.name)_$(scenario.name).mp4"

    numstillframes = stilltime * fps / realtimerate
    numsimframes = simtime * fps / realtimerate
    ts = linspace(0, simtime, numsimframes)

    results = SimulationResults(model, state0, ts)

    ioff()
    fig = SimulationFigure(results, fig_size, font_size, model_only, restrict_ztraj,
        show_region, show_region2, show_icp_line, show_orbital_energy, show_gravity, show_leg_force, show_ballistic)

    FFMpegWriter = get(anim.writers, "ffmpeg")
    writer = FFMpegWriter(fps = fps, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])

    writer[:setup](fig.fig, filename, dpi)

    for i in 1 : length(ts)
        fig.model_ax.show_qdot = i == 1
        update(fig, results, 1 : i, i)

        numframes = i == 1 ? numstillframes : 1
        for j in 1 : numframes
            writer[:grab_frame]()
        end
    end

    if do_export
        writer[:saving]()
    end
    writer[:finish]()

    plt[:close](fig.fig)
    plt[:rcdefaults]()

    if do_export && any(map(d -> isa(d, IJulia.InlineDisplay), Base.Multimedia.displays))
        display_mp4(filename)
    end

    nothing
end
