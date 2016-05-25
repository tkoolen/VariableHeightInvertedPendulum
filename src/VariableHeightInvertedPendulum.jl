module VariableHeightInvertedPendulum

using PyCall
using PyPlot
using ODE

@pyimport matplotlib.patches as patches
@pyimport matplotlib.animation as anim

include("plotutil.jl")
getomega(g, z) = √(g / z)
include("lipm.jl")
include("com_symbol.jl")


function sim_movie(model, state0, filename::ASCIIString; fps = 30., realtimerate = 0.25, dpi = 100, simtime = 1., test = false)
    ztraj_color = "#364b9f" # blue
    force_color = "#f37620" # orange

    # test ? ion() : ioff()
    ioff()
    fig = figure()

    x0 = getx(model, state0)
    xd0 = getxd(model, state0)
    z0 = getz(model, state0)

    xrange = [-abs(x0); abs(x0)]
    xdrange = [-abs(xd0); abs(xd0)]

    ts = linspace(0, simtime, 100)
    numframes = simtime * fps / realtimerate
    ω = getomega(model.g, z0)

    E = orbital_energy(model, state0)

    ts, states = ode45((t, state) -> odefun(model, state), state0, ts, points=:specified)
    xs = [getx(model, state)::Float64 for state in states]
    zs = [getz(model, state)::Float64 for state in states]
    xds = [getxd(model, state)::Float64 for state in states]
    zds = [getzd(model, state)::Float64 for state in states]
    leg_forces = [normalized_leg_force(model, state) for state in states]

    let
        ax = subplot(1, 2, 1)
        xlim(xrange)
        ylim(xdrange)

        # axis lines
        axhline(linewidth = 1, color = "black")
        axvline(linewidth = 1, color = "black")

        # eigenvectors
        plot(xrange,  ω * xrange, color = "black", linewidth = 2.0)
        plot(xrange, -ω * xrange, color = "black", linewidth = 2.0)

        # ICP line annotation
        xarrowtip = xrange[2] - 0.1 * diff(xrange)
        arrowtip = [xarrowtip; -ω * xarrowtip]
        offset = [-diff(xrange) / 2; 0]
        annotate(L"x + \sqrt{\frac{z}{g}} \dot{x} = 0", xy = arrowtip, xytext = arrowtip + offset, xycoords="data", arrowprops = Dict("facecolor"=>"black"))

        # axis labels
        xlabel(L"x" * " (m)")
        ylabel(L"\dot x" * " (m/s)")
    end

    # model
    com_radius = 0.05
    xrange_model = [xrange[1] - 1.5 * com_radius; xrange[2] + 1.5 * com_radius]
    zrange_model = [-0.1 * z0; 1.5 * z0]
    force_scale = (diff(zrange_model)[1] / 4) / max(map(f -> f[2], leg_forces)...)

    let
        ax = subplot(1, 2, 2)
        xlim(xrange_model)
        ylim(zrange_model)
        xlabel(L"x" * " (m)")
        ylabel(L"z" * " (m)")
        plt[:axis]("equal")

        # z trajectory
        xs_linspace = linspace(xrange_model[1], xrange_model[2], 100)
        ztraj = [getztraj(model, x)::Float64 for x in xs_linspace]
        plot(xs_linspace, ztraj, color = ztraj_color)

        # foot
        plot([0], [0], "ko", zorder = 4)
    end

    # state space plot dynamic objects
    subplot(1, 2, 1)
    traj = plot([], [], color = "red")[1]
    point = plot([], [], "ro")[1]

    # model plot dynamic objects
    ax = subplot(1, 2, 2)
    leg = plot([], [], color = "black")[1]
    com = CoMSymbol(ax, com_radius, 2.; zorder = 4)
    function add_force()
        arrow = patches.FancyArrowPatch([0.; 0.], [0.; 0.], arrowstyle = "simple", fc = force_color, ec = force_color, mutation_scale = 20, zorder = 3)
        ax[:add_patch](arrow)
    end
    fgrav = add_force()
    fgravtext = text(0., 0., L"m \mathbf{g}")
    fleg = add_force()
    flegtext = text(0., 0., L"\mathbf{f}_{\mathrm{leg}}")

    if !test
        FFMpegWriter = get(anim.writers, "ffmpeg")
        writer = FFMpegWriter(fps = fps, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
        writer[:setup](fig, filename, dpi)
    else
        ts = ts[end]
    end

    for i in 1 : length(ts)
        x = xs[i]
        z = zs[i]
        xd = xds[i]
        zd = zds[i]

        # title
        suptitle("Orbital energy: " * (@sprintf "%0.2f" E) * ", " * L"z = " * "$z, " * L"\dot{z} = " * "$zd")

        # state space plot
        traj[:set_data](xs[1 : i], xds[1 : i])
        point[:set_data](xs[i], xds[i])

        # model plot
        leg[:set_data]([0; x], [0; z])
        set_center(com, [x; z])

        force_text_offset = [0.03; 0.]

        fgrav_tip = [x; z] + force_scale * [0; -model.g]
        fgrav[:set_positions]([x; z], fgrav_tip)
        fgravtext[:set_position](fgrav_tip + force_text_offset)

        fleg_tip = force_scale * leg_forces[i]
        fleg[:set_positions]([0; 0], fleg_tip)
        flegtext[:set_position](fleg_tip + force_text_offset)

        if !test
            writer[:grab_frame]()
        end
    end

    if test
        ion()
        figure()
        show()
    else
        writer[:saving]()
        writer[:finish]()
    end
end

export
    LIPM,
    sim_movie,
    # create_mp4,
    display_mp4
end
