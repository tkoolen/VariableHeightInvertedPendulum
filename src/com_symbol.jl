type CoMSymbol
    circle
    wedges

    function CoMSymbol(ax, radius, linewidth; zorder = 1)
        circle = patches.Circle([0.; 0.], radius, color = "black", fill = false, linewidth = linewidth, zorder = zorder)
        ax[:add_artist](circle)
        wedges = []
        for i = 1 : 4
            color = isodd(i) ? "black" : "white"
            wedge = patches.Wedge([0.; 0.], radius, (i - 1) * 90, i * 90, fc = color, zorder = zorder)
            ax[:add_artist](wedge)
            push!(wedges, wedge)
        end
        new(circle, wedges)
    end
end

function set_center(com::CoMSymbol, center::Vector)
    com.circle[:center] = center
    for wedge in com.wedges
        wedge[:set_center](center)
    end
end
