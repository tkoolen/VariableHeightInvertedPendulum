# from https://github.com/JuliaLang/IJulia.jl/issues/107

function create_mp4(fig, images, filename; fps = 30)
    @pyimport matplotlib.animation as anim
    ani = anim.ArtistAnimation(fig, images, interval=1000/fps, blit=true, repeat_delay=1000)
    ani[:save](filename, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
end

function display_mp4(filename)
    display("text/html", string("""<video autoplay controls><source src="data:video/x-m4v;base64,""",
    base64encode(open(readbytes,filename)),"""" type="video/mp4"></video>"""))
end
