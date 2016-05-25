module VariableHeightInvertedPendulum

using PyCall
using PyPlot
using ODE
using Polynomials

@pyimport matplotlib.patches as patches
@pyimport matplotlib.animation as anim

getomega(g, z) = √(g / z)

include("scenario.jl")
include("plotutil.jl")
include("lipm.jl")
include("polynomial_variable_height.jl")
include("com_symbol.jl")
include("sim_movie.jl")

export
    LIPM,
    PolynomialVariableHeightModel,
    Scenario,
    velocity_given_orbital_energy,
    sim_movie,
    # create_mp4,
    display_mp4
end
