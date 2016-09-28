module VariableHeightInvertedPendulum

using PyCall
using PyPlot
using ODE
using Polynomials
using ForwardDiff

@pyimport matplotlib.patches as patches
@pyimport matplotlib.animation as anim

getomega(g, z) = √(g / z)

include("scenario.jl")
include("plotutil.jl")
include("lipm.jl")
include("polynomial_variable_height.jl")
include("cubic_controlled_variable_height_inverted_pendulum.jl")
include("com_symbol.jl")
include("sim_movie.jl")

export
    LIPM,
    PolynomialVariableHeightModel,
    CubicControlledVariableHeightInvertedPendulum,
    Scenario,
    velocity_given_orbital_energy,
    sim_movie,
    zero_orbital_energy_trajectory,
    is_force_always_nonnegative
    # create_mp4,
    # display_mp4
end
