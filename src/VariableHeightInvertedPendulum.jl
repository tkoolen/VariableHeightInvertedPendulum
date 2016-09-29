module VariableHeightInvertedPendulum

using PyCall
using PyPlot
using ODE
using Polynomials
using ForwardDiff

@pyimport matplotlib.patches as patches
@pyimport matplotlib.animation as anim

# TODO: move
const ztraj_color = "#364b9f" # blue
const force_color = "#f37620" # orange
const velocity_color = "#be2025" # red

getomega(g, z) = √(g / z)

include("scenario.jl")
include("plotutil.jl")
include("lipm.jl")
include("polynomial_variable_height.jl")
include("cubic_controlled_variable_height_inverted_pendulum.jl")
include("com_symbol.jl")
include("simulation_results.jl")
include("state_space_axes.jl")
include("force_axes.jl")
include("model_axes.jl")
include("simulation_figure.jl")
include("static_sim_figure.jl")
include("sim_movie.jl")

export
    LIPM,
    PolynomialVariableHeightModel,
    CubicControlledVariableHeightInvertedPendulum,
    Scenario,
    velocity_given_orbital_energy,
    sim_movie,
    sim_figure,
    zero_orbital_energy_trajectory,
    is_force_always_nonnegative
end
