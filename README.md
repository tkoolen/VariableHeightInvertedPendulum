# VariableHeightInvertedPendulum

This repository contains code associated with

```tex
@inproceedings{Koolen2016balance,
  title={Balance control using center of mass height variation: limitations imposed by unilateral contact},
  author={Koolen, Twan and Posa, Michael and Tedrake, Russ}
  booktitle={Humanoid Robots (Humanoids), 2016 IEEE-RAS 16th International Conference on, accepted},
  year={2016},
  month={Nov}
}

```

## Mathematica
The Mathematica notebook in the `mathematica` directory contains checks/derivations of all of the main results in the paper, as well as the application of the CAD algorithm. It also contains code for generating two of the plots in the paper.

The Mathematica code has only been tested in Mathematica 10.4.1.0.

## Julia
The `src` and `notebook` directories contain Julia code used to run the simulations and create related plots presented in the paper. The `notebook` directory contains IJulia notebooks that can be easily modified to visualize additional initial conditions. Visualizing these simulations using animations is also possible.

To install the Julia code, perform the following steps:

1. Install Julia 0.5 from http://julialang.org/downloads/.
1. (Optional, but recommended) create a separate Julia package directory (which will contain the code in this repo and all dependencies). If you skip this step, everything will get installed into the default, global julia package directory, which may be alright. On OSX/Linux:
``mkdir ~/VariableHeightInvertedPendulum && cd ~/VariableHeightInvertedPendulum && export JULIA_PKGDIR=`pwd``` 
1. Run `julia` from the command line.
1. Enter `Pkg.init()` to initialize the package directory.
1. Enter `Pkg.clone("https://github.com/tkoolen/VariableHeightInvertedPendulum.git")` to get this package and all of its dependencies.

After this installation procedure, the IJulia notebooks can be run using the following commands:

1. `julia using VariableHeightInvertedPendulum`
1. `IJulia.notebook(dir = Pkg.dir("VariableHeightInvertedPendulum") * "/notebook")`

A browser window should pop up, showing various notebooks. Click one to open; the rest should speak for itself.

If you're using a new terminal and didn't skip the optional installation step, make sure you set `JULIA_PKGDIR` before starting Julia (``cd ~/VariableHeightInvertedPendulum && export JULIA_PKGDIR=`pwd```)

To create animations, you may have to install `ffmpeg` (`sudo apt-get install ffmpeg` on Ubuntu, `brew install ffmpeg` on OSX).

