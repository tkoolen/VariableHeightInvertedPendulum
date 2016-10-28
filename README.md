# VariableHeightInvertedPendulum
[![Build Status](https://travis-ci.org/tkoolen/VariableHeightInvertedPendulum.svg?branch=master)](https://travis-ci.org/tkoolen/VariableHeightInvertedPendulum)

This repository contains code associated with the paper "Balance control using center of mass height variation: limitations imposed by unilateral contact".

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
The Mathematica notebook in the `mathematica` directory contains checks/derivations of all of the main results in the paper (including the the application of the CAD algorithm). It also contains code for generating two of the plots in the paper.

The Mathematica code has only been tested in Mathematica 10.4.1.0.

## Julia
The `src` directory contains Julia code for running simulations of the variable-height inverted pendulum and creating the simulation plots presented in the paper. Animation of the simulation results is also available. The `notebook` directory contains IJulia notebooks that visualize simulation results for a few initial conditions. 

Non-interactive online preview of the notebooks:
* [figures from paper](http://nbviewer.jupyter.org/github/tkoolen/VariableHeightInvertedPendulum/blob/master/notebook/PaperFigures.ipynb)
* [additional figures and animations](http://nbviewer.jupyter.org/github/tkoolen/VariableHeightInvertedPendulum/blob/master/notebook/AdditionalFiguresAndAnimations.ipynb)

The notebooks can be easily modified to visualize additional initial conditions. Interacting with the notebooks requires installing Julia as well as this package. The Julia code has been tested on Ubuntu 14.04 and 16.04, and on OSX 10.11.6. To use the Julia code, perform the following steps:

1. Install Julia 0.5 from http://julialang.org/downloads/.
1. (Optional, but recommended) create a separate Julia package directory (which will contain the code in this repo and all dependencies). If you skip this step, everything will get installed into the default, global julia package directory, which is alright if you're a casual user but can potentially lead to dependency version clashes if you install other packages. On OSX/Linux:
``mkdir ~/VariableHeightInvertedPendulum && cd ~/VariableHeightInvertedPendulum && export JULIA_PKGDIR=`pwd` `` 
1. Run `julia` from the command line.
1. Enter `Pkg.init()` to initialize the package directory.
1. Enter `Pkg.clone("https://github.com/tkoolen/VariableHeightInvertedPendulum.git")` to get this package and all of its dependencies.

After this installation procedure, the IJulia notebooks can be run using the following commands from the Julia command line:

1. `using VariableHeightInvertedPendulum`
1. `IJulia.notebook(dir = Pkg.dir("VariableHeightInvertedPendulum") * "/notebook")`

A browser window should pop up, showing various notebooks. Click one to open; the rest should speak for itself (shift + Enter is the shortcut to run a cell and advance to the next).

If you didn't skip the optional installation step, make sure you set `JULIA_PKGDIR` before starting Julia (``cd ~/VariableHeightInvertedPendulum && export JULIA_PKGDIR=`pwd` ``) any time you open a new terminal.

To create animations, you may have to install `ffmpeg` (`sudo apt-get install ffmpeg` on Ubuntu 16.04, `brew install ffmpeg` on OSX). For Ubuntu 14.04 see the `.travis.yml` for how to set up ffmpeg.

Feel free to open an issue if there are any problems.

