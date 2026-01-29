# This folder contains 3 examples for how to work with the `SNSensitivityEstimate.jl` package

## The prerequisites to using this package are:
- Have julia installed (versions 1.10+ should be fine) or use `module load julia` if working on `CC` 
- To install the required dependencies:
  - `cd PATH/TO/EXAMPLES/` (must be in the path where `Project.toml` and `Manifest.toml` files are located)
  - Open julia REPL by typing `julia` in bash
  - When in REPL type `import Pkg; Pkg.activate("."); Pkg.instantiate()` which will activate the local environment and automatically load the required dependencies 
  - To run any particular script `cd` to where the files are located, i.e. `cd examples` and run in bash `julia --project="PATH/TO/Project.toml" ex1D.jl` 
  - That's it
  - (RECOMMENDED WORKFLOW, but not needed) While it is possible to use julia in "scripting fashion" (by running .jl scripts in bash) it is not the best approach. You should really consider using the workflow of VSCode + julia extension and work in REPL itself with one session. 

# Use of the examples

## ex1D.jl
This is a very simple example of calculating the 1D sensitivity ROI for a fake signal (gaussian) and fake backrgound (exponential) 
It uses only the sumE channel to find the simple 1D roi and calculates the sensitivity. 
To use the 