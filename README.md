# SNSensitivityEstimate

[![Build Status](https://github.com/shoram444/SNSensitivityEstimate.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/shoram444/SNSensitivityEstimate.jl/actions/workflows/CI.yml?query=branch%3Amain)

# SNSensitivityEstimate.jl
Julia package for estimating sensitivity in SuperNEMO

# Installation
### Install Julia v1.10 with juliaup

1. Install juliaup
```bash
curl -fsSL https://install.julialang.org | sh
```
On Windows, install it via the Microsoft Store: [Juliaup on Microsoft Store](https://aka.ms/julia-ms-store)

2. Restart your shell (or run source ~/.bashrc, source ~/.zshrc, etc.) This ensures julia is available on your PATH.

3. Install Julia 1.10:

```bash
juliaup add 1.10
```

4. Set Julia 1.10 as the default version:

```bash
juliaup default 1.10
```

5. Launch Julia:
```bash
julia
```

# Install SNSensitivityEstimate.jl 
1. Open up julia
```bash
julia
```
2. Open package manager and install `SNSensitivityEstimate` from github: 
```bash
import Pkg
Pkg.add("https://github.com/Shoram444/SNSensitivityEstimate.jl")
```
3. Use `SNSensitivityEstimate`
```bash
using SNSensitivityEstimate
```

# Examples
There's a small example in the directory `examples` which show how one can calculate the sensitivity to a mock process with gaussian signal and exponential background. The example relies on the package `BlackBoxOptim` for finding the NDimensional ROI to calcualte the sensitivity for.
