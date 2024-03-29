A mishmash of Julia implementations for uncertainty quantification in geophysical and climate models, as I've been using them in my postgraduate research.

## Why?
This is for my own purposes - I reuse a lot of code across various projects, and so including it all in a single, tested Julia package is convenient.
This is also an oppourtunity to work on maintaining something like a software package, with documentation and CI, which are skills I am not working on elsewhere.
I don't expect this to be useful to anyone else, but I am keeping it on Github for the sheer sake of it.

## What is here?
Within the main module `MPhil`, everything is loosely organised into submodules:

- `StochasticSensitivity`: Computation of the stochastic sensitivity tools, first introduced by [Balasuriya (2020)](https://epubs.siam.org/doi/10.1137/18M1222922) and extended by [Blake, Maclean and Balasuriya (2023)](), for uncertainty quantification in differential equation models.

- `Numerics`: Custom implementations of finite-difference approximations and numerical solvers of integrals and differential equations.

- `OceanData`: Manipiulation and processing of oceanography data, e.g. geostrophic velocity data.

- `Visualisation`: Tools for various visualisation of results with [`Makie.jl`](https://docs.makie.org/stable/).

See within each subdirectory of `src` for more details.

