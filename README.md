# ExtremalOptimization

[![Build Status](https://github.com/francescoalemanno/ExtremalOptimization.jl/workflows/CI/badge.svg)](https://github.com/francescoalemanno/ExtremalOptimization.jl/actions)

## API:

```julia
optimize(f, s, N; β=1.0, γ=1.2, rng=Random.GLOBAL_RNG, reps_per_particle=100, callback=state->nothing)
```

- `f` : cost function to minimize, whose argument is either a scalar or a vector, must returns a scalar value.
- `s` : function whose input is the particle number and output is a random initial point to be ranked by `f`.
- `N` : number of particles to use, choose a number greater than `d+4` where `d` is the number of dimensions.
- `reps_per_particle` : total number of iterations per particle.

## Usage example:

```julia
rosenbrock2d(x) = (x[1]-1)^2+(x[2]-x[1]^2)^2
initpoint(i) = randn(2)
optimize(rosenbrock2d, initpoint, 50)
```
output
```
(x = [0.9999999999993153, 0.9999999999997592], fx = 1.742385603613753e-24)
```
as expected the algorithm has found the optimum at `(1, 1)`.