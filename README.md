# ExtremalOptimization

for functions of continuous variables.

[![Build Status](https://github.com/francescoalemanno/ExtremalOptimization.jl/workflows/CI/badge.svg)](https://github.com/francescoalemanno/ExtremalOptimization.jl/actions)

This package implements the basic mechanism of Extremal Optimization (τ-EO) as described in [Boettcher, Stefan; Percus, Allon G. (2001-06-04). "Optimization with Extremal Dynamics". Physical Review Letters](https://arxiv.org/pdf/cond-mat/0010337.pdf).

The only twist w.r.t. classical EO is an affine invariant update equation for the worst performing solutions,

<!-- $
X \to X_1 + \frac{\beta}{\sqrt{2}} \cdot Z \cdot (X_2 - X_3), \,\,\,\, Z \sim \mathcal{N}(0,1)
$ --> <img style="transform: translateY(0.25em);" src="svg/eq.svg"/>

where X₁, X₂, X₃ are chosen random inside the pool of candidate solutions, this update mechanism allows EO to work on continuous spaces, and be invariant w.r.t. affine transformations of X and monotonous tranformations of the cost function.

## API:
```julia
function optimize(
    f,
    s,
    N;
    reps_per_particle = 100,
    β = 1.5,
    A = 1.0,
    atol = 0.0,
    rtol = sqrt(eps(1.0)),
    f_atol = 0.0,
    f_rtol = sqrt(eps(1.0)),
    verbose = false,
    rng = Random.GLOBAL_RNG,
    callback = state -> nothing,
)
```

- `f` : cost function to minimize, whose argument is either a scalar or a vector, must returns a scalar value.
- `s` : function whose input is the particle number and output is a random initial point to be ranked by `f`.
- `N` : number of particles to use, choose a number greater than `d+4` where `d` is the number of dimensions.
- `reps_per_particle` : maximum number of iterations per particle.

## Usage example:

```julia
using ExtremalOptimization
rosenbrock2d(x) = (x[1]-1)^2+(x[2]-x[1]^2)^2
initpoint(i) = randn(2)
optimize(rosenbrock2d, initpoint, 15)
```
output
```
(x = [0.99999998, 0.99999987], fx = 7.61e-15, f_nevals = 1172)
```
as expected the algorithm has found the optimum at `(1, 1)`, up to the specified tolerance.
