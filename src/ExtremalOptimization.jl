module ExtremalOptimization
using Random, StatsBase, LinearAlgebra

struct EOState{vP,vC,vO,tW,tN}
    P::Vector{vP}
    C::Vector{vC}
    order::Vector{vO}
    W::tW
    N::tN
end

struct EOProblem{tF,tS}
    f::tF
    s::tS
end

function approxequality(a, b, atol, rtol)
    Δ = norm(a - b)
    M = hypot(norm(a), norm(b))
    Δ <= max(atol, M * rtol)
end

function hasconverged(state::EOState, atol, rtol, f_atol, f_rtol)
    f_a = state.C[state.order[1]]
    f_b = state.C[state.order[2]]
    sol_a = state.P[state.order[1]]
    sol_b = state.P[state.order[2]]
    approxequality(f_a, f_b, f_atol, f_rtol) || approxequality(sol_a, sol_b, atol, rtol)
end

function eostate(problem::EOProblem, N; τ, kw...)
    P = [problem.s(i) for i = 1:N]
    C = [problem.f(P[i]) for i = 1:N]
    for i in eachindex(P)
        isfinite(C[i]) && continue
        error("The cost function has produced non-finite values for the initial point: $(P[i])")
    end
    order = sortperm(C)
    wS = N^(-τ)
    W = Weights([i^(-τ) - wS for i = N:-1:1])
    best = order[1]
    P[best], C[best], EOState(P, C, order, W, N)
end

function sample_distinct(rng, N, set, notset)
    N <= 0 && return ()
    @label resample
    s = rand(rng, set)
    (s ∈ notset) && @goto resample
    return (s, sample_distinct(rng, N - 1, set, (notset..., s))...)
end

function eostate(problem::EOProblem, S::EOState; β, rng, verbose, kw...)
    i = S.order[sample(rng, S.W)]
    @label regen
    j, k, q = sample_distinct(rng, 3, 1:S.N, (i,))
    nP = S.P[j] + (β/sqrt(2)) * randn(rng) * (S.P[k] - S.P[q])
    nC = problem.f(nP)
    isfinite(nC) || begin
        verbose && println("Failed function evaluation at $nP")
        @goto regen
    end
    S.P[i] = nP
    S.C[i] = nC
    sortperm!(S.order, S.C)
    best = S.order[1]
    S.P[best], S.C[best], EOState(S.P, S.C, S.order, S.W, S.N)
end
"""
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
optimize(rosenbrock2d, initpoint, 50)
```
output
```
(x = [1.000000001, 1.000000004], fx = 4.0e-18, f_nevals = 2726)
```
as expected the algorithm has found the optimum at `(1, 1)`, up to the specified tolerance.
"""
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
    problem = EOProblem(f, s)
    τ = 1 + A / log(N)
    best, cost, state = eostate(problem, N; τ)
    callback(state)
    fnevals = N
    for i = 1:(N*reps_per_particle)
        best, cost, state = eostate(problem, state; β, rng, verbose)
        fnevals += 1
        callback(state)
        verbose && println(i, " ", cost)
        hasconverged(state, atol, rtol, f_atol, f_rtol) && break
    end
    (x = best, fx = cost, f_nevals = fnevals)
end

export optimize
end
