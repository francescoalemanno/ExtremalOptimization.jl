module ExtremalOptimization
    using Random, StatsBase

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

    function eostate(problem::EOProblem,N; γ, kw...)
        P = [problem.s(i) for i=1:N]
        C = [problem.f(P[i]) for i=1:N]
        order = sortperm(C)
        W = Weights([i^(-γ) for i=N:-1:1])
        best=order[1]
        P[best], C[best],  EOState(P,C,order,W,N)
    end

    function sample_distinct(rng,N,set,notset)
        N <= 0 && return ()
        @label resample
        s = rand(rng,set)
        (s ∈ notset) && @goto resample
        return (s, sample_distinct(rng,N-1,set,(notset...,s))...)
    end    

    function eostate(problem::EOProblem,S::EOState; β, rng, kw...)
        α = β*randn(rng)
        i = S.order[sample(rng,S.W)]
        j, k, q = sample_distinct(rng,3,1:S.N,(i,))
        nP = S.P[j] + α * (S.P[k] - S.P[q])
        nC = problem.f(nP)
        S.P[i] = nP
        S.C[i] = nC
        sortperm!(S.order,S.C)
        best = S.order[1]
        S.P[best], S.C[best], EOState(S.P,S.C,S.order,S.W,S.N)
    end

    function optimize(f,s,N;β=1.0,γ=1.2,rng=Random.GLOBAL_RNG, reps_per_particle=100, callback=state->nothing)
        problem=EOProblem(f,s)
        best, cost, state=eostate(problem,N;γ)
        callback(state)
        for i in 1:(N*reps_per_particle)
            best, cost, state=eostate(problem,state; β,rng)
            callback(state)
        end
        (x=best, fx=cost)
    end

    export optimize
end


