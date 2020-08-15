module ExtremalOptimization
    using Random, StatsBase

    function optimize(f,s,N;γ=1.2,β=1.0,rng=Random.GLOBAL_RNG,reps_per_particle=100)
        P = [s(i) for i=1:N]
        C = [f(P[i]) for i=1:N]
        order = sortperm(C)
        W = Weights([i^(-γ) for i=N:-1:1])
        for reps = 1:(reps_per_particle*N)
            i = order[sample(rng,W)]
            j = i
            while j==i; j=rand(rng,1:N); end
            α=(1 + β * abs2(randn(rng)))/2
            nP=P[i]+α*(P[j]-P[i])
            nC=f(nP)
            P[i]=nP
            C[i]=nC
            sortperm!(order,C)
        end
        best=order[begin]
        (x=P[best],cost=C[best],pop=P[order],costs=C[order])
    end
    export optimize
end
