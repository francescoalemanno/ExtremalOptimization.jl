using ExtremalOptimization
using Test, Random

@testset "ExtremalOptimization.jl" begin
    sol = optimize(x -> (x - 1)^2, i -> randn(), 50)
    @test sol.x ≈ 1 rtol = 1e-4
    sol = optimize(
        x -> (x[1] - 1)^2 + (x[2] - x[1]^2)^2,
        i -> randn(2),
        50,
        reps_per_particle = 200,
    )
    map(sol.x) do x
        @test x ≈ 1 rtol = 1e-4
    end
end
