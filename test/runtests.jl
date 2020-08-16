using ExtremalOptimization
using Test, Random

@testset "ExtremalOptimization.jl" begin
    sol = optimize(x->(x-1)^2,i->randn(),50)
    @test sol.x ≈ 1.0
    sol = optimize(x->(x[1]-1)^2+(x[2]-x[1]^2)^2,i->randn(2),50,reps_per_particle=200)
    @test all(sol.x .≈ 1)
end
