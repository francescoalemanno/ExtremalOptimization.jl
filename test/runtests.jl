using ExtremalOptimization
using Test

@testset "ExtremalOptimization.jl" begin
    sol = optimize(x->(x-1)^2,i->randn(),50)
    @test sol.x ≈ 1.0
    @test sol.cost ≈ 0.0
end
