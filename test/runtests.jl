using SomeStructuredMatrices
using Test
using Aqua

@testset "SomeStructuredMatrices.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(SomeStructuredMatrices, ambiguities=false)
    end
    # Write your tests here.
end
