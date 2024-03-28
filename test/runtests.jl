using SomeStructuredMatrices
using Test
using Aqua
using LinearAlgebra
using Random

@testset "SomeStructuredMatrices.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(SomeStructuredMatrices, ambiguities=false)
    end
	include("implicit_product.jl")
	@testset "Trivial Operations" begin
		@test optimise(I) == I
		@test optimize(I) == I
	end
end
