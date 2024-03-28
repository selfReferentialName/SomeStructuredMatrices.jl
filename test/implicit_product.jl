using Random

@testset "ImplicitProduct" begin
	A = rand(10, 10)
	B = rand(10, 100)
	C = rand(100, 10)
	ABCt = A * B * C
	ABCi = ImplicitProduct([A, B, C])
	ABCo = optimise(ABCi)
	x = rand(10)
	@test ABCt*x ≈ A*(B*(C*x))
	@test ABCi*x == A*(B*(C*x))
	@test ABCo*x ≈ A*(B*(C*x))
	@test size(ABCi) == (10,10)
	@test size(ABCi,1) == 10
	@test size(ABCi,2) == 10
end
