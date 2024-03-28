@testset "ImplicitProduct" begin
	A = rand(10, 10)
	B = rand(10, 100)
	C = rand(100, 10)
	ABCt = A * B * C
	ABCi = ImplicitProduct([A, B, C]) :: ImplicitProduct
	ABCo = optimise(ABCi)
	ABCm = A * (ImplicitProduct([B]) * C)
	ABCm2 = A * ImplicitProduct([B]) * C
	ABCm3 = ImplicitProduct([A]) * ImplicitProduct([B]) * ImplicitProduct([C])
	ABCc = ImplicitProduct([ImplicitProduct([A]), B, C])
	x = rand(10)
	@test ABCt*x ≈ A*(B*(C*x))
	@test ABCi*x == A*(B*(C*x))
	@test ABCo*x ≈ A*(B*(C*x))
	@test ABCi' * x ≈ ABCt' * x
	@test transpose(ABCi) * x ≈ transpose(ABCt) * x
	@test size(ABCi) == (10,10)
	@test size(ABCi,1) == 10
	@test size(ABCi,2) == 10
	@test ABCo ≈ optimise(ABCm)
	@test ABCo ≈ optimise(ABCm2)
	@test ABCo ≈ optimise(ABCc)
	@test ABCo ≈ optimise(ABCm3)
	@test eltype(ABCi) == Float64

	nbig = Int(2e3)
	rbig = Int(2e3)
	Abig = rand(nbig, rbig)
	Bbig = rand(rbig, nbig)
	ABbig = ImplicitProduct([Abig, Bbig])
	xbig = rand(nbig)
	@test typeof(optimise(ABbig)) == ImplicitProduct
	@test optimise(ABbig) * xbig ≈ ABbig * xbig
end

@testset "ImplicitProduct throws exceptions correctly" begin
	@test_throws DomainError ImplicitProduct([])
	@test_throws BoundsError size(ImplicitProduct([rand(10,20)]), 3)
	@test_throws BoundsError size(ImplicitProduct([rand(10,20)]), 0)
	@test_throws BoundsError size(ImplicitProduct([rand(10,20)]), -1)
	@test_throws DimensionMismatch ImplicitProduct([rand(2,3), rand([2,3])])
end
