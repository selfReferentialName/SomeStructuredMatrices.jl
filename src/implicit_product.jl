"""
	ImplicitProduct

Lazy matrix-matrix multiplication. Handles an arbitrary number of matrices
multiplied together, allowing for matrix-vector multiplication by doing a
sequence of multiplications with the contained vectors. Useful for low rank
approximations.
"""
struct ImplicitProduct
	mats :: Vector
	function ImplicitProduct(mats :: Vector)
		if length(mats) < 1
			throw(DomainError("Cannot make a matrix out of zero multiplications."))
		end
		for i=2:length(mats)
			if size(mats[i-1],2) != size(mats[i],1)
				throw(DimensionMismatch("Cannot implicitly multiply matrices of size " * string(size(mats[i-1])) * " and " * string(size(mats[i]))))
			end
		end
		if all(A->typeof(A)!=ImplicitProduct, mats)
			return new(mats)
		else
			rmats = Vector{LinearMap{eltype(mats[1])}}(undef, 0)
			for i=1:length(mats)
				if typeof(mats[i]) == ImplicitProduct
					for mat in mats[i].mats
						push!(rmats, mat)
					end
				else
					push!(rmats, mats[i])
				end
			end
			return new(rmats)
		end
	end
end

function Base.:*(A :: ImplicitProduct, x :: AbstractVector)
	y = x
	for i=length(A.mats):-1:1
		y = A.mats[i] * y
	end
	return y
end
function Base.:*(A :: ImplicitProduct ...)
	get_mats(M) = M.mats
	return ImplicitProduct(vcat(get_mats.(A)...))
end
function Base.:*(A :: ImplicitProduct, B...)
	return ImplicitProduct(vcat(A.mats, B))
end
function Base.:*(A, B :: ImplicitProduct)
	return ImplicitProduct(vcat([A], B.mats))
end

function LinearAlgebra.transpose(A :: ImplicitProduct)
	return ImplicitProduct(transpose.(A.mats))
end

function Base.eltype(A :: ImplicitProduct)
	if length(A.mats) == 0
		return Any
	else
		return Union{eltype.(A)...}
	end
end

function Base.size(A :: ImplicitProduct)
	return (size(A.mats[1],1), size(A.mats[end],2))
end
function Base.size(A :: ImplicitProduct, dim :: Int)
	if dim == 1
		return size(A.mats[1],1)
	elseif dim == 2
		return size(A.mats[end],2)
	else
		throw(BoundsError("Dimension " * string(dim) * " is not a dimension matrices have."))
	end
end

"""
	optimise(A :: ImplicitProduct; max_flops_per_mul=1000000000, min_stage_size=50)

Try to reduce the number of matrices in order to make matrix-vector multiplication take as few
floating point operations as possible. It does this by matrix multiplications, which it
assumes take ``nmk`` (FMA) FLOPS to multiply a n x m matrix by a m x k one (with vectors having ``k=1``).
It does no matrix multiplications predicted to take over `max_flops_per_mul` steps, and will
always do matrix multiplications if the resulting matrix is smaller than or the same size as
`min_stage_size` x `min_stage_size`.

As `ImplicitProduct`s are immutable, this operates out of place, returning the optimised matrix.
"""
function optimise(A :: ImplicitProduct; max_flops_per_mul=1000000000, min_stage_size=50)
	# step 1: figure out optimal matrix multiplication sequence, assuming we multiply everything out
	# https://en.wikipedia.org/wiki/Matrix_chain_multiplication
	cost = fill(Inf64, (length(A.mats), length(A.mats))) # cost[i,j] is min flop cost of A.mats[i]*...*A.mats[j]
	splitat = fill(Int64(0), (length(A.mats), length(A.mats))) # splitat[i,j] = k means do (A.mats[i]*...*A.mats[k]) * (A.mats[k+1]*...*A.mats[j])
	for i=1:length(A.mats)
		cost[i,i] = 0
	end
	for len=2:length(A.mats)
		for i=1:(length(A.mats)-len+1)
			j = i + len - 1
			for k=i:(j-1)
				my_cost = cost[i,k] + cost[k+1,j]
				if my_cost < cost[i,j]
					cost[i,j] = my_cost
					splitat[i,j] = k
				end
			end
		end
	end

	# step 2: carry out the multiplications
	do_mult(A, B) = (size(A,1) * size(A,2) * size(B,2) <= max_flops_per_mul &&
			size(A,1)*size(B,2)^2 <= size(A,1)*size(A,2)^2 + size(A,2)*size(B,2)^2) ||
			(size(A,1) <= min_stage_size && size(B,2) <= min_stage_size)
	function multiply_through(i, j)
		if i == j
			return A.mats[i]
		else
			my_split = splitat[i, j]
			L = multiply_through(i, my_split)
			R = multiply_through(my_split+1, j)
			if do_mult(L, R)
				return L*R
			else
				return ImplicitProduct([L, R])
			end
		end
	end
	return multiply_through(1, length(A.mats))
end
