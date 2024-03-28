struct ImplicitProduct{T} <: TransposeableLinearMap{T}
	mats :: Vector{LinearMap{T}}
	function ImplicitProduct(mats :: Vector{LinearMap})
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

function Base.*(A :: ImplicitProduct, x :: AbstractVector)
	y = x
	for i=1:length(A.mats)
		y = A.mats[i] * y
	end
	return y
end
function Base.*(A... :: Vector{ImplicitProduct})
	get_mats(M) = M.mats
	return ImplicitProduct(vcat(get_mats.(A)...))
end
function Base.*(A :: ImplicitProduct, B... :: Vector{LinearMap})
	return ImplicitProduct(vcat(A.mats, B))
end
function Base.*(A :: LinearMap, B :: ImplicitProduct)
	return ImplicitProduct(vcat([A], B.mats))
end

function LinearAlgebra.transpose(A :: LinearMap)
	return ImplicitProduct(transpose.(A.mats))
end

function Base.eltype(A :: LinearMap{T})
	return T
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
	splitat = fill(Int64(0), (length(A.mats))) # splitat[i,j] = k means do (A.mats[i]*...*A.mats[k]) * (A.mats[k+1]*...*A.mats[j])
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
	whereright = BitVector(undef, length(A.mats))
	whereright .= false
	splits = Vector{Int64}(undef, length(A.mats))
	splits[1] = splitat[1,end]
	bounds = Vector{Tuple{Int64,Int64}}(undef, length(A.mats))
	bounds[1] = (1, length(A.mats))
	depth = 1
	going_up = false
	cur_mat = nothing
	left_mats = Vector{LinearMap{eltype}}
	do_mult(A, B) = (size(A,1) * size(A,2) * size(B,2) <= max_flops_per_mul) ||
			(size(A,1) <= min_stage_size && size(B,2) <= min_stage_size)
	while true
		if splits[depth] == bounds[depth][whereright[depth] ? 2 : 1] # at an end
			cur_mat = A.mats[spltis[depth]]
			depth -= 1
			going_up = true
		elseif depth == 0 # gone up to the top
			return cur_mat
		elseif going_up && whereright[depth]
			if do_mult(left_mats[depth], cur_mat)
				cur_mat = left_mats[depth] * cur_mat
			else
				cur_mat = ImplicitProduct([left_mats, cur_mat])
			end
			depth -= 1
			going_up = true
		elseif going_up # && !whereright[depth]
			left_mats[depth] = cur_mat
			depth += 1
			whereright[depth] = true
			going_up = false
		else # going down, not at end
			if whereright[depth]
				splits[depth+1] = splitat[splits[depth], bounds[depth][2]]
				bounds[depth+1] = (splits[depth], bounds[depth][2])
			else
				splits[depth+1] = splitat[bounds[depth][1], splits[depth]]
				bounds[depth+1] = (bounds[depth][1], splits[depth])
			end
			depth += 1
			whereright[depth] = false
			going_up = false
		end
	end
end
