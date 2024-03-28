module SomeStructuredMatrices

# export LinearMap
#export optimise!, optimize!
export optimise, optimize
export ImplicitProduct

using LinearAlgebra

#"""
#	LinearMap{T}
#
#A generalisation of matrices on ``T^n``, defined by supporting multiplication.
#"""
#abstract type LinearMap{T} end

#AbstractMatrix{T} <: LinearMap{T} where T
#AbstractQ{T} <: LinearMap{T} where T

include("implicit_product.jl")

"""
	optimise(A; kwargs...)

Copy A. This is here to help make optimise(A) safe to call for all A.
"""
function optimise(A; kwargs...)
	return copy(A)
end

"""
	optimize(A; kwargs...)

An alias for optimise.
"""
optimize = optimise

end
