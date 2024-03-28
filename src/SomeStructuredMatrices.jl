module SomeStructuredMatrices

# export LinearMap
#export optimise!, optimize!
export optimise, optimize

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

optimize = optimise

end
