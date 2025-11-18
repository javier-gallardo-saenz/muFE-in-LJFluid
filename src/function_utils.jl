using LinearAlgebra
using Distances

"""
Parameters for Lennard-Jones potential
"""
struct LJ_params
    σ::Float64
    ϵ::Float64
end


"""
Lennard-Jones potential evaluation
"""
function LJ_pot(r::R, params::LJ_params) where {R<:Real}
    aux = r/params.σ
    return 4*params.ϵ*((aux)^12 - (aux)^6)
end


"""
Softened interaction potential
"""
function soft_pot(r::R, params::LJ_params, λ::R) where {R<:Real}
    A = 0.5*(1-λ) + (r/params.σ)^6
    return 4*params.ϵ*(A^(-2) - A^(-1))
end


"""
Compute pairwise distances from a matrix R^{d x N} in which each column encodes the position of a particle 
Very memory efficient, but Distance.jl is better performance-wise
"""
function pairwise_dist(q::AbstractMatrix{R}) where {R<:Real}
    d, N = size(q)
    T = typeof(norm(@view q[:,1]))
    D = Vector{T}(undef, Int(N*(N-1)/2))
    idx = 1
    @inbounds for i in 1:N
        qi = @view q[:,i]
        for j in i+1:N
            qj = @view q[:,j]
            s = zero(T)
            @simd for k in 1:d
                aux = qi[k] - qj[k]
                s += aux * aux
            end
            D[idx] = sqrt(s)
            idx += 1
        end
    end
    return D 
end


"""
Compute pairwise distances from a matrix R^{d x N} in which each column encodes the position of a particle,
 given a certain metric (using Distances.jl)
"""
function pairwise_dist_perf(metric::Metric, q::AbstractMatrix{R}) where {R <: Real}
    N = size(q, 2)
    q1 = @view q[:,1]
    T = eltype(colwise(metric, q1, @view q[:, 2:2]))
    D = Vector{T}(undef, Int(N*(N-1)/2)) 
    idx = 1
    @inbounds for i in 1:N-1
        qi = @view q[:, i]
        aux_dists = colwise(metric, qi, @view q[:, i+1:end])
        D[idx:idx + length(aux_dists) - 1] = aux_dists
        idx += length(aux_dists)
    end

    return D
end


"""
Particle instertion Hamiltonian, vectorized for efficiency
Inserted particle is fixed at the origin during the process
If stor_over_perf = true then pairwise distances are computed using the storage-optimized pairwise_dist
Otherwise pairwise distances are computed using Distances.jl for optimized performance
Periodic boundary conditions with the standard cutoff L/2 can be imposed using the PeriodicEuclidean(period) metric.
"""
function H_pi(p::AbstractMatrix{R}, q::AbstractMatrix{R}, m::R, params::LJ_params, λ::R,
     stor_over_perf::Bool = false, metric::Metric = Euclidean()) where {R<:Real}
    if stor_over_perf
        pw_dist = pairwise_dist(q)
    else
        pw_dist = pairwise_dist_perf(metric, q, dims=2)
    end
    K = sum((norm.(eachcol(p))).^2)/(2*m)
    V_orig = sum(LJ_pot.(pw_dist, Ref(params)))  
    V_inserted = λ*sum(LJ_pot.(norm.(eachcol(q)), Ref(params))) 
    return K + V_orig + V_inserted
end


"""
Softened particle insertion Hamiltonian, vectorized for efficiency
Inserted particle is fixed at the origin during the process
If stor_over_perf = true then pairwise distances are computed using the storage-optimized pairwise_dist
Otherwise pairwise distances are computed using Distances.jl for optimized performance
Periodic boundary conditions with the standard cutoff L/2 can be imposed using the PeriodicEuclidean(period) metric.
"""
function H_spi(p::AbstractMatrix{R}, q::AbstractMatrix{R}, m::R, params::LJ_params, λ::R, 
    stor_over_perf::Bool = false, metric::Metric = Euclidean()) where {R<:Real}
    if stor_over_perf
        pw_dist = pairwise_dist(q)
    else
        pw_dist = pairwise(metric, q, dims=2)
    end
    K = sum((norm.(eachcol(p))).^2)/(2*m)
    V_orig = sum(LJ_pot.(pw_dist, Ref(params)))  
    V_inserted = λ*sum(soft_pot.(norm.(eachcol(q)), Ref(params), Ref(λ))) 
    return K + V_orig + V_inserted
end
