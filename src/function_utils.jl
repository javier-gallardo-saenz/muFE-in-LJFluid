using LinearAlgebra

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
Tried to make it efficient
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
Particle instertion Hamiltonian, vectorized for efficiency
Inserted particle is fixed at the origin during the process
"""
function H_pi(p::AbstractMatrix{R}, q::AbstractMatrix{R}, m::R, params::LJ_params, λ::R) where {R<:Real}
    pw_dist = pairwise_dist(q)
    K = sum((norm.(eachcol(p))).^2)/(2*m)
    V_orig = sum(LJ_pot.(pw_dist, params))  
    V_inserted = λ*sum(LJ_pot.(norm.(eachcol(q)), params)) 
    return K + V_orig + V_inserted
end


"""
Softened particle insertion Hamiltonian, vectorized for efficiency
Inserted particle is fixed at the origin during the process
"""
function H_spi(p::AbstractMatrix{R}, q::AbstractMatrix{R}, m::R, params::LJ_params, λ::R) where {R<:Real}
    pw_dist = pairwise_dist(q)
    K = sum((norm.(eachcol(p))).^2)/(2*m)
    V_orig = sum(LJ_pot.(pw_dist, params))  
    V_inserted = λ*sum(soft_pot.(norm.(eachcol(q)), params)) 
    return K + V_orig + V_inserted
end