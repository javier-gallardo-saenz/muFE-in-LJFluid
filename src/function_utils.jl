using LinearAlgebra
using Distances
using StaticArrays

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
    return 4*params.ϵ*((aux)^(-12) - (aux)^(-6))
end


"""
Lennard-Jones potential derivative evaluation
"""
function LJ_pot_der(r::R, params::LJ_params) where {R<:Real}
    aux = (params.σ/r)^6
    return 24*params.ϵ*(-2*aux^(2) + aux)/r
end


"""
Softened interaction potential evaluation
"""
function soft_pot(r::R, params::LJ_params, λ::R) where {R<:Real}
    A = 0.5*(1-λ) + (r/params.σ)^6
    return 4*params.ϵ*(A^(-2) - A^(-1))
end


"""
Softened interaction potential derivative wrt r evaluation
"""
function soft_pot_dr(r::R, params::LJ_params, λ::R) where {R<:Real}
    aux = (r/params.σ)^6
    A = 0.5*(1-λ) + aux
    return 24*params.ϵ*aux*(-2*A^(-3) + A^(-2))/r
end


"""
Softened interaction potential derivative wrt λ evaluation
"""
function soft_pot_dλ(r::R, params::LJ_params, λ::R) where {R<:Real}
    aux = (r/params.σ)^6
    A = 0.5*(1-λ) + aux
    return 2*params.ϵ*(2*A^(-3) - A^(-2))
end


"""
Compute pairwise distances from a matrix R^{d x N} in which each column encodes the position of a particle 
Very memory efficient, but Distance.jl is better performance-wise if only distances are needed
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
Compute distances and distance vectors in a same pass
"""
function pairwise_dist_vec(q::AbstractMatrix{R}) where {R<:Real}
    d, N = size(q)
    T = typeof(norm(@view q[:,1]))
    num_pairs = Int(N*(N-1)/2)

    D = Vector{T}(undef, num_pairs)
    Δr = Vector{SVector{d,T}}(undef, num_pairs)

    idx = 1
    @inbounds for i in 1:N-1
        qi = @view q[:,]
        for j in i+1:N
            qj = @view q[:,j]
            Δr[idx] = qi .- qj
            D[idx] = sqrt(s)
            idx += 1
        end
    end
    return D 
end


"""
Particle instertion Hamiltonian. Versatile formulation that can be integrated with ForwardDiff, but not the most efficient. 
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
q-gradient of the particle insertion Hamiltonian.
Assumes q is given as Vector{SVector{d, R}} with length(q) = N
Periodic boundary conditions with the standard cutoff L/2 are imposed.
"""
function dH_pi_dq(q::Vector{SVector{R}}, params::LJ_params, λ::R, L::SVector{R}) where {R<:Real}
    d, N = length(q[1]), length(q)
    @assert d == length(L) "The boxes dimensions must match the dimension of q"
    T = typeof(norm(@view q[:,1]))

    dHdq = [@MVector zeros(T, d) for _ in 1:N]
    dq = MVector{d, T}(undef)

    @inbounds for i in 1:N-1
        #compute derivative of interaction potential wrt the inserted particle, fixed at the origin
        ri = sqrt(sum(abs2, q[i])) #no need to worry about minimum image convention here
        dHdq[i] .+= λ*LJ_pot_der(ri, params)*q[i]/ri

        for j in i+1:N
            #compute radius and distance vector between particles i and j
            dq .= q[i] .- q[j]
            dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
            r  = sqrt(sum(abs2, dq))

            #compute gradient
            Fij = LJ_pot_der(r, params)*dq/r
            dHdq[i] .+= Fij
            dHdq[j] .-= Fij
        end
    end

    rN = sqrt(sum(abs2, q[N]))
    dHdq[N] .+= λ*LJ_pot_der(rN, params)*q[N]/rN

    return dHdq
    
end


"""
q-gradient and λ-gradient of the particle insertion Hamiltonian.
Assumes q is given as Vector{SVector{d, R}} with length(q) = N
Periodic boundary conditions with the standard cutoff L/2 are imposed.
"""
function dH_pi_dq_and_dλ(q::Vector{SVector{R}}, params::LJ_params, λ::R, L::SVector{R}) where {R<:Real}
    d, N = length(q[1]), length(q)
    @assert d == length(L) "The boxes dimensions must match the dimension of q"
    T = typeof(norm(@view q[:,1]))

    dHdq = [@MVector zeros(T, d) for _ in 1:N]
    dHdλ = 0
    dq = MVector{d, T}(undef)

    @inbounds for i in 1:N-1
        #compute derivative of interaction potential with the inserted particle, fixed at the origin
        #and derivative of interaction potential wrt λ
        ri = sqrt(sum(abs2, q[i])) #no need to worry about minimum image convention here
        dHdq[i] .+= λ*LJ_pot_der(ri, params)*q[i]/ri
        dHdλ += LJ_pot(ri, params)

        for j in i+1:N
            #compute radius and distance vector between particles i and j
            dq .= q[i] .- q[j]
            dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
            r  = sqrt(sum(abs2, dq))

            #compute gradient
            Fij = LJ_pot_der(r, params)*dq/r
            dHdq[i] .+= Fij
            dHdq[j] .-= Fij
        end
    end

    rN = sqrt(sum(abs2, q[N]))
    dHdq[N] .+= λ*LJ_pot_der(rN, params)*q[N]/rN
    dHdλ += LJ_pot(rN, params)

    return dHdq, dHdλ
    
end


"""
q-gradient and U(q) of the particle insertion Hamiltonian, H = T(p) + U(q)
Assumes q is given as Vector{SVector{d, R}} with length(q) = N
Periodic boundary conditions with the standard cutoff L/2 are imposed.
"""
function u_and_dH_pi_dq(q::Vector{SVector{R}}, params::LJ_params, λ::R, L::SVector{R}) where {R<:Real}
    d, N = length(q[1]), length(q)
    @assert d == length(L) "The boxes dimensions must match the dimension of q"
    T = typeof(norm(@view q[:,1]))

    dHdq = [@MVector zeros(T, d) for _ in 1:N]
    U = 0
    dq = MVector{d, T}(undef)

    @inbounds for i in 1:N-1
        #compute derivative of interaction potential with the inserted particle, fixed at the origin
        #and contribution of the λ-term to V
        ri = sqrt(sum(abs2, q[i])) #no need to worry about minimum image convention here
        dHdq[i] .+= λ*LJ_pot_der(ri, params)*q[i]/ri
        U += λ*LJ_pot(ri, params)

        for j in i+1:N
            #compute radius and distance vector between particles i and j
            dq .= q[i] .- q[j]
            dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
            r  = sqrt(sum(abs2, dq))

            #compute gradient and contribution to V
            Fij = LJ_pot_der(r, params)*dq/r
            dHdq[i] .+= Fij
            dHdq[j] .-= Fij
            U += LJ_pot(r, params)
        end
    end

    rN = sqrt(sum(abs2, q[N]))
    dHdq[N] .+= λ*LJ_pot_der(rN, params)*q[N]/rN
    U += λ*LJ_pot(rN, params)

    return dHdq, U
    
end


"""
q-gradient, λ-gradient and U(q) of the particle insertion Hamiltonian, H = T(p) + U(q)
Assumes q is given as Vector{SVector{d, R}} with length(q) = N
Periodic boundary conditions with the standard cutoff L/2 are imposed.
"""
function U_and_dH_pi_dq_and_dλ(q::Vector{SVector{R}}, params::LJ_params, λ::R, L::SVector{R}) where {R<:Real}
    d, N = length(q[1]), length(q)
    @assert d == length(L) "The boxes dimensions must match the dimension of q"
    T = typeof(norm(@view q[:,1]))

    dHdq = [@MVector zeros(T, d) for _ in 1:N]
    dHdλ = 0
    U = 0
    dq = MVector{d, T}(undef)

    @inbounds for i in 1:N-1
        #compute derivative of interaction potential with the inserted particle, fixed at the origin
        #and contribution of the λ-term to V
        ri = sqrt(sum(abs2, q[i])) #no need to worry about minimum image convention here
        dHdq[i] .+= λ*LJ_pot_der(ri, params)*q[i]/ri
        U += λ*LJ_pot(ri, params)
        dHdλ += LJ_pot(ri, params)

        for j in i+1:N
            #compute radius and distance vector between particles i and j
            dq .= q[i] .- q[j]
            dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
            r  = sqrt(sum(abs2, dq))

            #compute gradient and contribution to V
            Fij = LJ_pot_der(r, params)*dq/r
            dHdq[i] .+= Fij
            dHdq[j] .-= Fij
            U += LJ_pot(r, params)
        end
    end

    rN = sqrt(sum(abs2, q[N]))
    dHdq[N] .+= λ*LJ_pot_der(rN, params)*q[N]/rN
    U += λ*LJ_pot(rN, params)
    dHdλ += LJ_pot(rN, params)

    return dHdq, dHdλ, U
    
end




"""
Softened particle insertion Hamiltonian. Versatile implementation that can be integrated with ForwardDiff, but not the most efficient.
Inserted particle is fixed at the origin during the process.
If stor_over_perf = true then pairwise distances are computed using the storage-optimized pairwise_dist.
Otherwise pairwise distances are computed using Distances.jl for optimized performance.
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


"""
q-gradient of the softened particle insertion Hamiltonian.
Assumes q is given as Vector{SVector{d, R}} with length(q) = N
Periodic boundary conditions with the standard cutoff L/2 are imposed.
"""
function dH_spi_dq(q::Vector{SVector{R}}, params::LJ_params, λ::R, L::SVector{R}) where {R<:Real}
    d, N = length(q[1]), length(q)
    @assert d == length(L) "The boxes dimensions must match the dimension of q"
    T = typeof(norm(@view q[:,1]))

    dHdq = [@MVector zeros(T, d) for _ in 1:N]
    dq = MVector{d, T}(undef)

    @inbounds for i in 1:N-1
        #compute derivative of interaction potential wrt the inserted particle, fixed at the origin
        r0 = sqrt(sum(abs2, q[i])) #no need to worry about minimum image convention here
        dHdq[i] .+= λ*soft_pot_dr(r0, params, λ)*q[i]/r0

        for j in i+1:N
            #compute radius and distance vector between particles i and j
            dq .= q[i] .- q[j]
            dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
            r  = sqrt(sum(abs2, dq))

            #compute gradient
            Fij = LJ_pot_der(r, params)*dq/r
            dHdq[i] .+= Fij
            dHdq[j] .-= Fij
        end
    end

    rN = sqrt(sum(abs2, q[N]))
    dHdq[N] .+= λ*soft_pot_dr(rN, params, λ)*q[N]/rN

    return dHdq
    
end


"""
q-gradient and λ-gradient of the softened particle insertion Hamiltonian.
Assumes q is given as Vector{SVector{d, R}} with length(q) = N
Periodic boundary conditions with the standard cutoff L/2 are imposed.
"""
function dH_spi_dq_and_dλ(q::Vector{SVector{R}}, params::LJ_params, λ::R, L::SVector{R}) where {R<:Real}
    d, N = length(q[1]), length(q)
    @assert d == length(L) "The boxes dimensions must match the dimension of q"
    T = typeof(norm(@view q[:,1]))

    dHdq = [@MVector zeros(T, d) for _ in 1:N]
    dHdλ = 0
    dq = MVector{d, T}(undef)

    @inbounds for i in 1:N-1
        #compute derivative of interaction potential with the inserted particle, fixed at the origin
        #and derivative of interaction potential wrt λ
        ri = sqrt(sum(abs2, q[i])) #no need to worry about minimum image convention here
        dHdq[i] .+= λ*soft_pot_dr(ri, params, λ)*q[i]/r0
        dHdλ += soft_pot(ri, params, λ) + soft_pot_dλ(ri, params, λ)

        for j in i+1:N
            #compute radius and distance vector between particles i and j
            dq .= q[i] .- q[j]
            dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
            r  = sqrt(sum(abs2, dq))

            #compute gradient
            Fij = LJ_pot_der(r, params)*dq/r
            dHdq[i] .+= Fij
            dHdq[j] .-= Fij
        end
    end

    rN = sqrt(sum(abs2, q[N]))
    dHdq[N] .+= λ*soft_pot_dr(rN, params, λ)*q[N]/rN
    dHdλ += soft_pot(rN, params, λ) + soft_pot_dλ(rN, params, λ)

    return dHdq, dHdλ
    
end


"""
q-gradient and U(q) of the softened particle insertion Hamiltonian, H = T(p) + U(q)
Assumes q is given as Vector{SVector{d, R}} with length(q) = N
Periodic boundary conditions with the standard cutoff L/2 are imposed.
"""
function U_and_dH_spi_dq(q::Vector{SVector{R}}, params::LJ_params, λ::R, L::SVector{R}) where {R<:Real}
    d, N = length(q[1]), length(q)
    @assert d == length(L) "The boxes dimensions must match the dimension of q"
    T = typeof(norm(@view q[:,1]))

    dHdq = [@MVector zeros(T, d) for _ in 1:N]
    U = 0
    dq = MVector{d, T}(undef)

    @inbounds for i in 1:N-1
        #compute derivative of interaction potential with the inserted particle, fixed at the origin
        #and contribution of the λ-term to V
        ri = sqrt(sum(abs2, q[i])) #no need to worry about minimum image convention here
        dHdq[i] .+= λ*soft_pot_dr(ri, params, λ)*q[i]/ri
        U += λ*soft_pot(r0, params, λ)

        for j in i+1:N
            #compute radius and distance vector between particles i and j
            dq .= q[i] .- q[j]
            dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
            r  = sqrt(sum(abs2, dq))

            #compute gradient and contribution to V
            Fij = LJ_pot_der(r, params)*dq/r
            dHdq[i] .+= Fij
            dHdq[j] .-= Fij
            U += LJ_pot(r, params)
        end
    end

    rN = sqrt(sum(abs2, q[N]))
    dHdq[N] .+= λ*soft_pot_dr(rN, params, λ)*q[N]/rN
    U += λ*soft_pot(rN, params, λ)

    return dHdq, U
    
end


"""
q-gradient, λ-gradient and U(q) of the softened particle insertion Hamiltonian, H = T(p) + U(q)
Assumes q is given as Vector{SVector{d, R}} with length(q) = N
Periodic boundary conditions with the standard cutoff L/2 are imposed.
"""
function U_and_dH_pi_dq_and_dλ(q::Vector{SVector{R}}, params::LJ_params, λ::R, L::SVector{R}) where {R<:Real}
    d, N = length(q[1]), length(q)
    @assert d == length(L) "The boxes dimensions must match the dimension of q"
    T = typeof(norm(@view q[:,1]))

    dHdq = [@MVector zeros(T, d) for _ in 1:N]
    dHdλ = 0
    U = 0
    dq = MVector{d, T}(undef)

    @inbounds for i in 1:N-1
        #compute derivative of interaction potential with the inserted particle, fixed at the origin
        #and contribution of the λ-term to V
        ri = sqrt(sum(abs2, q[i])) #no need to worry about minimum image convention here
        dHdq[i] .+= λ*soft_pot_dr(ri, params, λ)*q[i]/ri
        dHdλ += soft_pot(ri, params, λ) + soft_pot_dλ(ri, params,λ)
        U += λ*soft_pot(r0, params, λ)

        for j in i+1:N
            #compute radius and distance vector between particles i and j
            dq .= q[i] .- q[j]
            dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
            r  = sqrt(sum(abs2, dq))

            #compute gradient and contribution to V
            Fij = LJ_pot_der(r, params)*dq/r
            dHdq[i] .+= Fij
            dHdq[j] .-= Fij
            U += LJ_pot(r, params)
        end
    end

    rN = sqrt(sum(abs2, q[N]))
    dHdq[N] .+= λ*soft_pot_dr(rN, params, λ)*q[N]/rN
    dHdλ += soft_pot(rN, params, λ) + soft_pot_dλ(rN, params,λ)
    U += λ*soft_pot(rN, params, λ)

    return dHdq, dHdλ, U
    
end



