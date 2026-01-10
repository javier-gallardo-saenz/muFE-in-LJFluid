using LinearAlgebra
using Distances
using StaticArrays

"""
Ideal term of the chemical energy. NEEDS DIMENSION-FULL UNITS
"""
function ideal_mu_term(kBT::Float64, m::Float64, L::SVector{d,Float64}, N::Int64; h = 6.62607015e-34) where {d}
    V = prod(L)
    Λ = h / sqrt(2π * m * kBT)
    return -kbt * log((V / Λ^d) / (N + 1))
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
λ*Lennard Jones
For use in Hamiltonian evaluation functions
"""
function λLJ(r::R, params::LJ_params, λ::R) where {R<:Real}
    return λ*LJ_pot(r, params)
end

"""
Derivative of λ*Lennard Jones wrt r
For use in Hamiltonian evaluation functions
"""
function dλLJ_dr(r::R, params::LJ_params, λ::R) where {R<:Real}
    return λ*LJ_pot_der(r, params)
end


"""
Derivative of λ*Lennard Jones wrt λ
For use in Hamiltonian evaluation functions
"""
function dλLJ_dλ(r::R, params::LJ_params, λ::R) where {R<:Real}
    return LJ_pot(r, params)
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
λ*Softened potential
For use in Hamiltonian evaluation functions
"""
function λsoft_pot(r::R, params::LJ_params, λ::R) where {R<:Real}
    return λ*soft_pot(r,params,λ)
end


"""
Derivative of λ*Softened potential wrt r 
For use in Hamiltonian evaluation functions
"""
function dλsoft_pot_dr(r::R, params::LJ_params, λ::R) where {R<:Real}
    return λ*soft_pot_dr(r, params, λ)
end


"""
Derivative of λ*Softened potential wrt λ
For use in Hamiltonian evaluation function
"""
function dλsoft_pot_dλ(r::R, params::LJ_params, λ::R) where {R<:Real}
    return soft_pot(r, params, λ) + λ*soft_pot_dλ(r, params, λ)
end


"""
Tail correction for μ_ext for a Lennard Jones system
"""
function tailcorr_μex_LJ(ρ::Real, rc::Real, params::LJ_params)
    #tail contribution is (2π*ρ)*int_{rc}^{∞}{4ϵ*((σ/r)^12 - (σ/r)^6)*r^2}
    int_r12 = (params.σ)^(12)*rc^(-9)/9.0
    int_r6 = (params.σ)^(6)*rc^(-3)/3.0
    return 8π*ρ*params.ϵ*(int_r12 - int_r6)
end


"""
Tail correction for average potential energy in a Lennard Jones system
We assume g(r) = 1 ∀r > rc
"""
function tailcorr_U_LJ(ρ::Real, N::Real, rc::Real, params::LJ_params)
    return N*tailcorr_μex_LJ(ρ, rc, params)/2.0
end



"""
q-gradient of the particle insertion Hamiltonian.
Assumes q is given as Vector{SVector{d, R}} with length(q) = N
Periodic boundary conditions with the standard cutoff L/2 are imposed.
Input functions V, dVdr must have signature r, LJ_params
Input functions λV, dλVdr and dλVdλ must have signature r, LJ_params, λ
"""
function dHpi_dq!(dHdq::Vector{MVector{d,Float64}}, q::Vector{<:MVector{d,Float64}}, params::LJ_params, λ::Real, L::SVector{d,Float64},
    dVdr::Function, dλVdr::Function) where {d}

    N = length(q)
    @assert N == length(dHdq) "q and dHdq must have the same number of particles"

    T = eltype(dHdq[1])
    dq = MVector{d, T}(undef)
    @inbounds for i in 1:N
        dHdq[i] .= zero(T)
    end
    
    @inbounds for i in 2:N
        #compute derivative of interaction potential wrt the inserted particle, which is labeled by 1
        dq .= q[1] .- q[i]
        dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
        r  = sqrt(sum(abs2, dq))
        F1i = dλVdr(r, params, λ)*dq/r
        dHdq[1] .+= F1i
        dHdq[i] .-= F1i

        for j in i+1:N
            #compute radius and distance vector between particles i and j
            dq .= q[i] .- q[j]
            dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
            r  = sqrt(sum(abs2, dq))

            #compute gradient
            Fij = dVdr(r, params)*dq/r
            dHdq[i] .+= Fij
            dHdq[j] .-= Fij
        end
    end

    return nothing
end


"""
q-gradient and λ-gradient of the softened particle insertion Hamiltonian, H = T(p) + U(q)
Assumes q is given as Vector{MVector{d, R}} with length(q) = N
Periodic boundary conditions with the standard cutoff L/2 are imposed.
Input functions V, dVdr must have signature r, LJ_params
Input functions λV, dλVdr and dλVdλ must have signature r, LJ_params, λ
"""
function dHpi_dq_and_dλ!(dHdq::Vector{MVector{d,Float64}}, q::Vector{MVector{d, Float64}}, params::LJ_params, λ::Real, L::SVector{d, Float64},
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function) where {d}

    N = length(q)
    @assert N == length(dHdq) "q and dHdq must have the same number of particles"

    T = eltype(dHdq[1])
    dq = MVector{d, T}(undef)
    @inbounds for i in 1:N
        dHdq[i] .= zero(T)
    end
    dHdλ = 0

    @inbounds for i in 2:N
        #compute derivative of interaction potential wrt the inserted particle, which is labeled by 1
        dq .= q[1] .- q[i]
        dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
        r  = sqrt(sum(abs2, dq))
        F1i = dλVdr(r, params, λ)*dq/r
        dHdq[1] .+= F1i
        dHdq[i] .-= F1i
        dHdλ += dλVdλ(r, params, λ)

        for j in i+1:N
            #compute radius and distance vector between particles i and j
            dq .= q[i] .- q[j]
            dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
            r  = sqrt(sum(abs2, dq))

            #compute gradient and contribution to V
            Fij = dVdr(r, params)*dq/r
            dHdq[i] .+= Fij
            dHdq[j] .-= Fij
        end
    end

    return dHdλ
end


"""
q-gradient and λU(q)|λ,λ_next of the softened particle insertion Hamiltonian, H = T(p) + U(q) + λU(q)
Assumes q is given as Vector{MVector{d, R}} with length(q) = N
Periodic boundary conditions with the standard cutoff L/2 are imposed.
Input functions V, dVdr must have signature r, LJ_params
Input functions λV, dλVdr and dλVdλ must have signature r, LJ_params, λ
"""
function twoU_and_dHpi_dq!(dHdq::Vector{MVector{d,Float64}}, q::Vector{MVector{d, Float64}}, params::LJ_params, λ::Real, λ_next::Real, L::SVector{d, Float64},
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function) where {d}

    N = length(q)
    @assert N == length(dHdq) "q and dHdq must have the same number of particles"

    T = eltype(dHdq[1])
    dq = MVector{d, T}(undef)
    @inbounds for i in 1:N
        dHdq[i] .= zero(T)
    end
    λU = 0
    λnextU = 0

    @inbounds for i in 2:N
        #compute derivative of interaction potential wrt the inserted particle, which is labeled by 1
        dq .= q[1] .- q[i]
        dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
        r  = sqrt(sum(abs2, dq))
        F1i = dλVdr(r, params, λ)*dq/r
        dHdq[1] .+= F1i
        dHdq[i] .-= F1i
        λU += λV(r, params, λ)
        λnextU += λV(r, params, λ_next)

        for j in i+1:N
            #compute radius and distance vector between particles i and j
            dq .= q[i] .- q[j]
            dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
            r  = sqrt(sum(abs2, dq))

            #compute gradient and contribution to V
            Fij = dVdr(r, params)*dq/r
            dHdq[i] .+= Fij
            dHdq[j] .-= Fij
        end
    end

    return λU, λnextU
end


"""
q-gradient and λU(q)|all λ of the softened particle insertion Hamiltonian, H = T(p) + U(q) + λU(q)
Assumes q is given as Vector{MVector{d, R}} with length(q) = N
Periodic boundary conditions with the standard cutoff L/2 are imposed.
Input functions V, dVdr must have signature r, LJ_params
Input functions λV, dλVdr and dλVdλ must have signature r, LJ_params, λ
"""
function allU_and_dHpi_dq!(dHdq::Vector{MVector{d,Float64}}, q::Vector{MVector{d, Float64}}, params::LJ_params, λ::Real, λ_schedule::Vector{Float64}, L::SVector{d, Float64},
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function) where {d}

    N = length(q)
    @assert N == length(dHdq) "q and dHdq must have the same number of particles"
    N = length(q)

    T = eltype(dHdq[1])
    dq = MVector{d, T}(undef)
    @inbounds for i in 1:N
        dHdq[i] .= zero(T)
    end
    λU = zeros(length(λ_schedule))

    @inbounds for i in 2:N
        #compute derivative of interaction potential wrt the inserted particle, which is labeled by 1
        dq .= q[1] .- q[i]
        dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
        r  = sqrt(sum(abs2, dq))
        F1i = dλVdr(r, params, λ)*dq/r
        dHdq[1] .+= F1i
        dHdq[i] .-= F1i
        λU .+= λV.(Ref(r), Ref(params), λ_schedule)

        for j in i+1:N
            #compute radius and distance vector between particles i and j
            dq .= q[i] .- q[j]
            dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
            r  = sqrt(sum(abs2, dq))

            #compute gradient and contribution to V
            Fij = dVdr(r, params)*dq/r
            dHdq[i] .+= Fij
            dHdq[j] .-= Fij
        end
    end

    return λU
end


"""
p-gradient of any Hamiltonian with the structure H = T(p) + U(q)
"""
function dHpi_dp!(∇_p_Hpi ::Vector{MVector{d, Float64}}, p::Vector{SVector{d, Float64}}, m::Real) where {d}
    @inbounds for i in 1:length{∇_p_Hpi}
        @. ∇_p_Hpi = p[i]/m        
    end
    return nothing
end


"""
Particle instertion Hamiltonian. Versatile formulation that can be integrated with ForwardDiff, but this is not efficient. 
Inserted particle is fixed at the origin during the process
Periodic boundary conditions with the standard cutoff L/2 are imposed.
Input functions V, dVdr must have signature r, LJ_params
Input functions λV, dλVdr and dλVdλ must have signature r, LJ_params, λ
"""
function Hpi(p::Vector{MVector{d,R}}, q::Vector{MVector{d,R}}, m::Real, params::LJ_params, λ::Real, L::SVector{d,R},
    V::Function, λV::Function) where {d,R}

    N = length(q)
    T = typeof(sqrt(sum(abs2,q[1])))

    U = 0
    K = 0
    dq = MVector{d, T}(undef)

    K += sum(abs2, p[1])

    @inbounds for i in 2:N
        dq .= q[1] .- q[i]
        dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
        r  = sqrt(sum(abs2, dq))
        U += λV(r, params, λ)
        K += sum(abs2, p[i])

        for j in i+1:N
            #compute radius and distance vector between particles i and j
            dq .= q[i] .- q[j]
            dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
            r  = sqrt(sum(abs2, dq))
            #compute contribution to U
            U += V(r, params)
        end
    end

    return U+K
end


"""
Particle instertion potential (U).
Inserted particle is fixed at the origin during the process
Periodic boundary conditions with the standard cutoff L/2 are imposed.
Input functions V, dVdr must have signature r, LJ_params
Input functions λV, dλVdr and dλVdλ must have signature r, LJ_params, λ
"""
function Upi(q::Vector{MVector{d,R}}, params::LJ_params, λ::Real, L::SVector{d,R},
    V::Function, λV::Function) where {d,R}

    N = length(q)
    T = typeof(sqrt(sum(abs2,q[1])))

    U = 0
    dq = MVector{d, T}(undef)

    @inbounds for i in 2:N
        dq .= q[1] .- q[i]
        dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
        r  = sqrt(sum(abs2, dq))
        U += λV(r, params, λ)

        for j in i+1:N
            #compute radius and distance vector between particles i and j
            dq .= q[i] .- q[j]
            dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
            r  = sqrt(sum(abs2, dq))
            #compute contribution to U
            U += V(r, params)
        end
    end

    return U   
end


"""
λ-dependent part of the particle instertion potential (U). 
Inserted particle is fixed at the origin during the process
Periodic boundary conditions with the standard cutoff L/2 are imposed.
Input functions V, dVdr must have signature r, LJ_params
Input functions λV, dλVdr and dλVdλ must have signature r, LJ_params, λ
"""
function λUpi(q::Vector{MVector{d,R}}, params::LJ_params, λ::Real, L::SVector{d,R},
    V::Function, λV::Function) where {d,R}

    N = length(q)
    T = typeof(sqrt(sum(abs2,q[1])))

    U = 0
    dq = MVector{d, T}(undef)

    @inbounds for i in 2:N
        #compute derivative of interaction potential wrt the inserted particle, which is labeled by 1
        dq .= q[1] .- q[i]
        dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
        r  = sqrt(sum(abs2, dq))
        U += λV(r, params, λ)
    end

    return U   
end


"""
λ-gradient of the softened particle insertion Hamiltonian, H = T(p) + U(q)
Assumes q is given as Vector{MVector{d, R}} with length(q) = N
Periodic boundary conditions with the standard cutoff L/2 are imposed.
Input functions V, dVdr must have signature r, LJ_params
Input functions λV, dλVdr and dλVdλ must have signature r, LJ_params, λ
"""
function dHpi_dλ(q::Vector{MVector{d, Float64}}, params::LJ_params, λ::Real, dλVdλ::Function) where {d}

    N = length(q)
    dHdλ = 0

    @inbounds for i in 2:N
        dq .= q[1] .- q[i]
        dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
        r  = sqrt(sum(abs2, dq))
        dHdλ += dλVdλ(r, params, λ)
    end

    return dHdλ
end




















######################
# DEPRECATED METHODS #
######################

"""
Compute pairwise distances from a matrix R^{d x N} in which each column encodes the position of a particle 
Very memory efficient, but Distance.jl is better performance-wise if only distances are needed
"""
function pairwise_dist(q::Vector{SVector{R}}) where {R<:Real}
    d, N = length(q[1]), length(q)
    T = typeof(norm(@view q[1]))
    D = Vector{T}(undef, Int(N*(N-1)/2))
    dq = MVector{d, T}(undef)

    idx = 1
    @inbounds for i in 1:N
        for j in i+1:N
            dq .= q[i] .- q[j]
            D[idx] = sqrt(sum(abs2, dq))
            idx += 1
        end
    end
    return D 
end


"""
q-gradient and λ-gradient of the softened particle insertion Hamiltonian.
Assumes q is given as Vector{SVector{d, R}} with length(q) = N
Periodic boundary conditions with the standard cutoff L/2 are imposed.
Input functions V, dVdr must have signature r, LJ_params
Input functions λV, dλVdr and dλVdλ must have signature r, LJ_params, λ
"""
function dHpi_dq_and_dλ(q::Vector{SVector{d,R}}, params::LJ_params, λ::R, L::SVector{d,R},
    dVdr::Function, dVdλ::Function, dλVdr::Function, dλVdλ::Function) where {d,R<:Real}

    N = length(q)
    T = typeof(sqrt(sum(abs2,q[1])))

    dHdq = [@MVector zeros(T, d) for _ in 1:N]
    dHdλ = 0
    dq = MVector{d, T}(undef)

    @inbounds for i in 1:N-1
        #compute derivative of interaction potential with the inserted particle, fixed at the origin
        #and derivative of interaction potential wrt λ
        ri = sqrt(sum(abs2, q[i])) #no need to worry about minimum image convention here
        dHdq[i] .+= dλVdr(ri, params, λ)*q[i]/r0
        dHdλ += dλVdλ(ri, params, λ)

        for j in i+1:N
            #compute radius and distance vector between particles i and j
            dq .= q[i] .- q[j]
            dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
            r  = sqrt(sum(abs2, dq))

            #compute gradient
            Fij = dVdr(r, params)*dq/r
            dHdq[i] .+= Fij
            dHdq[j] .-= Fij
        end
    end

    rN = sqrt(sum(abs2, q[N]))
    dHdq[N] .+= dλVdr(rN, params, λ)*q[N]/rN
    dHdλ += dλVdλ(rN, params, λ)

    return dHdq, dHdλ
    
end


"""
q-gradient and U(q) of the softened particle insertion Hamiltonian, H = T(p) + U(q)
Assumes q is given as Vector{SVector{d, R}} with length(q) = N
Periodic boundary conditions with the standard cutoff L/2 are imposed.
Input functions V, dVdr must have signature r, LJ_params
Input functions λV, dλVdr and dλVdλ must have signature r, LJ_params, λ
"""
function U_and_dHpi_dq(q::Vector{SVector{d,R}}, params::LJ_params, λ::R, L::SVector{d,R},
    V::Function, dVdr::Function, λV::Function, dλVdr::Function) where {d,R<:Real}

    N = length(q)
    T = typeof(sqrt(sum(abs2,q[1])))

    dHdq = [@MVector zeros(T, d) for _ in 1:N]
    U = 0
    dq = MVector{d, T}(undef)

    @inbounds for i in 1:N-1
        #compute derivative of interaction potential with the inserted particle, fixed at the origin
        #and contribution of the λ-term to V
        ri = sqrt(sum(abs2, q[i])) #no need to worry about minimum image convention here
        dHdq[i] .+= dλVdr(ri, params, λ)*q[i]/ri
        U += λV(ri, params, λ)

        for j in i+1:N
            #compute radius and distance vector between particles i and j
            dq .= q[i] .- q[j]
            dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
            r  = sqrt(sum(abs2, dq))

            #compute gradient and contribution to V
            Fij = dVdr(r, params)*dq/r
            dHdq[i] .+= Fij
            dHdq[j] .-= Fij
            U += V(r, params)
        end
    end

    rN = sqrt(sum(abs2, q[N]))
    dHdq[N] .+= dλVdr(rN, params, λ)*q[N]/rN
    U += λV(rN, params, λ)

    return dHdq, U
    
end


"""
Compute pairwise distances from a matrix R^{d x N} in which each column encodes the position of a particle,
 given a certain metric (using Distances.jl)
"""
function pairwise_dist_perf(metric::Union{Metric, PreMetric}, q::Vector{SVector{d, R}}) where {d, R <: Real}
    N = length(q)
    T = eltype(metric(q[1], q[2]))
    D = Vector{T}(undef, Int(N*(N-1)/2)) 
    idx = 1
    @inbounds for i in 1:N-1
        for j in i+1:N
            D[idx] = metric(q[i], q[j])
            idx += 1
        end
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


