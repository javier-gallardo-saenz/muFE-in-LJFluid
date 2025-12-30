using Random

"""
Parameters for Lennard-Jones potential
"""
struct LJ_params
    σ::Float64
    ϵ::Float64
end


"""
Boltzmann Constant in J/K
"""
kB = 1.3806449e-23


"""
Data storage for FEP computation
"""
mutable struct fep_data
    λ1_U::Vector{Float64}
    λ2_U::Vector{Float64}
end

"""
Data storage for μ_ext analysis with BAR and MBAR
"""
mutable struct mu_data_MBAR
    λ_U::Matrix{Float64}
end


#########################################
#
#          UNIT CONVERSION
#
#########################################

"""
Given a number of particles and a desired density, compute the required box length to achieve it
If we use real rho then we get real length
If we use reduced rho then we get reduced length
"""
function box_length(N, rho)
    return (N/rho)^(1/3)
end


"""
Convert reduced temperature Tstar to real T (Kelvin)
"""
function T_real(Tstar, p_real::LJ_params)
    return Tstar*p_real.ϵ/kB
end

"""
Convert reduced density rho star to real density (m^-3)
"""
function rho_real(rho_star, p_real::LJ_params)
    return rho_star / p_real.σ^3
end


"""
Convert reduced chemical potential mu star to real chemical potential (J)
"""
function mu_real(mu_star, p_real::LJ_params)
    return mu_star * p_real.ϵ
end


"""
Get reduced time constant τ
"""
function reduced_time_unit(m_real, p_real::LJ_params)
    return p.σ*sqrt(m/p.ϵ)
end


"""
Convert reduced time units to real time (s)
"""
function t_real(t_star, m_real, p_real::LJ_params)
    return t_star * reduced_time_unit(m_real, p_real)
end




#########################################
#                                   
#         STATISTICAL TOOLS
#
#########################################
"""
Compute autocorrelation time of the observations in the given vector
"""
function autocorr_time(v::Vector{Float64}, maxlag::Int)
    n = length(v)
    μ = mean(v)
    C0 = sum((v .- μ).^2)/n

    τ = 1
    for k in 1:maxlag
        Ck = sum((v[1:end-k] .- μ).*(v[k+1:end] .- μ))/(n-k)
        τ += 2*Ck/C0
    end
    return τ 
end

"""
Compute block average and standard deviation of the observations in the given vector
"""
function block_average(v::Vector{Float64}, block_length::Int)
    n = length(v)
    #get quotient of n/block_length
    #the remainder indicates the number of samples we won't be able to use in our analysis
    num_blocks = div(n, block_length) 
    #put the original vector in blocks, obtaining matrix of size (blocksize, nb)
    blocks = reshape(v[1:num_blocks*block_length], block_length, num_blocks) 
    μ_per_block = vec(mean(blocks, dims=1))
    μ = mean(μ_per_block)
    σ = std(μ_per_block) / sqrt(num_blocks)
    return μ, σ
end


"""
Bootstrap the observations in the given vector
f is the function applied to the samples to obtain our final estimate (scalar)
"""
function bootstrap(v::Vector{Float64}, f::Function, n_boot::Int)
    n = length(v)
    bootstrapped_obs = zeros(n_boot)
    for i in 1:n_boot
        resample = samples[rand(1:n, n)]
        bootstrapped_obs[i] = f(resample)
    end
    μ = mean(bootstrapped_obs)
    σ = std(bootstrapped_obs)

    return μ, σ
end

