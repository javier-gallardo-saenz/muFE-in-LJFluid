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
