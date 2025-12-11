"""
Data storage for μ_ext analysis
"""
struct mu_data
    λ1_U::Vector{Float64}
    λ2_U::Vector{Float64}
end


struct mu_data_MBAR
    λ_U::Vector{Vector{Float64}}
end