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
function LJ_pot(r::Float64, params::LJ_params)
    aux = r/params.σ
    return 4*params.ϵ*((aux)^12 - (aux)^6)
end


"""
Softened interaction potential
"""
function soft_pot(r::Float64, params::LJ_params, λ::Float64)
    A = 0.5*(1-λ) + (r/params.σ)^6
    return 4*params.ϵ*(A^(-2) - A^(-1))
end

"""
Particle instertion Hamiltonian
"""
function H_pi(p, q, m, params, λ)
end


"""
Softened particle insertion Hamiltonian
"""
function H_spi(p, q, m, params, λ)
end