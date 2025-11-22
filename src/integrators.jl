using StaticArrays
using Random

using function_utils
using auxiliars



"""
Langevin thermostat with BAOAB integrator.
Does not compute any averages or store values while it runs
"""
function BAOAB!(p::Vector{<:MVector{d,<:Real}}, q::Vector{<:MVector{d,<:Real}}, m::Real, params::LJ_params, λ::Real, L::SVector{d,R},
    dVdr::Function, dλVdr::Function,
      γ::Real, kbT::Real, δt::Real, steps::Int, write_every::Int, save_as::String = "BAOAB_sim",
      load_file::String = "", load_start::Bool = false) where {d::Int, R<:Real}

    if load_start
        step₀, m₀, γ₀, kbT₀, λ₀, δt₀, L₀, params₀, functions₀, p₀, q₀ = load_simulation(load_file)

        step = step₀
        m = m₀
        γ = γ₀
        kbT = kbT₀
        λ = λ₀
        δt = δt₀
        L = L₀
        params = params₀
        p .= p₀
        q .= q₀
        
        dVdr = resolve_function(functions₀.dVdr)
        dλVdr = resolve_function(funcions₀.dλVdr)

        initial_step = step₀ + 1

    else
        initial_step = 1
    end
    
    function_names = H_eval_fun_names()
    H_eval_fun_names.dVdr = string(nameof(dVdr))
    H_eval_fun_names.dλVdr = string(nameof(dλVdr))  

    N = length(q)
    T = typeof(sqrt(sum(abs2,q[1])))

    dHdq = [@MVector zeros(T, d) for _ in 1:N] 
    ξ = @MVector zeros(T,d)

    #Langevin Parameters
    c = exp(-γ*δt)
    σ = sqrt(m*kbT*(1-c^2))

    @inbounds for step in initial_step:steps
        #step B: p -> p + 0.5*δt*F(q)
        dHpi_dq!(dHdq, q, params, λ, L, dVdr, dλVdr)
        @inbounds for i in 1:N
            p[i] .-= 0.5*δt*dHdq[i]
        end

        #step A: q -> q + 0.5*δt*dH/dp 
        @inbounds for i in 1:N
            q[i] .+= 0.5*δt*p[i]/m
        end

        #step O: Langevin Thermostatting
        @inbounds for i in 1:N
            @inbounds for j in 1:d
                ξ[i] = randn()
            end
            p[i] .= c.*p[i] .+ σ.*ξ
        end

        #step A: q -> q + 0.5*δt*dH/dp
        @inbounds for i in 1:N
            q[i] .+= 0.5*δt*p[i]/m
        end

        #step B: p -> p + 0.5*δt*F(q)
        dHpi_dq!(dHdq, q, params, λ, L, dVdr, dλVdr)
        @inbounds for i in 1:N
            p[i] .-= 0.5*δt*dHdq[i]
        end

        if steps%write_every == 0
            save_simulation_state(save_as, p, q, m, params, λ, L, function_names, γ, kbT, δt, step)
        end

    end
end


