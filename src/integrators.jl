using StaticArrays
using Random

using function_utils
using auxiliars



"""
Langevin thermostat with BAOAB integrator.
Does not compute any averages or store values while it runs
"""
function BAOAB!(p::Vector{MVector{d,Float64}}, q::Vector{MVector{d,Float64}}, m::Real, params::LJ_params, λ::Real, L::SVector{d,Float64},
    dVdr::Function, dλVdr::Function,
    γ::Real, kbT::Real, δt::Real, steps::Int, write_every::Int, save_as::String = "BAOAB_sim",
    load_file::String = "", load_start::Bool = false) where {d}

    if load_start
        step₀, m₀, γ₀, kbT₀, λ₀, δt₀, L₀, params₀, functions₀, p₀, q₀, _ = load_simulation(load_file)

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

    @inbounds for step in initial_step:initial_step+steps
        #step B: p -> p + 0.5*δt*F(q)
        dHpi_dq!(dHdq, q, params, λ, L, dVdr, dλVdr)
        @inbounds for i in 1:N
            p[i] .-= 0.5*δt*dHdq[i]
        end

        #step A: q -> q + 0.5*δt*dH/dp (+pbcs)
        @inbounds for i in 1:N
            q[i] .+= 0.5*δt*p[i]/m .- L .* round.(dq./ L)
        end

        #step O: Langevin Thermostatting
        @inbounds for i in 1:N
            @inbounds for j in 1:d
                ξ[i] = randn()
            end
            p[i] .= c.*p[i] .+ σ.*ξ
        end

        #step A: q -> q + 0.5*δt*dH/dp (+ pbcs)
        @inbounds for i in 1:N
            q[i] .+= 0.5*δt*p[i]/m .- L .* round.(dq./ L)
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

    return expectancies
end



"""
Langevin thermostat with BAOAB integrator.
Computes the expectation value of ∂H/∂λ on the go
Compatible with restart
"""
function BAOAB_ti!(p::Vector{MVector{d,Float64}}, q::Vector{MVector{d,Float64}}, m::Real, params::LJ_params, λ::Real, L::SVector{d,Float64}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real, steps::Int, write_every::Int, save_as::String = "BAOAB_E_sim", 
    load_file::Union{nothing, String} = nothing) where {d}

    N = length(q)
    T = typeof(sqrt(sum(abs2,q[1])))

    if load_file !== nothing
        step₀, m₀, γ₀, kbT₀, λ₀, δt₀, L₀, params₀, functions₀, p₀, q₀, expectancy₀ = load_simulation(load_file)

        m = m₀
        γ = γ₀
        kbT = kbT₀
        λ = λ₀
        δt = δt₀
        L = L₀
        params = params₀
        p .= p₀
        q .= q₀
        
        if functions₀.V !== ""
            V = resolve_function(functions₀.V)
        end

        if functions₀.dVdr !== ""
        dVdr = resolve_function(functions₀.dVdr)
        end
        
        if functions₀.λV !== ""
            λV = resolve_function(functions₀.λV)
        end

        if functions₀.dλVdr !== ""
            dλVdr = resolve_function(funcions₀.dλVdr)
        end

        if functions₀.dλVdλ !== ""
            dλVdλ = resolve_function(funcions₀.dλVdλ)
        end

        if expectancy₀ !== nothing
            expectancy = step₀*expectancies₀
        else
            expectancy = 0
        end

        initial_step = step₀ + 1

    else
        initial_step = 1
        #Let expectancy =  <∂H/∂λ>
        expectancy = 0

    end
    
    #store the names of the specific potential functions that are being used in this simulation
    function_names = H_eval_fun_names()
    H_eval_fun_names.dVdr = string(nameof(dVdr))
    H_eval_fun_names.dλVdr = string(nameof(dλVdr))  
    H_eval_fun_names.V = string(nameof(V))
    H_eval_fun_names.λV = string(nameof(λV))
    H_eval_fun_names.dλVdλ = string(nameof(dλVdλ))


    dHdq = [@MVector zeros(T, d) for _ in 1:N] 
    ξ = @MVector zeros(T,d)


    #Langevin Parameters
    c = exp(-γ*δt)
    σ = sqrt(m*kbT*(1-c^2))

    @inbounds for step in initial_step:initial_step+steps
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
        ∂H∂λ = dHpi_dq_and_dλ!(dHdq, q, params, λ, L, V, dVdr, λV, dλVdr, dλVdλ)
        @inbounds for i in 1:N
            p[i] .-= 0.5*δt*dHdq[i]
        end

        #update expectancy
        expectancy += ∂H∂λ

        if steps%write_every == 0
            save_simulation_state(save_as, p, q, m, params, λ, L, function_names, γ, kbT, δt, step, expectancy/step)
        end

    end

    expectancies /= initial_step+steps

    return nothing  
end



"""
Langevin thermostat with BAOAB integrator.
Collects data for FEP and BAR
Only half compatible with restart, collected data for FEP and BAR needs to be independently processsed
"""
function BAOAB_twoλ!(p::Vector{MVector{d,Float64}}, q::Vector{MVector{d,Float64}}, m::Real, params::LJ_params, λ::Real, λ_next::Real, L::SVector{d,Float64}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real, steps::Int, write_every::Int, save_as::String = "BAOAB_E_sim", 
    load_file::Union{nothing, String} = nothing) where {d}

    N = length(q)
    T = typeof(sqrt(sum(abs2,q[1])))

    if load_file !== nothing
        step₀, m₀, γ₀, kbT₀, λ₀, δt₀, L₀, params₀, functions₀, p₀, q₀, _ = load_simulation(load_file)

        m = m₀
        γ = γ₀
        kbT = kbT₀
        λ = λ₀
        δt = δt₀
        L = L₀
        params = params₀
        p .= p₀
        q .= q₀
        
        if functions₀.V !== ""
            V = resolve_function(functions₀.V)
        end

        if functions₀.dVdr !== ""
        dVdr = resolve_function(functions₀.dVdr)
        end
        
        if functions₀.λV !== ""
            λV = resolve_function(functions₀.λV)
        end

        if functions₀.dλVdr !== ""
            dλVdr = resolve_function(funcions₀.dλVdr)
        end

        if functions₀.dλVdλ !== ""
            dλVdλ = resolve_function(funcions₀.dλVdλ)
        end

        initial_step = step₀ + 1

    else
        initial_step = 1
    end
    
    #store the names of the specific potential functions that are being used in this simulation
    function_names = H_eval_fun_names()
    H_eval_fun_names.dVdr = string(nameof(dVdr))
    H_eval_fun_names.dλVdr = string(nameof(dλVdr))  
    H_eval_fun_names.V = string(nameof(V))
    H_eval_fun_names.λV = string(nameof(λV))
    H_eval_fun_names.dλVdλ = string(nameof(dλVdλ))

    #initiate vectors for updates 
    dHdq = [@MVector zeros(T, d) for _ in 1:N] 
    ξ = @MVector zeros(T,d)

    #Langevin Parameters
    c = exp(-γ*δt)
    σ = sqrt(m*kbT*(1-c^2))

    @inbounds for step in initial_step:initial_step+steps
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
        λU, λUnext = twoU_and_dHpi_dq!(dHdq, q, params, λ, λ_next, L, V, dVdr, λV, dλVdr, dλVdλ)
        K = 0
        @inbounds for i in 1:N
            p[i] .-= 0.5*δt*dHdq[i]
        end


        if steps%write_every == 0
            save_simulation_state(save_as, p, q, m, params, λ, L, function_names, γ, kbT, δt, step)
        end

    end

    expectancies /= initial_step+steps

    return nothing  
end


