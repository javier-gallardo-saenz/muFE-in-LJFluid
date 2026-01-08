using StaticArrays
using Random
using ProgressBars

"""
1 step of the Langevin thermostat with BAOAB integrator
For efficient multistep integration, see next function.
"""
function onestep_BAOAB!(p::Vector{MVector{d,Float64}}, q::Vector{MVector{d,Float64}}, m::Real, params::LJ_params, λ::Real, L::SVector{d,Float64},
    dVdr::Function, dλVdr::Function, γ::Real, kbT::Real, δt::Real) where {d}

    N = length(q)
    T = typeof(sqrt(sum(abs2,q[1])))

    dHdq = [@MVector zeros(T, d) for _ in 1:N] 
    ξ = @MVector zeros(T,d)

    #Langevin Parameters
    c = exp(-γ*δt)
    σ = sqrt(m*kbT*(1-c^2))

    #step B: p -> p + 0.5*δt*F(q)
    dHpi_dq!(dHdq, q, params, λ, L, dVdr, dλVdr)
    @inbounds for i in 1:N
        p[i] .-= 0.5*δt*dHdq[i]
    end

    #step A: q -> q + 0.5*δt*dH/dp (+pbcs)
    @inbounds for i in 1:N
        q[i] .+= 0.5*δt*p[i]/m
        q[i] .-= L .* round.(q[i]./ L)
    end

    #step O: Langevin Thermostatting
    @inbounds for i in 1:N
        @inbounds for j in 1:d
            ξ[j] = randn()
        end
        p[i] .= c.*p[i] .+ σ.*ξ
    end

    #step A: q -> q + 0.5*δt*dH/dp (+ pbcs)
    @inbounds for i in 1:N
        q[i] .+= 0.5*δt*p[i]/m 
        q[i] .-= L .* round.(q[i]./ L)
    end

    #step B: p -> p + 0.5*δt*F(q)
    dHpi_dq!(dHdq, q, params, λ, L, dVdr, dλVdr)
    @inbounds for i in 1:N
        p[i] .-= 0.5*δt*dHdq[i]
    end

    return nothing

end


"""
Langevin thermostat with BAOAB integrator.
Does not compute any averages or store values while it runs
It supports loading
"""
function BAOAB_load_and_save!(p::Vector{MVector{d,Float64}}, q::Vector{MVector{d,Float64}}, m::Real, params::LJ_params, λ::Real, L::SVector{d,Float64},
    dVdr::Function, dλVdr::Function,
    γ::Real, kbT::Real, δt::Real, steps::Int, write_every::Int, save_as::String = "BAOAB_sim",
    load_file::Union{Nothing, String} = nothing, load_start::Bool = false) where {d}

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
    
    function_names = H_eval_fun_names("","","","","")
    function_names.dVdr = string(nameof(dVdr))
    function_names.dλVdr = string(nameof(dλVdr))  

    N = length(q)
    T = typeof(sqrt(sum(abs2,q[1])))

    dHdq = [@MVector zeros(T, d) for _ in 1:N] 
    ξ = @MVector zeros(T,d)

    #Langevin Parameters
    c = exp(-γ*δt)
    σ = sqrt(m*kbT*(1-c^2))

    println("Running BAOAB Langevin")
    @inbounds for step in ProgressBar(initial_step:initial_step+steps-1)
        #step B: p -> p + 0.5*δt*F(q)
        dHpi_dq!(dHdq, q, params, λ, L, dVdr, dλVdr)
        @inbounds for i in 1:N
            p[i] .-= 0.5*δt*dHdq[i]
        end

        #step A: q -> q + 0.5*δt*dH/dp (+pbcs)
        @inbounds for i in 1:N
            q[i] .+= 0.5*δt*p[i]/m
            q[i] .-= L .* round.(q[i]./ L)
        end

        #step O: Langevin Thermostatting
        @inbounds for i in 1:N
            @inbounds for j in 1:d
                ξ[j] = randn()
            end
            p[i] .= c.*p[i] .+ σ.*ξ
        end

        #step A: q -> q + 0.5*δt*dH/dp (+ pbcs)
        @inbounds for i in 1:N
            q[i] .+= 0.5*δt*p[i]/m 
            q[i] .-= L .* round.(q[i]./ L)
        end

        #step B: p -> p + 0.5*δt*F(q)
        dHpi_dq!(dHdq, q, params, λ, L, dVdr, dλVdr)
        @inbounds for i in 1:N
            p[i] .-= 0.5*δt*dHdq[i]
        end

        if step%write_every == 0
            save_simulation_state(save_as, p, q, m, params, λ, L, function_names, γ, kbT, δt, step)
        end

    end

    return nothing
end

"""
Langevin thermostat with BAOAB integrator.
Does not compute any averages or store values while it runs
"""
function BAOAB!(p::Vector{MVector{d,Float64}}, q::Vector{MVector{d,Float64}}, m::Real, params::LJ_params, λ::Real, L::SVector{d,Float64},
    dVdr::Function, dλVdr::Function,
    γ::Real, kbT::Real, δt::Real, steps::Int) where {d}

    N = length(q)
    T = typeof(sqrt(sum(abs2,q[1])))

    dHdq = [@MVector zeros(T, d) for _ in 1:N] 
    ξ = @MVector zeros(T,d)

    #Langevin Parameters
    c = exp(-γ*δt)
    σ = sqrt(m*kbT*(1-c^2))

    println("Running BAOAB Langevin")
    @inbounds for step in ProgressBar(1:steps)
        #step B: p -> p + 0.5*δt*F(q)
        dHpi_dq!(dHdq, q, params, λ, L, dVdr, dλVdr)
        @inbounds for i in 1:N
            p[i] .-= 0.5*δt*dHdq[i]
        end

        #step A: q -> q + 0.5*δt*dH/dp (+pbcs)
        @inbounds for i in 1:N
            q[i] .+= 0.5*δt*p[i]/m
            q[i] .-= L .* round.(q[i]./ L)
        end

        #step O: Langevin Thermostatting
        @inbounds for i in 1:N
            @inbounds for j in 1:d
                ξ[j] = randn()
            end
            p[i] .= c.*p[i] .+ σ.*ξ
        end

        #step A: q -> q + 0.5*δt*dH/dp (+ pbcs)
        @inbounds for i in 1:N
            q[i] .+= 0.5*δt*p[i]/m 
            q[i] .-= L .* round.(q[i]./ L)
        end

        #step B: p -> p + 0.5*δt*F(q)
        dHpi_dq!(dHdq, q, params, λ, L, dVdr, dλVdr)
        @inbounds for i in 1:N
            p[i] .-= 0.5*δt*dHdq[i]
        end

    end

    return nothing
end


"""
Langevin thermostat with BAOAB integrator.
Computes the expectation value of ∂H/∂λ on the go
Compatible with restart
"""
function BAOAB_ti!(data_mu::mu_data, p::Vector{MVector{d,Float64}}, q::Vector{MVector{d,Float64}}, m::Real, params::LJ_params, λ::Real, L::SVector{d,Float64}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real, steps::Int, sample_every::Int) where {d}

    N = length(q)
    T = typeof(sqrt(sum(abs2,q[1])))

    dHdq = [@MVector zeros(T, d) for _ in 1:N] 
    ξ = @MVector zeros(T,d)

    #Langevin Parameters
    c = exp(-γ*δt)
    σ = sqrt(m*kbT*(1-c^2))

    @inbounds for step in ProgressBar(1:steps)
        #step B: p -> p + 0.5*δt*F(q)
        dHpi_dq!(dHdq, q, params, λ, L, dVdr, dλVdr)
        @inbounds for i in 1:N
            p[i] .-= 0.5*δt*dHdq[i]
        end

        #step A: q -> q + 0.5*δt*dH/dp 
        @inbounds for i in 1:N
            q[i] .+= 0.5*δt*p[i]/m
            q[i] .-= L .* round.(q[i]./ L)
        end

        #step O: Langevin Thermostatting
        @inbounds for i in 1:N
            @inbounds for j in 1:d
                ξ[j] = randn()
            end
            p[i] .= c.*p[i] .+ σ.*ξ
        end

        #step A: q -> q + 0.5*δt*dH/dp
        @inbounds for i in 1:N
            q[i] .+= 0.5*δt*p[i]/m
            q[i] .-= L .* round.(q[i]./ L)
        end

        #step B: p -> p + 0.5*δt*F(q)
        if step%sample_every == 0
            ∂H∂λ = dHpi_dq_and_dλ!(dHdq, q, params, λ, L, V, dVdr, λV, dλVdr, dλVdλ)
            push!(data_mu.w, ∂H∂λ)
        else
            dHpi_dq!(dHdq, q, params, λ, L, dVdr, dλVdr)
        end
        @inbounds for i in 1:N
            p[i] .-= 0.5*δt*dHdq[i]
        end
    end

    return nothing
end



"""
Langevin thermostat with BAOAB integrator.
Collects data for FEP and BAR
Only half compatible with restart, collected data for FEP and BAR needs to be independently processsed
"""
function BAOAB_twoλ!(ΔU::Vector{Float64},
    p::Vector{MVector{d,Float64}}, q::Vector{MVector{d,Float64}}, m::Real, params::LJ_params, λ::Real, λ_next::Real, L::SVector{d,Float64}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real, steps::Int, sample_every::Int) where {d}

    if length(ΔU) != div(steps,sample_every)
        throw(DimensionMismatch("The length of the auxiliary vector ($(length(ΔU))) should match the number of samples that will be taken ($(div(steps,sample_every)))"))
    end

    N = length(q)
    T = typeof(sqrt(sum(abs2,q[1])))

    initial_step = 1

    #initiate vectors for updates 
    dHdq = [@MVector zeros(T, d) for _ in 1:N] 
    ξ = @MVector zeros(T,d)

    #Langevin Parameters
    c = exp(-γ*δt)
    σ = sqrt(m*kbT*(1-c^2))

    @inbounds for step in ProgressBar(initial_step:initial_step+steps-1)
        #step B: p -> p + 0.5*δt*F(q)
        dHpi_dq!(dHdq, q, params, λ, L, dVdr, dλVdr)
        @inbounds for i in 1:N
            p[i] .-= 0.5*δt*dHdq[i]
        end

        #step A: q -> q + 0.5*δt*dH/dp 
        @inbounds for i in 1:N
            q[i] .+= 0.5*δt*p[i]/m
            q[i] .-= L .* round.(q[i]./ L)
        end

        #step O: Langevin Thermostatting
        @inbounds for i in 1:N
            @inbounds for j in 1:d
                ξ[j] = randn()
            end
            p[i] .= c.*p[i] .+ σ.*ξ
        end

        #step A: q -> q + 0.5*δt*dH/dp
        @inbounds for i in 1:N
            q[i] .+= 0.5*δt*p[i]/m
            q[i] .-= L .* round.(q[i]./ L)
        end

        #step B: p -> p + 0.5*δt*F(q)
        if step % sample_every == 0
            index = div(step, sample_every)
            U1, U2 = twoU_and_dHpi_dq!(dHdq, q, params, λ, λ_next, L, V, dVdr, λV, dλVdr, dλVdλ)
            ΔU[index] = U2 - U1
        else 
            dHpi_dq!(dHdq, q, params, λ, L, dVdr, dλVdr)
        end

        @inbounds for i in 1:N
            p[i] .-= 0.5*δt*dHdq[i]
        end
        
    end

    return nothing  
end


"""
Langevin thermostat with BAOAB integrator.
Collects data for MBAR
Only half compatible with restart, collected data for FEP and BAR needs to be independently processsed
"""
function BAOAB_allλ!(data_mu::mu_data_MBAR,
    p::Vector{MVector{d,Float64}}, q::Vector{MVector{d,Float64}}, m::Real, params::LJ_params, λ::Real, λ_schedule::Vector{Float64}, L::SVector{d,Float64}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real, steps::Int, sample_every::Int) where {d}

    S, _ = size(data_mu.λ_U)
    exp_S = div(steps, sample_every)
    if S != exp_S
        throw(DimensionMismatch("The provided struct for MBAR data storage first dimension ($(S)) must match the number of samples that the integrator will collect ($(exp_S))"))
    end
    current_sample = 1

    N = length(q)
    T = typeof(sqrt(sum(abs2,q[1])))

    #initiate vectors for updates 
    dHdq = [@MVector zeros(T, d) for _ in 1:N] 
    ξ = @MVector zeros(T,d)

    #Langevin Parameters
    c = exp(-γ*δt)
    σ = sqrt(m*kbT*(1-c^2))

    @inbounds for step in ProgressBar(1:steps)
        #step B: p -> p + 0.5*δt*F(q)
        dHpi_dq!(dHdq, q, params, λ, L, dVdr, dλVdr)
        @inbounds for i in 1:N
            p[i] .-= 0.5*δt*dHdq[i]
        end

        #step A: q -> q + 0.5*δt*dH/dp 
        @inbounds for i in 1:N
            q[i] .+= 0.5*δt*p[i]/m 
            q[i] .-= L .* round.(q[i]./ L)
        end

        #step O: Langevin Thermostatting
        @inbounds for i in 1:N
            @inbounds for j in 1:d
                ξ[j] = randn()
            end
            p[i] .= c.*p[i] .+ σ.*ξ
        end

        #step A: q -> q + 0.5*δt*dH/dp
        @inbounds for i in 1:N
            q[i] .+= 0.5*δt*p[i]/m 
            q[i] .-= L .* round.(q[i]./ L)
        end

        #step B: p -> p + 0.5*δt*F(q)
        if step % sample_every == 0
            allλU = allU_and_dHpi_dq!(dHdq, q, params, λ, λ_schedule, L, V, dVdr, λV, dλVdr, dλVdλ)
            data_mu.λ_U[current_sample, :] = allλU
            current_sample += 1
        else
            dHpi_dq!(dHdq, q, params, λ, L, dVdr, dλVdr)
        end

        @inbounds for i in 1:N
            p[i] .-= 0.5*δt*dHdq[i]
        end

    end

    return nothing  
end


"""
Langevin thermostat with BAOAB integrator.
Computes radial distribution function
"""
function BAOAB_rdf!(hist::Vector{Float64}, dr::Float64, p::Vector{MVector{d,Float64}}, q::Vector{MVector{d,Float64}}, m::Real, params::LJ_params, λ::Real, L::SVector{d,Float64},
    dVdr::Function, dλVdr::Function,
    γ::Real, kbT::Real, δt::Real, steps::Int, sample_every::Int) where {d}

    N = length(q)
    T = typeof(sqrt(sum(abs2,q[1])))

    dHdq = [@MVector zeros(T, d) for _ in 1:N] 
    ξ = @MVector zeros(T,d)
    dq = @MVector zeros(T,d)

    #Langevin Parameters
    c = exp(-γ*δt)
    σ = sqrt(m*kbT*(1-c^2))

    println("Running BAOAB Langevin")
    @inbounds for step in ProgressBar(1:steps)
        #step B: p -> p + 0.5*δt*F(q)
        dHpi_dq!(dHdq, q, params, λ, L, dVdr, dλVdr)
        @inbounds for i in 1:N
            p[i] .-= 0.5*δt*dHdq[i]
        end

        #step A: q -> q + 0.5*δt*dH/dp (+pbcs)
        @inbounds for i in 1:N
            q[i] .+= 0.5*δt*p[i]/m
            q[i] .-= L .* round.(q[i]./ L)
        end

        #step O: Langevin Thermostatting
        @inbounds for i in 1:N
            @inbounds for j in 1:d
                ξ[j] = randn()
            end
            p[i] .= c.*p[i] .+ σ.*ξ
        end

        #step A: q -> q + 0.5*δt*dH/dp (+ pbcs)
        @inbounds for i in 1:N
            q[i] .+= 0.5*δt*p[i]/m 
            q[i] .-= L .* round.(q[i]./ L)
        end

        #step B: p -> p + 0.5*δt*F(q)
        dHpi_dq!(dHdq, q, params, λ, L, dVdr, dλVdr)
        @inbounds for i in 1:N
            p[i] .-= 0.5*δt*dHdq[i]
        end

        if step%sample_every == 0
            @inbounds for i in 1:N-1
                for j in i+1:N
                    dq .= q[i] .- q[j]
                    dq .= dq .- L .* round.(dq./ L)
                    r  = sqrt(sum(abs2, dq))
                    k = Int(floor(r/dr)) + 1
                    if 1 ≤ k ≤ length(hist)
                        hist[k] += 1
                    end
                end
            end
        end

    end

    return nothing
end


