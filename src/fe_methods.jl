using StaticArrays
using Random
using ProgressBars
using Roots
using Statistics


"""
Widom method
NOT OPTIMIZED
"""
function widom_method(sample_every::Int, n_trial_insertions::Int, n_tries_per_insertion::Int,
    N::Int, m::Real, params::LJ_params, L::SVector{d,R}, 
    V::Function, dVdr::Function, dλVdr::Function,
    γ::Real, kbT::Real, δt::Real, write_every::Int = 10000) where {d, R<:Real}

    #sample initial momentums from the Boltzmann distribution
    p = [MVector{d, R}(randn(d) .* sqrt(m * kbT)) for _ in 1:N]

    #fix intial positions uniformly on the box
    n = ceil(Int, N^(1/d)) #grid size
    grid = range(0,L[1];length=n+1)[1:end-1] 
    coords = collect(Iterators.product(ntuple(_->grid,d)...))

    q = Vector{MVector{d,R}}(undef,N)
    for i in 1:N
        q[i] = MVector{d,R}(Float64.(coords[i]))
    end

    #initial equilibration, longer to get rid of effects of unphysical starting state 
    println("Initial equilibration/relaxation")
    BAOAB!(p, q, m, params, 0.0, L, dVdr, dλVdr, γ, kbT, δt, 15000, write_every)

    avg_exp_minusbetadU = 0.0
    dq = MVector{d, R}(undef)
    dUs = zeros(n_trial_insertions*n_tries_per_insertion)

    for i in ProgressBar(1:n_trial_insertions)
        BAOAB!(p, q, m, params, 0.0, L, dVdr, dλVdr, γ, kbT, δt, sample_every, write_every)
        
        for k in 1:n_tries_per_insertion
            dU = 0
            rand_pos = MVector{d, Float64}(rand(d).*L)
            @inbounds for i in 1:N 
                #compute radius and distance vector between particles i and the inserted one
                dq .= q[i] .- rand_pos
                dq .= dq .- L .* round.(dq./ L) #pbc + minimum image convention
                r  = sqrt(sum(abs2, dq))
                #compute contribution to ΔV
                dU += V(r, params)
                dUs[(i-1)*n_tries_per_insertion + k] = dU
            end
            #handle numerical explosions if insertion was too close to a particle
            if isnan(dU)||isinf(dU)
                aux = 0
            else
                aux = exp(-dU/kbT)
            end  
            avg_exp_minusbetadU += aux   
        end
    end
    avg_exp_minusbetadU /= n_tries_per_insertion*n_trial_insertions
    println(avg_exp_minusbetadU)
    return -kbT*log(avg_exp_minusbetadU), dUs
end



"""
Perform thermodynamic integration (TI)
"""
function thermodynamic_integration(initial_eq_steps::Int, eq_steps::Int, prod_steps::Int, N::Int, m::Real, params::LJ_params, λ_steps::Int, L::SVector{d,R}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real, write_every::Int) where {d, R<:Real}

    #define λ schedule
    λ_schedule = range(0, 1; length = λ_steps+1)

    #sample initial momentums from the Boltzmann distribution
    p = [MVector{d, R}(randn(d) .* sqrt(m * kbT)) for _ in 1:N]

    #fix intial positions uniformly on the box
    n = ceil(Int, N^(1/d)) #grid size
    grid = range(0,L[1];length=n+1)[1:end-1] 
    coords = collect(Iterators.product(ntuple(_->grid,d)...))

    q = Vector{MVector{d,R}}(undef,N)
    for i in 1:N
        q[i] = MVector{d,R}(Float64.(coords[i]))
    end

    #initial equilibration, longer to get rid of effects of unphysical starting state 
    println("Initial equilibration/relaxation")
    BAOAB!(p, q, m, params, λ_schedule[1], L, dVdr, dλVdr, γ, kbT, δt, initial_eq_steps, write_every)

    #vectors to store <∂H/∂λ>
    dHdλ_avgs = zeros(length(λ_schedule))

    println("Running MD for the different λ values")
    for i in ProgressBar(eachindex(λ_schedule))
        #equilibrate
        println("Equilibrating")
        BAOAB!(p, q, m, params, λ_schedule[i], L, dVdr, dλVdr, γ, kbT, δt, eq_steps, write_every)
        #production run
        message = "Production run for λ=$(λ_schedule[i])"
        println(message)
        dHdλ_avgs[i] = BAOAB_ti!(p, q, m, params, λ_schedule[i], L, V, dVdr, λV, dλVdr, dλVdλ,
         γ, kbT, δt, prod_steps, write_every)
        println(dHdλ_avgs[i])
    end
    println(dHdλ_avgs)
    
    
    #for visualization
    μ_ext_ti = zeros(length(λ_schedule))
    for i in 2:(λ_steps + 1)
        δλ = λ_schedule[i] - λ_schedule[i-1]
        μ_ext_ti[i] = μ_ext_ti[i-1] + 0.5*δλ*(dHdλ_avgs[i] + dHdλ_avgs[i-1])
    end
    
    return μ_ext_ti
end


"""
Collect data for FEP
"""
function collect_fep(initial_eq_steps::Int, eq_steps::Int, prod_steps::Int, N::Int, m::Real, params::LJ_params, λ_steps::Int, L::SVector{d,R}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real, write_every::Int) where {d, R<:Real}
    
    data = [fep_data(Vector{Float64}(), Vector{Float64}()) for _ in 1:λ_steps]
    #define λ schedule
    λ_schedule = range(0, 1; length = λ_steps+1)

    #sample initial momentums from the Boltzmann distribution
    p = [MVector{d, R}(randn(d) .* sqrt(m * kbT)) for _ in 1:N]

    #fix intial positions uniformly on the box
    n = ceil(Int, N^(1/d)) #grid size
    grid = range(0,L[1];length=n+1)[1:end-1] 
    coords = collect(Iterators.product(ntuple(_->grid,d)...))

    q = Vector{MVector{d,R}}(undef,N)
    for i in 1:N
        q[i] = MVector{d,R}(Float64.(coords[i]))
    end

    #initial equilibration, longer to get rid of effects of unphysical starting state 
    println("Initial equilibration/relaxation")
    BAOAB!(p, q, m, params, λ_schedule[1], L, dVdr, dλVdr, γ, kbT, δt, initial_eq_steps, write_every)

    for i in ProgressBar(1:λ_steps)
        #equilibrate
        println("Equilibrating")
        BAOAB!(p, q, m, params, λ_schedule[i], L, dVdr, dλVdr, γ, kbT, δt, eq_steps, write_every)
        #production run
        message = "Production run for λ=$(λ_schedule[i])"
        println(message)
        BAOAB_twoλ!(data[i], p, q, m, params, λ_schedule[i], λ_schedule[i+1], L, V, dVdr, λV, dλVdr, dλVdλ,
         γ, kbT, δt, prod_steps)
    end

    return data
end


"""
Auxiliary function for FEP computation
"""
function fep_aux(d::fep_data, kbT::Real)
    aux = mean(exp.((d.λ1_U .- d.λ2_U)./kbT))
    return -kbT*log(aux)
end


"""
Compute the free energy difference using FEP
"""
function fep(initial_eq_steps::Int, eq_steps::Int, prod_steps::Int, N::Int, m::Real, params::LJ_params, λ_steps::Int, L::SVector{d,R}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real, write_every::Int) where {d, R<:Real}

    #generate the necessary data for both FEP and BAR computation of the free energy
    data = collect_fep(initial_eq_steps, eq_steps, prod_steps, N, m, params, λ_steps, L, V, dVdr, λV, dλVdr, dλVdλ, 
    γ, kbT, δt, write_every)

    ΔF_fep = cumsum(fep_aux.(data, kbT))

    return ΔF_fep
end


"""
Collect data for MBAR (used for BAR as well)
"""
function collect_mbar(initial_eq_steps::Int, eq_steps::Int, prod_steps::Int, N::Int, m::Real, params::LJ_params, λ_steps::Int, L::SVector{d,R}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real, write_every::Int) where {d, R<:Real}
    
    data = [mu_data_MBAR(Vector{Vector{Float64}}()) for _ in 1:(λ_steps+1)]
    #define λ schedule
    λ_schedule = range(0, 1; length = λ_steps+1)
    λ_vector = Vector{Float64}(λ_schedule)

    #sample initial momentums from the Boltzmann distribution
    p = [MVector{d, R}(randn(d) .* sqrt(m * kbT)) for _ in 1:N]

    #fix intial positions uniformly on the box
    n = ceil(Int, N^(1/d)) #grid size
    grid = range(0,L[1];length=n+1)[1:end-1] 
    coords = collect(Iterators.product(ntuple(_->grid,d)...))

    q = Vector{MVector{d,R}}(undef,N)
    for i in 1:N
        q[i] = MVector{d,R}(Float64.(coords[i]))
    end

    #initial equilibration, longer to get rid of effects of unphysical starting state 
    println("Initial equilibration/relaxation")
    BAOAB!(p, q, m, params, λ_schedule[1], L, dVdr, dλVdr, γ, kbT, δt, initial_eq_steps, write_every)

    for i in ProgressBar(1:λ_steps+1)
        #equilibrate
        println("Equilibrating")
        BAOAB!(p, q, m, params, λ_schedule[i], L, dVdr, dλVdr, γ, kbT, δt, eq_steps, write_every)
        #production run
        message = "Production run for λ=$(λ_schedule[i])"
        println(message)
        BAOAB_allλ!(data[i], p, q,  m, params, λ_schedule[i], λ_vector, L, V, dVdr, λV, dλVdr, dλVdλ,
         γ, kbT, δt, prod_steps)
    end

    return data
end


"""
Auxiliary function for BAR computation
"""
function BAR_aux(ΔλU_12::Vector{Float64}, ΔλU_21::Vector{Float64}, kbT::Real)
    #define function that needs to be zero
    f(ΔF) = mean(1 ./(1 .+ exp.((ΔλU_12 .- ΔF)./kbT))) - mean(1 ./(1 .+ exp.(-(ΔλU_21 .+ ΔF)./kbT)))
    ΔF_bar = find_zero(f, 0.0) 
    return ΔF_bar
end


"""
Auxiliary function for MBAR computation
"""
function MBAR_aux(ΔλU_12::Vector{Float64})
end


"""
Compute the free energy difference using BAR and MBAR
"""
function BAR(initial_eq_steps::Int, eq_steps::Int, prod_steps::Int, N::Int, m::Real, params::LJ_params, λ_steps::Int, L::SVector{d,R}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real, write_every::Int) where {d, R<:Real}

    #generate the necessary data for both FEP and BAR computation of the free energy
    data = collect_mbar(initial_eq_steps, eq_steps, prod_steps, N, m, params, λ_steps, L, V, dVdr, λV, dλVdr, dλVdλ, 
    γ, kbT, δt, write_every)

    local_ΔF = zeros(λ_steps)

    for i in 1:λ_steps
            ΔλU_12 = data[i].λ_U[:, i+1] .- data[i].λ_U[:, i]
            ΔλU_21 = data[i+1].λ_U[:,i] .- data[i+1].λ_U[:, i+1] 
            local_ΔF[i] = BAR_aux(ΔλU_12, ΔλU_21, kbT)
    end
        
    return cumsum(local_ΔF)       
end



