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
    γ::Real, kbT::Real, δt::Real) where {d, R<:Real}

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
    BAOAB!(p, q, m, params, 0.0, L, dVdr, dλVdr, γ, kbT, δt, 15000)

    avg_exp_minusbetadU = 0.0
    dq = MVector{d, R}(undef)
    dUs = zeros(n_trial_insertions*n_tries_per_insertion)

    for i in ProgressBar(1:n_trial_insertions)
        BAOAB!(p, q, m, params, 0.0, L, dVdr, dλVdr, γ, kbT, δt, sample_every)
        
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
    γ::Real, kbT::Real, δt::Real, sample_every::Int) where {d, R<:Real}

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
    BAOAB!(p, q, m, params, λ_schedule[1], L, dVdr, dλVdr, γ, kbT, δt, initial_eq_steps)

    #vectors to store <∂H/∂λ>
    dHdλ_avgs = zeros(length(λ_schedule))

    println("Running MD for the different λ values")
    for i in ProgressBar(eachindex(λ_schedule))
        #equilibrate
        println("Equilibrating")
        BAOAB!(p, q, m, params, λ_schedule[i], L, dVdr, dλVdr, γ, kbT, δt, eq_steps)
        #production run
        message = "Production run for λ=$(λ_schedule[i])"
        println(message)
        dHdλ_avgs[i] = BAOAB_ti!(p, q, m, params, λ_schedule[i], L, V, dVdr, λV, dλVdr, dλVdλ,
         γ, kbT, δt, prod_steps, sample_every)
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
    γ::Real, kbT::Real, δt::Real) where {d, R<:Real}
    
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
    BAOAB!(p, q, m, params, λ_schedule[1], L, dVdr, dλVdr, γ, kbT, δt, initial_eq_steps)

    #compute autocorrelation time
    autocorr_comp_steps = 10000
    Us = zeros(autocorr_comp_steps)
    println("Collecting data to estimate the autocorrelation time of the potential energy")
    @inbounds for i in ProgressBar(1:autocorr_comp_steps)
        onestep_BAOAB!(p, q, m, params, λ_schedule[1], L, dVdr, dλVdr, γ, kbT, δt)
        Us[i] = Upi(q, params, 0.0, L, V, λV)
    end
    τ = autocorr_time(Us, 2000)
    println("The autocorrelation time of the potential energy in this steup is approximately $(τ)")
    sample_every = ceil(Int, τ)
    println("Sampling every $(sample_every) steps")

    for i in ProgressBar(1:λ_steps)
        #equilibrate
        println("Equilibrating")
        BAOAB!(p, q, m, params, λ_schedule[i], L, dVdr, dλVdr, γ, kbT, δt, eq_steps)
        #production run
        message = "Production run for λ=$(λ_schedule[i])"
        println(message)
        BAOAB_twoλ!(data[i], p, q, m, params, λ_schedule[i], λ_schedule[i+1], L, V, dVdr, λV, dλVdr, dλVdλ,
         γ, kbT, δt, prod_steps, sample_every)
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
    γ::Real, kbT::Real, δt::Real) where {d, R<:Real}

    #generate the necessary data for both FEP and BAR computation of the free energy
    data = collect_fep(initial_eq_steps, eq_steps, prod_steps, N, m, params, λ_steps, L, V, dVdr, λV, dλVdr, dλVdλ, 
    γ, kbT, δt)

    ΔF_fep = cumsum(fep_aux.(data, kbT))

    println("<λU> for λ_schedule[λ_steps] is:")
    println(mean(data[λ_steps].λ1_U))

    return ΔF_fep
end


"""
Collect data for MBAR (used for BAR as well)
"""
function collect_mbar(initial_eq_steps::Int, eq_steps::Int, prod_steps::Int, N::Int, m::Real, params::LJ_params, λ_steps::Int, L::SVector{d,R}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real) where {d, R<:Real}
    
    data = [mu_data_MBAR(zeros(prod_steps, λ_steps+1)) for _ in 1:(λ_steps+1)]
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
    BAOAB!(p, q, m, params, λ_schedule[1], L, dVdr, dλVdr, γ, kbT, δt, initial_eq_steps)

    #compute autocorrelation time
    autocorr_comp_steps = 10000
    Us = zeros(autocorr_comp_steps)
    println("Collecting data to estimate the autocorrelation time of the potential energy")
    @inbounds for i in ProgressBar(1:autocorr_comp_steps)
        onestep_BAOAB!(p, q, m, params, λ_schedule[1], L, dVdr, dλVdr, γ, kbT, δt)
        Us[i] = Upi(q, params, 0.0, L, V, λV)
    end
    τ = autocorr_time(Us, 2000)
    println("The autocorrelation time of the potential energy in this steup is approximately $(τ)")
    sample_every = ceil(Int, τ)
    println("Sampling every $(sample_every) steps")

    for i in ProgressBar(1:λ_steps+1)
        #equilibrate
        println("Equilibrating")
        BAOAB!(p, q, m, params, λ_schedule[i], L, dVdr, dλVdr, γ, kbT, δt, eq_steps)
        #production run
        message = "Production run for λ=$(λ_schedule[i])"
        println(message)
        BAOAB_allλ!(data[i], p, q,  m, params, λ_schedule[i], λ_vector, L, V, dVdr, λV, dλVdr, dλVdλ,
         γ, kbT, δt, prod_steps, sample_every)
    end

    return data
end


"""
Auxiliary function for BAR computation
"""
function BAR_aux(ΔλU_12::Vector{Float64}, ΔλU_21::Vector{Float64}, kbT::Real)
    #define function that needs to be zero
    f(ΔF) = mean(1 ./(1 .+ exp.((ΔλU_12 .- ΔF)./kbT))) - mean(1 ./(1 .+ exp.((ΔλU_21 .+ ΔF)./kbT)))
    ΔF_bar = find_zero(f, 0.0) 
    return ΔF_bar
end


"""
Compute the free energy difference using BAR
"""
function BAR(initial_eq_steps::Int, eq_steps::Int, prod_steps::Int, N::Int, m::Real, params::LJ_params, λ_steps::Int, L::SVector{d,R}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real) where {d, R<:Real}

    #generate the necessary data for both FEP and BAR computation of the free energy
    data = collect_mbar(initial_eq_steps, eq_steps, prod_steps, N, m, params, λ_steps, L, V, dVdr, λV, dλVdr, dλVdλ, 
    γ, kbT, δt)

    local_ΔF = zeros(λ_steps)

    for i in 1:λ_steps
            ΔλU_12 = data[i].λ_U[:,i+1] .- data[i].λ_U[:,i]
            ΔλU_21 = data[i+1].λ_U[:,i] .- data[i+1].λ_U[:,i+1] 
            local_ΔF[i] = BAR_aux(ΔλU_12, ΔλU_21, kbT)
    end
        
    return cumsum(local_ΔF)       
end


"""
Auxiliary function for computing the logarithm of the sum of the elements of a vector
while ensuring stability (thanks to the substraction of the maximum element in the argument of the exponential)
"""
function log_sum_exp(x::Vector{Float64}, weights::Union{Nothing, Vector{Float64}}=nothing)
    Mx = maximum(x)
    if weights === nothing
        return Mx + log(sum(exp.(x .- Mx)))
    else 
        if length(x) != length(weights)
            throw(DimensionMismatch("The dimensions of x and weights must match, but currently x has length $(length(x)) and weights has length $(length(weights))"))
        end
        return Mx + log(sum(weights .* exp.(x .- Mx)))
    end
end


"""
Auxiliary function for computation of MBAR overlap matrix
"""
function mbar_overlap(u_sk::Matrix{Float64}, f_k::Vector{Float64}, S_k::Vector{Float64}, log_denom::Vector{Float64})
    S_total, K = size(u_sk)
    C = zeros(K,K)

    for s in 1:S_total
        w_s = exp.(f_k .- u_sk[s,:] .- log_denom[s])
        for k in 1:K, l in 1:K
            C[k,l] += S_k[k]*w_s[k]*w_s[l]
        end
    end

    return C
end
    

"""
Auxiliary function for MBAR computation.
u_sk is the matrix of reduced lambda-dependent energies (i.e., lambda dependent energies divided by kbT)
u_sk must have shape (S_total, lambda_steps+1)
In our case S_total is always prod_steps*(lambda_steps+1).
S_k is the vector with the number of samples for each lambda value. 
In our case, S_k = [prod_steps ... prod_steps].
"""
function MBAR_aux(u_sk::Matrix{Float64}, S_k::Vector{Float64}; max_iter=500, tol=1e-12, α=0.2)
    S_total, K = size(u_sk)
    if S_total != sum(S_k)
        throw(DimensionMismatch("The sum of the elements in S_k must match the first dimension of u_sk"))
    end

    #f = F/kbT
    f_k = zeros(K)
    
    # warm start with picard's method (fixed point iteration)
    for iter in 1:50
        f_old = copy(f_k)

        #compute denominator in MBAR formula
        log_denom = [log_sum_exp(f_k .- u_sk[s,:], S_k) for s in 1:S_total]

        #update f_k in log_space (Picard's method)
        for k in 1:K
            f_k[k] = -log_sum_exp(-u_sk[:,k] .- log_denom)
        end

        #normalize to keep f_k[1] = 0
        f_k .-= f_k[1]

        #check convergence 
        if maximum(abs.(f_k[2:end] .- f_old[2:end])) < tol
            println("MBAR converged in $(iter) iterations")
            return f_k
        end
    end

    #Newton's method (inspired by pyMBAR)
    # R = [f_k[k] + log_sum_exp(-u_sk[:,k] .- log_denom) for k in 1:K]
    for iter in 1:max_iter
        #compute denominator in MBAR formula
        log_denom = [log_sum_exp(f_k .- u_sk[s,:], S_k) for s in 1:S_total]

        #compute residuals 
        R = [f_k[k] + log_sum_exp(-u_sk[:,k] .- log_denom) for k in 1:K]

        #compute overlap matrix
        C = mbar_overlap(u_sk, f_k, S_k, log_denom)

        #compute Jacobian
        J = I - C

        #we need to remove one state as we only have relative free energies
        #here we always remove lambda=0
        δf = zeros(K)
        δf[2:end] = -(J[2:end, 2:end]\ R[2:end])

        f_k .+= α .* δf

        #compute new denominator and residuals  
        log_denom = [log_sum_exp(f_k .- u_sk[s,:], S_k) for s in 1:S_total]
        R = [f_k[k] + log_sum_exp(-u_sk[:,k] .- log_denom) for k in 1:K]

        #check convergence 
        if maximum(abs.(R[2:end])) < tol
            println("MBAR converged in $(50 + iter) iterations")
            return f_k
        end
    end

    throw(Error("MBAR did not converge"))
end


"""
Compute the free energy difference using BAR and MBAR
"""
function MBAR(initial_eq_steps::Int, eq_steps::Int, prod_steps::Int, N::Int, m::Real, params::LJ_params, λ_steps::Int, L::SVector{d,R}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real) where {d, R<:Real}

    #generate the necessary data for both FEP and BAR computation of the free energy
    data = collect_mbar(initial_eq_steps, eq_steps, prod_steps, N, m, params, λ_steps, L, V, dVdr, λV, dλVdr, dλVdλ, 
    γ, kbT, δt)

    local_ΔF_BAR = zeros(λ_steps)

    for i in 1:λ_steps
            ΔλU_12 = data[i].λ_U[:,i+1] .- data[i].λ_U[:,i]
            ΔλU_21 = data[i+1].λ_U[:,i] .- data[i+1].λ_U[:,i+1] 
            local_ΔF_BAR[i] = BAR_aux(ΔλU_12, ΔλU_21, kbT)
    end
    ΔF_BAR = cumsum(local_ΔF_BAR)

    #flatten our data struct so that for every lambda we have all samples
    #from all simulations (ran at different lambdas) in the same column
    U_sk = vcat([d.λ_U for d in data]...)
    u_sk = U_sk ./kbT

    #the samples for every lambda are all the same in our case but
    #we assemble them in a vector for mbar implementation convenience and convert them to floats
    S_k = fill(prod_steps, λ_steps+1)
    S_k = Float64.(S_k)

    #obtain free energy differences between states using MBAR
    f_k = MBAR_aux(u_sk, S_k)
    ΔF_MBAR = kbT.*f_k

    return ΔF_BAR, ΔF_MBAR
end



