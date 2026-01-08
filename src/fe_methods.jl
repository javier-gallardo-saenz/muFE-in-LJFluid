using StaticArrays
using Random
using ProgressBars
using Roots
using Statistics
using LinearAlgebra


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
Collect data for thermodynamic integration
"""
function collect_ti(initial_eq_steps::Int, eq_steps::Int, prod_steps::Int, N::Int, m::Real, params::LJ_params, λ_steps::Int, L::SVector{d,R}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real) where {d, R<:Real}

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

    #data struct to store ∂H/∂λ for every sample in each lambda run
    data = [mu_data(Vector{Float64}(), 1.0) for _ in 1:λ_steps+1]

    println("Running MD for the different λ values")
    for i in ProgressBar(eachindex(λ_schedule))
        #equilibrate
        println("Equilibrating")
        BAOAB!(p, q, m, params, λ_schedule[i], L, dVdr, dλVdr, γ, kbT, δt, eq_steps)
        #production run
        println("Production run for λ=$(λ_schedule[i])")
        BAOAB_ti!(data[i], p, q, m, params, λ_schedule[i], L, V, dVdr, λV, dλVdr, dλVdλ, γ, kbT, δt, prod_steps, 1)
        # compute autocorrelation time for the current lambda
        data[i].τ = autocorr_time(data[i].w, 2000)
    end

    return data
end


"""
Perform thermodynamic integration (TI)
"""
function thermodynamic_integration(initial_eq_steps::Int, eq_steps::Int, prod_steps::Int, N::Int, m::Real, params::LJ_params, λ_steps::Int, L::SVector{d,R}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real, data::Union{Nothing, Vector{mu_data}}=nothing) where {d, R<:Real}

    if data === nothing
        #generate the necessary data for both FEP and BAR computation of the free energy
        data = collect_ti(initial_eq_steps, eq_steps, prod_steps, N, m, params, λ_steps, L, V, dVdr, λV, dλVdr, dλVdλ, 
        γ, kbT, δt)
        if length(data) != λ_steps+1
            throw(DimensionMismatch("The length of the data vector ($(length(data)) should match the number of λ values ($(λ_steps+1))"))
        end
    end

    #vector to store std(∂H/∂λ)
    dHdλ_std = zeros(λ_steps+1)
    for i in eachindex(data)
        _ , dHdλ_std[i] = block_average(data[i].w, 2*ceil(Int,data[i].τ))
    end

    #recover λ schedule
    λ_schedule = range(0, 1; length = λ_steps+1)

    #free energy computation with trapezoidal rule
    #uncertainty computation using error propagation for linear combination
    μ_ext = zeros(λ_steps+1)
    σ_μ_ext = zeros(λ_steps+1)
    for i in 2:(λ_steps + 1)
        δλ = λ_schedule[i] - λ_schedule[i-1]
        μ_ext[i] = μ_ext[i-1] + 0.5*δλ*(mean(data[i].w) + mean(data[i-1].w))
        σ_μ_ext[i] = σ_μ_ext[i-1] + (0.5*δλ)^(2)*(dHdλ_std[i]^2 + dHdλ_std[i-1]^2)
    end
    σ_μ_ext = sqrt.(σ_μ_ext)

    #recover autocorrelation times for each lambda
    τs = [data[i].τ for i in eachindex(data)]

    return μ_ext, σ_μ_ext, τs
end


"""
Collect data for FEP
"""
function collect_fep(initial_eq_steps::Int, eq_steps::Int, prod_steps::Int, N::Int, m::Real, params::LJ_params, λ_steps::Int, L::SVector{d,R}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real) where {d, R<:Real}
    
    data= [mu_data(Vector{Float64}(), 1.0) for _ in 1:λ_steps]

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

    #create auxiliary vectors to store samples of U at the different lambdas and to store the autocorrelation time for each lambda
    auxΔU = zeros(prod_steps)

    for i in ProgressBar(1:λ_steps)
        #equilibrate
        println("Equilibrating")
        BAOAB!(p, q, m, params, λ_schedule[i], L, dVdr, dλVdr, γ, kbT, δt, eq_steps)
        #production run
        message = "Production run for λ=$(λ_schedule[i])"
        println(message)
        BAOAB_twoλ!(auxΔU, p, q, m, params, λ_schedule[i], λ_schedule[i+1], L, V, dVdr, λV, dλVdr, dλVdλ, 
        γ, kbT, δt, prod_steps, 1)
        #compute autocorrelation time of our observable exp(-auxΔU./kbT)
        data[i].w = exp.(-auxΔU./kbT)
        data[i].τ = autocorr_time(data[i].w, 2000)
    end

    return data
end


"""
Compute the free energy difference using FEP
"""
function fep(initial_eq_steps::Int, eq_steps::Int, prod_steps::Int, N::Int, m::Real, params::LJ_params, λ_steps::Int, L::SVector{d,R}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real, data::Union{Nothing, Vector{mu_data}}=nothing) where {d, R<:Real}

    if data === nothing
        #generate the necessary data for both FEP and BAR computation of the free energy
        data = collect_fep(initial_eq_steps, eq_steps, prod_steps, N, m, params, λ_steps, L, V, dVdr, λV, dλVdr, dλVdλ, 
        γ, kbT, δt)
        if length(data) != λ_steps+1
            throw(DimensionMismatch("The length of the data vector ($(length(data)) should match the number of λ values ($(λ_steps+1))"))
        end
    end

    #define FEP function for bootstrapping
    f(w::Vector{Float64}) = -kbT*log(mean(w))
    
    #bootstrapping
    ΔF = zeros(λ_steps+1)
    ΔF_σ = zeros(λ_steps+1)
    for i in 1:λ_steps
        _, σ = bootstrap(block_averaged_samples(data[i].w, 2*ceil(Int, data[i].τ)), f, 1000)
        ΔF[i+1] = ΔF[i] - kbT*log(mean(data[i].w))
        ΔF_σ[i+1] = sqrt(ΔF_σ[i+1]^2 + σ^2)
    end

    #recover autocorrelation times for each lambda
    τs = [data[i].τ for i in eachindex(data)]

    return ΔF, ΔF_σ, τs
end


"""
Collect data for MBAR (used for BAR as well)
"""
function collect_mbar(initial_eq_steps::Int, eq_steps::Int, prod_steps::Int, N::Int, m::Real, params::LJ_params, λ_steps::Int, L::SVector{d,R}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real) where {d, R<:Real}
    
    data = [mu_data_MBAR(zeros(prod_steps, λ_steps+1), 1.0) for _ in 1:(λ_steps+1)]
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
    autocorr_comp_steps = 20000
    Us = zeros(autocorr_comp_steps)
    println("Collecting data to estimate the autocorrelation time of the potential energy for the λ=0 state")
    @inbounds for i in ProgressBar(1:autocorr_comp_steps)
        onestep_BAOAB!(p, q, m, params, 0.0, L, dVdr, dλVdr, γ, kbT, δt)
        Us[i] = Upi(q, params, 0.0, L, V, λV) 
    end
    τ0 = autocorr_time(Us, 2000)

    for i in ProgressBar(1:λ_steps+1)
        #equilibrate
        println("Equilibrating")
        BAOAB!(p, q, m, params, λ_schedule[i], L, dVdr, dλVdr, γ, kbT, δt, eq_steps)
        #production run
        message = "Production run for λ=$(λ_schedule[i])"
        println(message)
        BAOAB_allλ!(data[i], p, q,  m, params, λ_schedule[i], λ_vector, L, V, dVdr, λV, dλVdr, dλVdλ,
         γ, kbT, δt, prod_steps, 1)
        #store autocorrelation time, which we compute from the samples taken at the lambda value of each simulation (except for the case λ=0)
        if i ==1
            data[i].τ = τ0
        else
            data[i].τ = autocorr_time(data[i].λ_U[:,i], 2000) 
        end
    end

    return data
end


"""
Auxiliary function for BAR computation. Assumes independent samples
"""
function BAR_aux(ΔλU_12::Vector{Float64}, ΔλU_21::Vector{Float64}, kbT::Real)
    #define function that needs to be zero
    N1 = length(ΔλU_12)
    N2 = length(ΔλU_21)
    aux = kbT*log(N1/N2)
    f(ΔF) = mean(1 ./(1 .+ exp.((ΔλU_12 .- ΔF .+ aux)./kbT))) - mean(1 ./(1 .+ exp.((ΔλU_21 .+ ΔF .- aux)./kbT)))
    ΔF_bar = find_zero(f, 0.0) 
    return ΔF_bar
end


"""
Auxiliary function for BAR computation with bootstrapping
"""
function BAR_aux_bootstrap(λU_resampled::Vector{Matrix{Float64}}, i::Int, kbT::Real)
    #extract the data we need
    ΔλU_12 = λU_resampled[i][:,i+1] .- λU_resampled[i][:,i]
    ΔλU_21 = λU_resampled[i+1][:,i] .- λU_resampled[i+1][:,i+1]
    #define function that needs to be zero
    N1 = length(ΔλU_12)
    N2 = length(ΔλU_21)
    aux = kbT*log(N1/N2)
    f(ΔF) = mean(1 ./(1 .+ exp.((ΔλU_12 .- ΔF .+ aux)./kbT))) - mean(1 ./(1 .+ exp.((ΔλU_21 .+ ΔF .- aux)./kbT)))
    ΔF_bar = find_zero(f, 0.0) 
    return ΔF_bar
end


"""
Auxiliary function for BAR theoretical, asymptotic σ computation. Assumes independent samples.
"""
function BAR_σ(ΔF::Float64, ΔλU_12::Vector{Float64}, ΔλU_21::Vector{Float64}, kbT::Real)
    #define function that needs to be zero
    N1 = length(ΔλU_12)
    N2 = length(ΔλU_21)
    aux = kbT*log(N1/N2)
    w1 = 1 ./(1 .+ exp.((ΔλU_12 .- ΔF .+ aux)./kbT))
    τ_w1 = autocorr_time(w1, 2000)
    w2 = 1 ./(1 .+ exp.((ΔλU_21 .+ ΔF .- aux)./kbT))
    τ_w2 = autocorr_time(w2, 2000)
    var = 1 / (sum(w1 .* (1 .- w1))/(2*τ_w1) + sum(w2 .* (1 .- w2))/(2*τ_w2))
    return sqrt(var)
end


"""
Compute the free energy difference using BAR. Assumes independent samples.
"""
function BAR(initial_eq_steps::Int, eq_steps::Int, prod_steps::Int, N::Int, m::Real, params::LJ_params, λ_steps::Int, L::SVector{d,R}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real, data::Union{Nothing, Vector{mu_data_MBAR}}=nothing) where {d, R<:Real}

    if data === nothing
        #generate the necessary data for both FEP and BAR computation of the free energy
        data = collect_mbar(initial_eq_steps, eq_steps, prod_steps, N, m, params, λ_steps, L, V, dVdr, λV, dλVdr, dλVdλ, 
        γ, kbT, δt)
        if length(data) != λ_steps+1
            throw(DimensionMismatch("The length of the data vector ($(length(data)) should match the number of λ values ($(λ_steps+1))"))
        end
    end

    #compute data blocks for bootstrapping
    data_blocks = Vector{Vector{Matrix{Float64}}}(undef, λ_steps+1)
    for i in eachindex(data)
        data_blocks[i] = blocks_for_mbar(data[i].λ_U, 2*ceil(Int, data[i].τ))
    end

    ΔF = zeros(λ_steps+1)
    σ_ΔF = zeros(λ_steps+1)
    σexp_ΔF = zeros(λ_steps+1)

    for i in 1:λ_steps
        idx = i
        f(λU_resampled::Vector{Matrix{Float64}}) = BAR_aux_bootstrap(λU_resampled, idx, kbT)
        μ, σ = block_bootstrap_bar(data_blocks, f, 1000)
        ΔF[i+1] = ΔF[i] + μ
        σexp_ΔF[i+1] = sqrt(σexp_ΔF[i]^2 + σ^2)

        #compute asymptotic std
        ΔλU_12 = data[i].λ_U[:,i+1] .- data[i].λ_U[:,i]
        ΔλU_21 = data[i+1].λ_U[:,i] .- data[i+1].λ_U[:,i+1]
        σ_asymp = BAR_σ(μ, ΔλU_12, ΔλU_21, kbT)
        σ_ΔF[i+1] = sqrt(σ_ΔF[i]^2 + σ_asymp^2)
    end

    #recover autocorrelation times for each lambda
    τs = [data[i].τ for i in eachindex(λ_schedule)]
        
    return ΔF, σ_ΔF, τs       
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
Auxiliary function for MBAR computation. Assumes independent samples.
u_sk is the matrix of reduced lambda-dependent energies (i.e., lambda dependent energies divided by kbT)
u_sk must have shape (S_total, lambda_steps+1)
"""
function MBAR_aux(u_sk::Matrix{Float64}, S_k::Vector{Float64}; max_iter::Int=500, tol::Float64=1e-12, α::Float64=0.5, f_k_0::Union{Nothing, Vector{Float64}}=nothing)
    S_total, K = size(u_sk)
    if S_total != sum(S_k)
        throw(DimensionMismatch("The sum of the elements in S_k must match the first dimension of u_sk"))
    end

    #f = F/kbT
    if f_k_0 !== nothing
        if K != length(f_k_0)
            throw(DimensionMismatch("The length of the provided initial guess for the reduced free energies ($(length(f_k_0))) must match the number of λ states ($(K))"))
        end
        f_k = copy(f_k_0)
    else 
        f_k = zeros(K)
    end
    
    # warm start with picard's method (fixed point iteration)
    println("Warmstarting MBAR with Picard's method")
    for iter in ProgressBar(1:50)
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
    println("Solving MBAR with Newton's method")
    for iter in ProgressBar(1:max_iter)
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

    throw(ErrorException("MBAR did not converge"))
end


"""
Auxiliary function for computation of MBAR covariance matrix
Assumes independent samples
"""
function mbar_cov(u_sk::Matrix{Float64}, f_k::Vector{Float64}, S_k::Vector{Float64}, log_denom::Vector{Float64}, tol::Float64=1e-9)
    S_total, K = size(u_sk)
    W = zeros(S_total,K)

    for s in 1:S_total
        W[s,:] = exp.(f_k .- u_sk[s,:] .- log_denom[s])
    end

    N = Diagonal(S_k)

    # Eigen-decomposition of W'W
    F = eigen(Symmetric(W' * W))
    S2 = max.(F.values, 0.0)
    σ = sqrt.(S2)

    # Threshold small singular values
    σ .= ifelse.(σ .< tol, 0.0, σ)

    Σ = Diagonal(σ)
    V = F.vectors

    # A = I - Σ V' N V Σ
    A = I - Σ * (V' * N * V) * Σ
    A_pinv = pinv(A)

    # Final covariance
    Θ = V * Σ * A_pinv * Σ * V'

    return Θ
end


"""
Compute the free energy difference using BAR and MBAR
"""
function MBAR(initial_eq_steps::Int, eq_steps::Int, prod_steps::Int, N::Int, m::Real, params::LJ_params, λ_steps::Int, L::SVector{d,R}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real, data::Union{Nothing, Vector{mu_data_MBAR}}=nothing) where {d, R<:Real}

    if data === nothing
        #generate the necessary data for both FEP and BAR computation of the free energy
        data = collect_mbar(initial_eq_steps, eq_steps, prod_steps, N, m, params, λ_steps, L, V, dVdr, λV, dλVdr, dλVdλ, 
        γ, kbT, δt)
        if length(data) != λ_steps+1
            throw(DimensionMismatch("The length of the data vector ($(length(data)) should match the number of λ values ($(λ_steps+1))"))
        end
    end

    S_k_subsampled = [div(prod_steps, ceil(Int, 2*data[i].τ)) for i in eachindex(data)]
    indices = [range(1, step=ceil(Int, 2*data[i].τ), length=S_k_subsampled[i]) for i in eachindex(data)]

    #BAR
    ΔF_BAR = zeros(λ_steps+1)
    σ_ΔF_BAR = zeros(λ_steps+1)

    println("Running BAR at each lambda interval")
    for i in ProgressBar(1:λ_steps)
        ΔλU_12 = data[i].λ_U[:,i+1] .- data[i].λ_U[:,i]
        ΔλU_21 = data[i+1].λ_U[:,i] .- data[i+1].λ_U[:,i+1]

        local_ΔF = BAR_aux(ΔλU_12, ΔλU_21, kbT)
        ΔF_BAR[i+1] = ΔF_BAR[i] + local_ΔF

        #compute asymptotic std with subsampled data
        ΔλU_12_subsampled = data[i].λ_U[indices[i],i+1] .- data[i].λ_U[indices[i],i]
        ΔλU_21_subsampled = data[i+1].λ_U[indices[i+1],i] .- data[i+1].λ_U[indices[i+1],i+1]

        σ_asymp = BAR_σ(local_ΔF, ΔλU_12_subsampled, ΔλU_21_subsampled, kbT)
        σ_ΔF_BAR[i+1] = sqrt(σ_ΔF_BAR[i]^2 + σ_asymp^2)
    end

    #MBAR 
    u_sk = vcat([data[i].λ_U ./ kbT for i in eachindex(data)]...)
    S_total = sum(S_k_subsampled)
    S_k_subsampled = Float64.(S_k_subsampled)

    S_k = fill(prod_steps, λ_steps+1)
    S_k = Float64.(S_k)

    f_k = MBAR_aux(u_sk, S_k, max_iter=10000)
    ΔF_MBAR = kbT .* f_k

    #compute theoretical, asymptotic standard deviation using covariance matrix on susbsampled data (insufficient memory otherwise)
    u_sk = vcat([data[i].λ_U[indices[i],:] ./ kbT for i in eachindex(data)]...)
    log_denom = [log_sum_exp(f_k .- u_sk[s,:], S_k_subsampled) for s in 1:S_total]
    Θ =  mbar_cov(u_sk, f_k, S_k_subsampled, log_denom)
    #because we set f_0 = 0
    σ_ΔF_MBAR = kbT .* sqrt.(diag(Θ))

    #recover autocorrelation times for each lambda
    τs = [data[i].τ for i in eachindex(λ_schedule)]

    return ΔF_BAR, σ_ΔF_BAR, ΔF_MBAR, σ_ΔF_MBAR, τs
end


"""
Compute bootstrapped BAR and MBAR. Too computationally demanding, NOT USED for very fine lambda schedules.
"""
function bootstrapped_MBAR(initial_eq_steps::Int, eq_steps::Int, prod_steps::Int, N::Int, m::Real, params::LJ_params, λ_steps::Int, L::SVector{d,R}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real, data::Union{Nothing, Vector{mu_data_MBAR}}=nothing) where {d, R<:Real}

    if data === nothing
        #generate the necessary data for both FEP and BAR computation of the free energy
        data = collect_mbar(initial_eq_steps, eq_steps, prod_steps, N, m, params, λ_steps, L, V, dVdr, λV, dλVdr, dλVdλ, 
        γ, kbT, δt)
        if length(data) != λ_steps+1
            throw(DimensionMismatch("The length of the data vector ($(length(data)) should match the number of λ values ($(λ_steps+1))"))
        end
    end

    #BAR with bootstrapping

    #compute data blocks for bootstrapping
    data_blocks = Vector{Vector{Matrix{Float64}}}(undef, λ_steps+1)
    for i in eachindex(data)
        data_blocks[i] = blocks_for_mbar(data[i].λ_U, 2*ceil(Int, data[i].τ))
    end

    ΔF_BAR = zeros(λ_steps+1)
    σ_ΔF_BAR = zeros(λ_steps+1)
    σexp_ΔF_BAR = zeros(λ_steps+1)

    println("Running BAR at each lambda interval")
    for i in ProgressBar(1:λ_steps)
        idx = i
        f(λU_resampled::Vector{Matrix{Float64}}) = BAR_aux_bootstrap(λU_resampled, idx, kbT)
        μ, σ = block_bootstrap_bar(data_blocks, f, 1000)
        ΔF_BAR[i+1] = ΔF_BAR[i] + μ
        σexp_ΔF_BAR[i+1] = sqrt(σexp_ΔF_BAR[i]^2 + σ^2)

        #compute asymptotic std
        ΔλU_12 = data[i].λ_U[:,i+1] .- data[i].λ_U[:,i]
        ΔλU_21 = data[i+1].λ_U[:,i] .- data[i+1].λ_U[:,i+1]
        σ_asymp = BAR_σ(μ, ΔλU_12, ΔλU_21, kbT)
        σ_ΔF_BAR[i+1] = sqrt(σ_ΔF_BAR[i]^2 + σ_asymp^2)
    end


    #MBAR with small bootstrapping

    #compute data blocks for bootstrapping
    data_blocks = Vector{Vector{Matrix{Float64}}}(undef, λ_steps+1)
    for i in eachindex(data)
        data_blocks[i] = blocks_for_mbar(data[i].λ_U ./ kbT, 2*ceil(Int, data[i].τ))
    end
    
    #the number of samples that bootstrapping will see for each run is num_blocks*2τ, where we have div(prod_steps, 2τ) blocks 
    S_k = [length(data_blocks[i])*length(data_blocks[i][1][:,1]) for i in eachindex(data_blocks)]
    S_total = sum(S_k)
    S_k = Float64.(S_k)

    #bootstrapping as implemented in block_bootstrap_mbar, but adapted for advantageous initial guessing 
    n_boot = 30
    K = length(data_blocks)
    bootstrapped_obs = zeros(n_boot, K)
    aux_f_k = zeros(K)
    num_blocks = [length(data_blocks[i]) for i in eachindex(data_blocks)]
    
    println("Running MBAR bootstrapping")
    for i in ProgressBar(1:n_boot)
        #sample blocks 
        indexes = [rand(1:num_blocks[k], num_blocks[k]) for k in 1:K]
        resample = [vcat(data_blocks[k][indexes[k]]...) for k in 1:K]  #we obtain vector of K matrices of size (S_k, K)
        u_sk_resampled = vcat(resample...) #we obtain matrix of size (S_tot, K)
        bootstrapped_obs[i,:] = MBAR_aux(u_sk_resampled, S_k, f_k_0=aux_f_k)
        aux_f_k = copy(bootstrapped_obs[i,:])
    end

    f_k = [mean(bootstrapped_obs[:,i]) for i in 1:K]
    ΔF_MBAR = kbT .* f_k
    σexp_ΔF_MBAR = [kbT * std(bootstrapped_obs[:,i]) for i in 1:K]

    #compute theoretical, asymptotic standard deviation using covariance matrix on susbsampled data (insufficient memory otherwise)
    S_k = [div(prod_steps, ceil(Int, 2*data[i].τ)) for i in eachindex(data)]
    indices = [range(1, step=ceil(Int, 2*data[i].τ), length=S_k[i]) for i in eachindex(data)]
    u_sk = vcat([data[i].λ_U[indices[i],:] ./ kbT for i in eachindex(data)]...)
    S_total = sum(S_k)
    S_k = Float64.(S_k)

    log_denom = [log_sum_exp(f_k .- u_sk[s,:], S_k) for s in 1:S_total]
    Σ =  mbar_cov(u_sk, f_k, S_k, log_denom)
    σ_ΔF_MBAR = kbT .* sqrt.(diag(Σ))

    #recover autocorrelation times for each lambda
    τs = [data[i].τ for i in eachindex(λ_schedule)]

    return ΔF_BAR, σexp_ΔF_BAR, σ_ΔF_BAR, ΔF_MBAR, σexp_ΔF_MBAR,  σ_ΔF_MBAR, τs
end




