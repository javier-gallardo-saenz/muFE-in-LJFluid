using StaticArrays
using random

using integrators
using utils


"""
Perform thermodynamic integration (TI)
"""
function thermodynamic_integration(N::Int, m::Real, params::LJ_params, λ_steps::Int, L::SVector{d,R}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real, write_every::Int) where {d, R<:Real}

    #define λ schedule
    λ_schedule = range(0, 1; length = λ_steps+1)

    #sample initial momentums from the Boltzmann distribution
    p = [@MVector randn(d) .* sqrt(m * kbT) for _ in 1:N]

    #fix intial positions uniformly on the box
    n = ceil(Int, N^(1/d)) #grid size
    grid = range(0,L[1];length=n+1)[1:end-1] 
    coords = collect(Iterators.product(ntuple(_->grid,d)...))

    q = Vector{MVector{d,R}}(undef,N)
    for i in 1:N
        q[i] = @MVector(FLoat64.(coords[i]))
    end

    #number of steps we can be sure will be enough to ensure equilibrium is reached
    safe_no_steps_for_eq = 1000
    safe_no_steps_for_prod = 1000

    #vectors to store <∂H/∂λ>
    dHdλ_avgs = zeros(length(λ_schedule))

    for i in eachindex(λ_schedule)
        #equilibrate
        BAOAB!(p, q, m, params, λ_schedule[i], L, dVdr, dλVdr, γ, kbT, δt, safe_no_steps_for_eq, write_every)
        #production run
        eH_avgs[i], dHdλ_avgs[i] = HBAOAB_E!(p,q, params, λ_schedule[i], L, dVdr, dλVdr, γ, kbT, δt, safe_no_steps_for_prod, write_every)
    
    
    #for visualization
    μ_ext_ti = zeros(length(λ_schedule))
    for i in 2:(λ_steps + 1)
        δλ = λ_schedule[i] - λ_schedule[i-1]
        μ_ext_ti[i] = μ_ext_ti[i-1] + 0.5*δλ*(dHdλ_avgs[i] - dHdλ_avgs[i-1])
    
    return μ_ext_ti


"""
Collect data for FEP, BAR
"""
function collect_fep_bar(N::Int, m::Real, params::LJ_params, λ_steps::Int, L::SVector{d,R}, 
    V::Function, dVdr::Function, λV::Function, dλVdr::Function, dλVdλ::Function,
    γ::Real, kbT::Real, δt::Real, write_every::Int) where {d, R<:Real}



    

    
    









end