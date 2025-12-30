include("src/muFE_LJ.jl")

using .muFE_LJ
using StaticArrays
using ProgressBars
using Statistics
using Plots

#this script uses reduced units

params = LJ_params(1.0,1.0)
m = 1.0
γ = 1.0
T = 1.0
eval_steps = 20000
eq_steps = 50000
n_samples = 250
δt = 0.005
N = 500
ρ = 0.8
l = box_length(N, ρ)
rc = l/2.0
L = @SVector[l, l, l]
d = length(L)

println("Running simulation in a cubic box of reduced length $(l), at reduced temperature $(T) and reduced density $(ρ)")
println("The BAOAB integrators will have a timestep in reduced units of $(δt) and use the friction coefficient $(γ)")

#########################################################
#       Check we truly are in an NVT ensemble
#########################################################

#sample initial momentums from the Boltzmann distribution
p = [MVector{d, Float64}(randn(d) .* sqrt(m * T)) for _ in 1:N]

#fix intial positions uniformly on the box
n = ceil(Int, N^(1/d)) #grid size
grid = range(-L[1]/2, L[1]/2; length=n+1)[1:end-1]
coords = collect(Iterators.product(ntuple(_ -> grid, d)...))

q = Vector{MVector{d,Float64}}(undef, N)
for i in 1:N
    q[i] = MVector{d,Float64}(Float64.(coords[i]))
end

#equilibrate
BAOAB!(p, q, m, params, 0.0, L, LJ_pot_der, dλsoft_pot_dr, γ, T, δt, eq_steps)

#Compute correlation time for U
Us = zeros(eval_steps)
@inbounds for i in ProgressBar(1:eval_steps)
    onestep_BAOAB!(p, q, m, params, 0.0, L, LJ_pot_der, dλsoft_pot_dr, γ, T, δt)
    Us[i] = Upi(q, params, 0.0, L, LJ_pot, λsoft_pot)
end
τ = autocorr_time(Us, 2000)
println("The autocorrelation time of the potential energy in this steup is approximately $(τ)")
sample_every = ceil(Int, τ)
println("Sampling every $(sample_every) steps")


#compute U, T 
temps = zeros(n_samples)
Us = zeros(n_samples)
@inbounds for i in ProgressBar(1:n_samples)
    BAOAB!(p, q, m, params, 0.0, L, LJ_pot_der, dλsoft_pot_dr, γ, T, δt, sample_every)
    temps[i] = sum(sum(abs2, pi) for pi in p) / (d*N*m)
    Us[i] = Upi(q, params, 0.0, L, LJ_pot, λsoft_pot)
end

plot(range(1, n_samples), temps, 
    label="Kinetic Temperature",          # Label for the first line (legend entry)
    linewidth=3,                # Line thickness
    linestyle=:solid,            # Line style (e.g., :solid, :dot, :dash)
    color=:blue,                # Line color
    title="Temperature conservation check", # Main plot title
    xlabel="Time",   # X-axis label
    ylabel="T^{*}",             # Y-axis label
    legend=:bottomleft,         # Position of the legend
    #grid=true                   # Show the grid lines
)

savefig("Temperature_Conservation_Check_$(ρ)_$(T).png")


plot(range(1, n_samples), Us, 
    label="Internal Energy",          # Label for the first line (legend entry)
    linewidth=3,                # Line thickness
    linestyle=:solid,            # Line style (e.g., :solid, :dot, :dash)
    color=:blue,                # Line color
    title="Internal energy check", # Main plot title
    xlabel="Time",   # X-axis label
    ylabel="U",             # Y-axis label
    legend=:bottomleft,         # Position of the legend
    #grid=true                   # Show the grid lines
)

savefig("U_check_$(ρ)_$(T).png")

U_correction = tailcorr_U_LJ(ρ, N, rc, params)
println(mean(Us))
println(mean(Us) + U_correction)