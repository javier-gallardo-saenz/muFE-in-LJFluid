include("src/muFE_LJ.jl")

using .muFE_LJ

using StaticArrays
using ProgressBars
using Plots
using JLD2

###################################################################################
#                  CALCULATE μ_ext THROUGH TI SIMULATIONS
###################################################################################

#this script uses reduced units

params = LJ_params(1.0,1.0)
m = 1.0
γ = 1.0
T = 1.0
ρ = 0.9
δt = 0.002
N = 256

l = box_length(N, ρ)
L = @SVector[l, l, l]
d = length(L)
rc = l/2.0
println("Running simulation in a cubic box of reduced length $(l), at reduced temperature $(T) and reduced density $(ρ)")
println("The BAOAB integrators will have a timestep in reduced units of $(δt) and use the friction coefficient $(γ)")

λ_steps = 100
initial_eq_steps = 50000
eq_steps = 15000
prod_steps = 20000

#comment this if data has already been collected
#data_ti = collect_ti(initial_eq_steps, eq_steps, prod_steps, N, m, params, λ_steps, L, LJ_pot, LJ_pot_der, λsoft_pot, dλsoft_pot_dr, dλsoft_pot_dλ,
# γ, T, δt)
#@save "TI_data/data_$(ρ)_$(T)_$(λ_steps).jld2" data_ti

@load "TI_data/data_$(ρ)_$(T)_$(λ_steps).jld2" data_ti

ΔF, σ_ΔF, τs = thermodynamic_integration(initial_eq_steps, eq_steps, prod_steps, N+1, m, params, λ_steps, L, 
 LJ_pot, LJ_pot_der, λsoft_pot, dλsoft_pot_dr, dλsoft_pot_dλ, 
 γ, T, δt, data_ti)

plot(range(0, 1; length = λ_steps+1), ΔF, 
    label="Free energy difference",          # Label for the first line (legend entry)
    linewidth=3,                # Line thickness
    linestyle=:solid,            # Line style (e.g., :solid, :dot, :dash)
    color=:blue,                # Line color
    title="Free energy difference as a function of λ", # Main plot title
    xlabel="λ",   # X-axis label
    ylabel="ΔF",             # Y-axis label
    legend=:bottomleft,         # Position of the legend
    #grid=true                   # Show the grid lines
)

plot!(range(0, 1; length = λ_steps+1), ΔF .+ σ_ΔF, fillrange=ΔF .- σ_ΔF, linewidth=0, fillalpha=0.3, label="±1σ")

println("The free energy difference as predicted by TI is")
println(ΔF)
println("The std predicted by block averaging is")
println(σ_ΔF)
println("The autocorrelation times were")
println(τs)

filename = "TI_plots/TI_$(ρ)_$(T)_$(λ_steps).png"
savefig(filename)

