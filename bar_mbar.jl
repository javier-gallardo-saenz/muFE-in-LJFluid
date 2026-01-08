include("src/muFE_LJ.jl")

using .muFE_LJ

using StaticArrays
using Plots
using JLD2

###################################################################################
#                  CALCULATE μ_ext THROUGH BAR 
###################################################################################

#T = [] #values of T in reduced units
#ρ = [] #values of ρ in reduced units

#PLAN: try temperatures 0.6, 1.2, 1.5
#PLAN: try densities 0.1, 0.8, 0.9 

#this script uses reduced units

params = LJ_params(1.0,1.0)
m = 1.0
γ = 1.0
T = 1.0
ρ = 0.8
δt = 0.002
N = 256

l = box_length(N, ρ)
L = @SVector[l, l, l]
d = length(L)
rc = l/2.0
println("Running simulation in a cubic box of reduced length $(l), at reduced temperature $(T) and reduced density $(ρ)")
println("The BAOAB integrators will have a timestep in reduced units of $(δt) and use the friction coefficient $(γ)")

λ_steps = 50
initial_eq_steps = 50000
eq_steps = 15000
prod_steps = 20000

#comment this if data has already been collected

#data_mbar = collect_mbar(initial_eq_steps, eq_steps, prod_steps, N, m, params, λ_steps, L, LJ_pot, LJ_pot_der, λsoft_pot, dλsoft_pot_dr, dλsoft_pot_dλ,
# γ, T, δt)
#@save "MBAR_data/data_$(ρ)_$(T)_$(λ_steps).jld2" data_mbar

@load "MBAR_data/data_$(ρ)_$(T)_$(λ_steps).jld2" data_mbar

ΔF_BAR, σexp_ΔF_BAR, σ_ΔF_BAR, ΔF_MBAR, σ_ΔF_MBAR, σ_ΔF_MBAR, τs = MBAR(initial_eq_steps, eq_steps, prod_steps, N, m, params, λ_steps, L, LJ_pot, LJ_pot_der, λsoft_pot, dλsoft_pot_dr, dλsoft_pot_dλ,
γ, T, δt, data_mbar)

println("The free energy difference as predicted by BAR is")
println(ΔF_BAR)
println("The theoretical, asymptotic std is")
println(σ_ΔF_BAR)

println("The free energy difference as predicted by BAR is")
println(ΔF_MBAR)
println("The theoretical, asymptotic std is")
println(σ_ΔF_MBAR)

tail_correction_muext = tailcorr_μex_LJ(ρ, rc, params)
println("The tail correction is")
println(tail_correction_muext)

p_BAR = plot(range(0, 1; length = λ_steps+1), ΔF_BAR, 
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

plot!(range(0, 1; length = λ_steps+1), ΔF_BAR .+ σ_ΔF_BAR, fillrange=ΔF_BAR .- σ_ΔF_BAR, linewidth=0, fillalpha=0.3, label="±1σ (asymptotic)")


filename = "BAR_plots/bar_$(ρ)_$(T)_$(λ_steps).png"
savefig(p_BAR, filename)


p_MBAR = plot(range(0, 1; length = λ_steps+1), ΔF_MBAR, 
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

plot!(range(0, 1; length = λ_steps+1), ΔF_BAR .+ σ_ΔF_MBAR, fillrange=ΔF_BAR .- σ_ΔF_MBAR, linewidth=0, fillalpha=0.3, label="±1σ (asymptotic)")

filename = "MBAR_plots/mbar_$(ρ)_$(T)_$(λ_steps).png"
savefig(p_MBAR, filename)

println("The autocorrelation times of λU were")
println(τs)




