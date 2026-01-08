include("src/muFE_LJ.jl")

using .muFE_LJ

using StaticArrays
using ProgressBars
using Plots

###################################################################################
#                  CALCULATE μ_ext THROUGH WIDOM METHOD
###################################################################################

#T = [] #values of T in reduced units
#ρ = [] #values of ρ in reduced units

#PLAN: try temperatures 0.6, 0.9, 1.2
#PLAN: try densities 0.1, 0.5, 

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

sample_every = 1000
n_trial_insertions = 1000
n_tries_per_insertion = 500

println("Running simulation in a cubic box of reduced length $(l), at reduced temperature $(T) and reduced density $(ρ)")
println("The BAOAB integrators will have a timestep in reduced units of $(δt) and use the friction coefficient $(γ)")

ΔF, dUs = widom_method(sample_every, n_trial_insertions, n_tries_per_insertion, N, m, params, L, LJ_pot, LJ_pot_der, dλLJ_dr,
γ, T, δt)

println("The free energy difference predicted by Widom's method is")
println(ΔF)

message_label = "ΔU Distribution at T*=$(T), ρ*=$(ρ)"

# --- Plotting the histogram ---
plot_hist = histogram(dUs, 
    # Use a specific number of bins for resolution (e.g., 50 to 100)
    bins = 100, 
    # Label the plot and axes
    label = "message_label",
    xlabel = "Insertion Energy (ΔU/kBT)", 
    ylabel = "Count",
    # Recommended: Use a log-scale on the y-axis to see the small, important tails
    yguide = "Frequency (log scale)",
    yscale = :log10, 
    title = "Widom Insertion Energy Analysis"
)

savefig(plot_hist, "Widom_plots/dU_hist_$(T).png")

