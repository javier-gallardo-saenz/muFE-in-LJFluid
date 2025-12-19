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
T = 1.5
ρ = 0.1
δt = 0.001
N = 256

l = box_length(N, ρ)
L = @SVector[l, l, l]
d = length(L)

sample_every = 100
n_trial_insertions = 10000
n_tries_per_insertion = 500

write_every = 100000


ΔF, dUs = widom_method(sample_every, n_trial_insertions, n_tries_per_insertion, N, m, params, L, LJ_pot, LJ_pot_der, dλLJ_dr,
γ, T, δt, write_every)

println(maximum(dUs))

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

savefig(plot_hist, "dU_hist.png")

