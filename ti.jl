include("src/muFE_LJ.jl")

using .muFE_LJ

using StaticArrays
using ProgressBars
using Plots

###################################################################################
#                  CALCULATE μ_ext THROUGH TI SIMULATIONS
###################################################################################

#T = [] #values of T in reduced units
#ρ = [] #values of ρ in reduced units

#PLAN: try temperatures 0.6, 0.9, 1.2
#PLAN: try densities 0.1, 0.5, 

#this script uses reduced units

params = LJ_params(1.0,1.0)
m = 1.0
γ = 1.0
T = 1.2
ρ = 0.8
δt = 0.001
N = 256

l = box_length(N, ρ)
L = @SVector[l, l, l]
d = length(L)

λ_steps = 10
write_every = 100000


ΔF = thermodynamic_integration(N, m, params, λ_steps, L, LJ_pot, LJ_pot_der, λsoft_pot, dλsoft_pot_dr, dλsoft_pot_dλ,
 γ, T, δt, write_every)

plot(range(0, 1; length = λ_steps+1), ΔF, 
    label="Free energy difference",          # Label for the first line (legend entry)
    linewidth=3,                # Line thickness
    linestyle=:solid,            # Line style (e.g., :solid, :dot, :dash)
    color=:blue,                # Line color
    title="Free enrgy difference as a function of λ", # Main plot title
    xlabel="λ",   # X-axis label
    ylabel="ΔF",             # Y-axis label
    legend=:bottomleft,         # Position of the legend
    #grid=true                   # Show the grid lines
)

println(ΔF)

filename = "TI_$(ρ)_$(T)_$(λ_steps).png"
savefig(filename)

