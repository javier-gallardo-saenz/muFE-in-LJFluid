include("src/muFE_LJ.jl")

using .muFE_LJ

using StaticArrays
using Plots

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
T = 1.5
ρ = 0.1
δt = 0.001
N = 256

l = box_length(N, ρ)
L = @SVector[l, l, l]
d = length(L)

λ_steps = 20
write_every = 100000
initial_eq_steps = 15000
eq_steps = 10000
prod_steps = 15000


ΔF = BAR(initial_eq_steps, eq_steps, prod_steps, N, m, params, λ_steps, L, LJ_pot, LJ_pot_der, λsoft_pot, dλsoft_pot_dr, dλsoft_pot_dλ,
 γ, T, δt, write_every)

ΔF = cat(0.0, ΔF; dims=1)

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

println(ΔF)

filename = "bar_$(ρ)_$(T)_$(λ_steps).png"
savefig(filename)