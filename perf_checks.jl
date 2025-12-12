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
T = 1.5
steps = 20000
δt = 0.001
N = 256
ρ = 0.1
l = box_length(N, ρ)
L = @SVector[l, l, l]
d = length(L)
write_every = 200


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

temps = zeros(steps+1)
Us = zeros(steps+1)

#equilibrate
BAOAB!(p, q, m, params, 0.0, L, LJ_pot_der, dλsoft_pot_dr, γ, T, δt, 20000, write_every)

@inbounds for i in ProgressBar(1:steps)
    temps[i] = sum(sum(abs2, pi) for pi in p) / (d*N*m)
    Us[i] = Upi(q, params, 0.0, L, LJ_pot, λsoft_pot)
    BAOAB!(p, q, m, params, 0.0, L, LJ_pot_der, dλsoft_pot_dr, γ, T, δt, 1, write_every)
end
temps[steps+1] = sum(sum(abs2, pi) for pi in p) / (d*N*m)
Us[steps+1] = Upi(q, params, 0.0, L, LJ_pot, λsoft_pot)

plot(range(0, steps*δt, length=steps+1), temps, 
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

savefig("Temperature_Conservation_Check.png")


plot(range(0, steps*δt, length=steps+1), Us, 
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

savefig("U_check.png")

println(mean(Us))