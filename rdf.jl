include("src/muFE_LJ.jl")

using .muFE_LJ
using StaticArrays
using ProgressBars
using Plots

#this script uses reduced units

params = LJ_params(1.0,1.0)
m = 1.0
γ = 1.0
T = 1.0
δt = 0.001
N = 256
ρ = 0.8

l = box_length(N, ρ)
L = @SVector[l, l, l]
d = length(L)

sample_every = 10
write_every = 1000000
λ = 0.5
λ_steps = 5

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

#equilibrate at λ = 0
BAOAB!(p, q, m, params, 0.0, L, LJ_pot_der , dλsoft_pot_dr, γ, T, δt, 5000, write_every)
println(maximum(q))

#bring system to desired λ
λ_schedule = range(0, λ; length = λ_steps+1)
for i in ProgressBar(eachindex(λ_schedule))
    BAOAB!(p, q, m, params, λ_schedule[i], L, LJ_pot_der , dλsoft_pot_dr, γ, T, δt, 10000, write_every)
    println(maximum(q))
end

#define bins
n_bins = 100
hist = zeros(Float64, n_bins)
dr = sqrt(sum(abs2, L./2))/n_bins
BAOAB_rdf!(hist, dr, p, q, m, params, λ, L, LJ_pot_der , dλsoft_pot_dr, γ, T, δt, 15000, sample_every)
println(maximum(q))

#normalize
npairs = N*(N-1)/2
g = zeros(Float64, n_bins)
rs = [(k-0.5)*dr for k in 1:n_bins]  # midpoint of each bin

for k in 1:n_bins
    r = rs[k]
    shell_volume = (4π/3) * ((r+dr)^3 - r^3)
    expected = ρ * shell_volume * npairs * T
    g[k] = hist[k] / expected
end

plot(rs, g, 
    label="Radial Distribution Function",          # Label for the first line (legend entry)
    linewidth=3,                # Line thickness
    linestyle=:solid,            # Line style (e.g., :solid, :dot, :dash)
    color=:blue,                # Line color
    title="Radial distribution function", # Main plot title
    xlabel="r",   # X-axis label
    ylabel="g(r)",             # Y-axis label
    legend=:bottomleft,         # Position of the legend
    #grid=true                   # Show the grid lines
)

savefig("rdf.png")

