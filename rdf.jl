include("src/muFE_LJ.jl")

using .muFE_LJ
using StaticArrays
using ProgressBars
using Plots
using Trapz

#this script uses reduced units

params = LJ_params(1.0,1.0)
m = 1.0
γ = 1.0
T = 1.0
δt = 0.005
N = 256
ρ = 0.8

l = box_length(N, ρ)
rc = l/2.0
L = @SVector[l, l, l]
d = length(L)

println("Running simulation in a cubic box of reduced length $(l), at reduced temperature $(T) and reduced density $(ρ)")
println("The BAOAB integrators will have a timestep in reduced units of $(δt) and use the friction coefficient $(γ)")

sample_every = 100
λ = 0.5
λ_steps = 10

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

#equilibrate at λ = 1, which corresponds to having N standard LJ particles
BAOAB!(p, q, m, params, 1.0, L, LJ_pot_der , dλsoft_pot_dr, γ, T, δt, 50000)

##########################################################################
# UNCOMMENT FOR COMPUTATION OF RDF AT a λ ≠ 1

#equilibrate at λ = 0
#BAOAB!(p, q, m, params, 0.0, L, LJ_pot_der , dλsoft_pot_dr, γ, T, δt, 50000)

#bring system to desired λ
#λ_schedule = range(0, λ; length = λ_steps+1)
#for i in ProgressBar(eachindex(λ_schedule))
#    BAOAB!(p, q, m, params, λ_schedule[i], L, LJ_pot_der , dλsoft_pot_dr, γ, T, δt, 10000)
#end
###########################################################################

#define bins
n_bins = 100
hist = zeros(Float64, n_bins)
dr = rc/n_bins
BAOAB_rdf!(hist, dr, p, q, m, params, λ, L, LJ_pot_der , dλsoft_pot_dr, γ, T, δt, 20000, sample_every)

#normalize
g = zeros(Float64, n_bins)
rs = [(k-0.5)*dr for k in 1:n_bins]  # midpoint of each bin
n_samples = div(20000, sample_every)

for k in 1:n_bins
    r = rs[k]
    shell_volume = (4π/3) * ((r+ dr/2.0)^3 - (r- dr/2.0)^3)
    g[k] = 2*hist[k] / (N*ρ*shell_volume*n_samples)
end

#compute approximation to the potential energy per particle using the rdf
integrand = zeros(length(rs))
for i in 1:length(rs)
    integrand[i] = (rs[i]^2)*LJ_pot(rs[i], params)*g[i]
end
U_per_particle_no_corr = 2π*ρ*trapz(rs, integrand)
U_per_particle = U_per_particle_no_corr + tailcorr_U_LJ(ρ, N, rc, params)/N


println("The rdf predicts an approximate potential energy per particle of $(U_per_particle)")
println("Before correction it was $(U_per_particle_no_corr)")


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

vline!([rc], label="cutoff radius", linewidth=2, linestyle=:dash, color=:red)
hline!([1], label="expected asymptote", linewidth=1, linestyle=:dash, color=:gray)

savefig("rdf_plots/rdf_$(ρ)_$(T)_$(λ).png")

