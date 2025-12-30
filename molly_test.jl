using Molly
using Molly.Potentials
using Molly.Observables

# Simulation parameters
N = 256
rho = 0.45
T = 1.2
dt = 0.005
nsteps = 10000

# Box length
L = (N/rho)^(1/3)
box = CubicBox(L)

# Particles
particles = ParticleArray(N)
set_positions!(particles, SimpleCubic(L, N))
set_velocities!(particles, MaxwellBoltzmann(T))

# LJ 12-6 potential, truncated at half-box (minimum image)
lj = LennardJones(1.0, 1.0; cutoff=L/2)
forcefield = ForceField(lj)

# Integrator: BAOAB Langevin
integrator = LangevinBAOAB(dt, 1.0, T) # friction γ = 1.0

# Observables
energy_obs = EnergyObserver(forcefield)
energy_data = []

# Simulation loop
for step in 1:nsteps
    integrate!(particles, forcefield, integrator, box)
    push!(energy_data, total_energy(particles, forcefield))
end

println("Average energy = ", mean(energy_data))
