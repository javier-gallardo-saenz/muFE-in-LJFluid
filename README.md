# 🧪 Free energy alchemical calculation of chemical potential in a dense Lennard-Jones fluid

This project implements **alchemical particle insertion** in a dense Lennard-Jones (LJ) fluid to estimate the **chemical potential** (μ) using free-energy difference estimation methods.  
The approach gradually "turns on" the interactions of an additional particle through a coupling parameter **λ ∈ [0, 1]**, avoiding the low acceptance probability of direct insertion in crowded systems.

---

## 🔬 Physical Background

The system Hamiltonian is:

H(λ) = Σᵢ (pᵢ² / 2m) + Σ_{i<j} V(rᵢ − rⱼ) + λ Σᵢ V(r₀ − rᵢ)

where V(r) is the Lennard-Jones potential, and λ controls the interaction strength of the additional particle at position r₀.  

To improve numerical stability, the last term is replaced by a **soft-core** form:

λ Σᵢ 4ε (A⁻¹² − A⁻⁶), where A = 0.5 (1 − λ) + (r / σ)

so that interactions are regular even as r → 0.

The free-energy difference between the fully interacting and non-interacting states is:

ΔF = F(λ = 1) − F(λ = 0)

and the chemical potential follows from:

μ = ΔF + k_B T ln [ (N + 1) Λᵈ / Lᵈ ]

with Λ = h / √(2 π m k_B T) the thermal wavelength.

---

## ⚙️ Simulation and Estimation Methods

- **Ensemble:** Canonical (**NVT**) with Langevin thermostat  
- **Integrator:** BAOAB
- **Potential:** Lennard-Jones with optional (but necessary for particle insertion simulations to work) soft-core modification  
- **Free-energy estimators:**
  - Thermodynamic Integration (TI)
  - Free Energy Perturbation (FEP)
  - Bennett Acceptance Ratio (BAR)
  - Multistate Bennett Acceptance Ratio (MBAR)

---


## 📁 Project Structure

```markdown
.
├── src/
│   └── muFE_LJ.jl
│   └── [muFE_LJ_implementation].jl
├── ti.jl
├── fep.jl
├── bar_mbar.jl
├── widom.jl
├── rdf.jl
├── perf_checks.jl
├── [method]_data/
├── [method]_plots/
```


### `muFE_LJ` module

All core functionality is implemented in the Julia module **`muFE_LJ`**, located in `src/`.  
The module exports all user-facing routines required for free-energy estimation.

All top-level scripts import this module and are intended to be run directly.

---


## ▶️ How to Run

### 1️⃣ Requirements and installation

- Julia **≥ 1.9** (recommended)

1. Clone the repository

```bash
git clone https://github.com/yourname/muFE_LJ.git
cd muFE_LJ
```

2. Activate project environment and install all dependencies:

```bash
julia
] activate .
] instantiate
```

### 2️⃣ Running a free-energy estimator

Each script is self-contained and intended to be user-editable.
Typically, the user only needs to modify:

thermodynamic parameters (temperature, density, number of particles)

paths to the provided MD data

estimator-specific parameters (λ grid, block sizes, etc.)

#### Thermodynamic Integration
julia ti.jl

#### Free Energy Perturbation
julia fep.jl

#### BAR / MBAR
julia bar_mbar.jl
Additionally, if one wants to run only BAR it is enough to modify the name of the function in line 42 to BAR

#### Widom insertion
julia widom.jl

---

## 📊 Analysis and Validation

rdf.jl computes radial distribution functions for structural validation

perf_checks.jl computes the temperature evolution and internal energy of the system

---

## 📌 Notes

Pre-generated molecular dynamics trajectories are provided, allowing direct execution of the free-energy estimators without rerunning MD.

The scripts are designed to be readable and easily modifiable for testing alternative parameters or methods.

