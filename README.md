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
- **Integrator:** BAOAB Langevin
- **Potential:** Lennard-Jones with optional (but necessary for particle insertion simulations to work) soft-core modification  
- **Free-energy estimators:**
  - Thermodynamic Integration (TI)
  - Free Energy Perturbation (FEP)
  - Bennett Acceptance Ratio (BAR)
  - Multistate Bennett Acceptance Ratio (MBAR)

---


## 📁 Project Structure

.
├── src/
│ └── muFE_LJ.jl # Main Julia module
| └── [muFE_LJ_implementation].jl # muFE_LJ funtions implementation
│
├── ti.jl # Thermodynamic Integration
├── fep.jl # Free Energy Perturbation
├── bar_mbar.jl # BAR and MBAR estimators
├── widom.jl # Widom insertion (reference)
├── rdf.jl # Radial distribution function
├── perf_checks.jl # Performance benchmarks
│
├── [method]_data/ # Pre-generated MD data (provided, method in [TI, FEP, MBAR], where we use MBAR's data for BAR)


---


