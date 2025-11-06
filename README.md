# 🧪 Free energy alchemical calculation of chemical potential in a dense Lennard-Jones fluid

This project implements **alchemical particle insertion** in a dense Lennard-Jones (LJ) fluid to estimate the **chemical potential** (μ) using free-energy difference estimation methods.  
The approach gradually "turns on" the interactions of an additional particle through a coupling parameter **λ ∈ [0, 1]**, avoiding the low acceptance probability of direct insertion in crowded systems.

---

## 🔬 Physical Background

The system Hamiltonian is:

\[
H(λ) =
\sum_{i=1}^{N} \frac{p_i^2}{2m}
+ \sum_{i<j}^{N} V(r_i - r_j)
+ λ \sum_{i=1}^{N} V(r_0 - r_i)
\]

where \(V(r)\) is the Lennard-Jones potential, and \(λ\) controls the interaction strength of the additional particle at position \(r_0\).  
To improve numerical stability, the last term is replaced by a **soft-core** form:

\[
λ \sum_{i=1}^{N} 4ε(A^{-12} - A^{-6}), \quad
A = 0.5(1 - λ) + (r/σ)
\]

so that interactions are regular even as \(r \to 0\).

The free-energy difference between the fully interacting and non-interacting states is

\[
ΔF = F(λ=1) - F(λ=0)
\]

and the chemical potential follows from

\[
μ = ΔF + k_B T \ln\!\left[\frac{(N+1) Λ^d}{L^d}\right]
\]

with \(Λ = h / \sqrt{2π m k_B T}\) the thermal wavelength.

---

## ⚙️ Simulation Method

- **Ensemble:** canonical (**NVT**) using a thermostat.
- **Integrator:** velocity-Verlet (or BAOAB Langevin splitting).
- **Potential:** standard Lennard-Jones or soft-core modification. 
- **Available estimators:** Thermodynamic Integration (TI), Free Energy Perturbation (FEP), Bennett’s Acceptance Ratio (BAR), or WHAM/MBAR for post-processing.

---

