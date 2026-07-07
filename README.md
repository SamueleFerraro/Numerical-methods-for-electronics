# Numerical Methods for Electronics

MATLAB projects developed for **Computational Modeling in Electronics and Biomathematics** (096659), A.Y. 2025/2026, Politecnico di Milano — Prof. Luca Dedè, Prof. Stefano Micheletti.

The repository contains two independent projects on numerical methods for advection-diffusion-reaction (ADR) and drift-diffusion problems in electronics:

- **[Project 1](#project-1--primal-vs-mixed-fem-for-a-1d-diffusion-reaction-problem)** — Primal vs. mixed finite element formulations for a 1D diffusion-reaction problem with a rapidly oscillating diffusion coefficient.
- **[Project 3](#project-3--non-isothermal-energy-transport-drift-diffusion-model)** — Non-isothermal (energy-transport) drift-diffusion simulation of a unipolar n⁺-n-n⁺ silicon structure, solved via Gummel iteration with Scharfetter-Gummel discretization.

## Repository structure

```
Numerical-methods-for-electronics/
├── Project1/
│   ├── MATLAB Codes/
│   │   ├── mainProblem.m              # driver: solves, plots, and runs convergence analysis
│   │   ├── gfem_primal_dir_hom.m      # primal FEM formulation
│   │   ├── gfem_mixed_dir_hom.m       # mixed (dual) FEM formulation
│   │   └── supportfunctions/          # course-provided, see Authorship & attribution
│   │       ├── H1_error.m             # H1 norm/seminorm error (Gauss quadrature, P1/P2)
│   │       ├── L2_error.m             # L2 norm error (Gauss quadrature, P0/P1/P2)
│   │       ├── convergence_rate.m     # estimates order p and constant C from a mesh sequence
│   │       └── setfonts.m             # plot styling helper
│   └── Report Project 1.pdf
└── Project3/
    ├── MatlabScriptsAndData/
    │   ├── data.m                     # physical constants, device geometry, mesh, bias sweep
    │   ├── unipolar.m                 # driver: Gummel loop over the bias sweep
    │   ├── continuity_n.m             # SG-discretized continuity equation for n
    │   ├── electron_current_sg.m      # non-isothermal SG electron current J_n
    │   ├── energy_theta.m             # energy-balance equation for the electron thermal voltage ϑ
    │   ├── log_mean.m                 # logarithmic mean (singularity-safe)
    │   ├── bern.m / bern_mike.m       # Bernoulli function B(x) — course-provided, see Authorship & attribution
    │   ├── BernoulliF.m               # alternative Bernoulli function implementation — course-provided
    │   ├── energy_balance_diagnostics.m  # global energy-conservation identity check
    │   ├── isothermal_baseline.mat            # precomputed isothermal solution (continuation init. guess)
    │   └── isothermal_baseline_Assignment.mat # precomputed isothermal baseline (comparison reference)
    ├── Plot/                          # pre-rendered result figures (I-V curve, bias sweep, mesh refinement)
    └── Samuele_Ferraro_project3_Report.pdf
```

## Project 1 — Primal vs. mixed FEM for a 1D diffusion-reaction problem

Solves the boundary value problem

```
-(μ(x) u'(x))' + σ u(x) = f(x),   x ∈ (0,1),   u(0) = u(1) = 0
```

with a rapidly oscillating diffusion coefficient `μ(x) = 1.5 + sin(30πx)` (15 full periods over the domain) and a manufactured exact solution, used to compare two P1 Galerkin FEM formulations:

- **Primal formulation** (`gfem_primal_dir_hom.m`): standard displacement-based FEM in `u`; the flux `J = -μ u'` is recovered in post-processing. The stiffness matrix is assembled using the **arithmetic mean** of `μ` on each element (trapezoidal quadrature).
- **Mixed formulation** (`gfem_mixed_dir_hom.m`): `u` and the flux `J` are both treated as independent unknowns (saddle-point/Schur-complement system), which builds in local mass conservation by construction. The flux-flux block is assembled from the **harmonic mean** of `μ` (via the trapezoidal rule applied to `1/μ`).

`mainProblem.m` runs both formulations, plots `u(x)` and `J(x)` against the exact solution, computes `H1(u)` and `L2(J)` errors, and performs a mesh-refinement study (`N = 3⁵, …, 3¹⁰`) to estimate the convergence order of each formulation via `convergence_rate.m`.

**Notable results** (see the report for the full derivation):
- For the specific validation mesh (`N+1 = 30` elements), `μ` sampled at the mesh nodes is aliased to a constant value, which makes the arithmetic and harmonic means coincide — primal and mixed formulations then produce numerically identical errors, a mesh-specific degeneracy rather than a general property of the two formulations.
- Over the non-aliased refinement sequence, both formulations converge at first order in `H1(u)` and `L2(J)`, matching the a priori estimates for P1/P0 mixed FEM.

**Run:** open `Project1/MATLAB Codes/mainProblem.m` in MATLAB (adds `supportfunctions` to the path automatically) and run it.

## Project 3 — Non-isothermal (energy-transport) drift-diffusion model

Simulates a 1D unipolar (electron-only) n⁺-n-n⁺ silicon structure (`L = 3 μm`, smoothed tanh doping profile, `N_D` from `10²² m⁻³` to `10²⁴ m⁻³`) under an applied bias sweep `V_a = 0 → 0.5 V`, comparing an **isothermal** drift-diffusion baseline against a **non-isothermal energy-transport** model that adds the electron temperature as a third unknown to capture hot-electron effects.

The model couples three equations, solved self-consistently at each bias point with **Gummel's decoupled iteration**:

1. **Poisson**: `-ε φ'' = q(N_D - n)`, linearized à la Gummel and updated with damping.
2. **Continuity**: `J_n' = 0`, with `J_n` given by the **non-isothermal Scharfetter-Gummel flux** (`continuity_n.m`, `electron_current_sg.m`), using the logarithmic mean of the electron thermal voltage `ϑ` between nodes (`log_mean.m`) and the Bernoulli function `B(x) = x/(eˣ-1)` (`bern_mike.m`, evaluated with a Taylor expansion near `x = 0` for numerical stability).
3. **Energy balance**: an SG-discretized equation for `ϑ(x) = k_B T_n(x)/q` (`energy_theta.m`), relaxing towards the lattice thermal voltage `V_T` with relaxation time `τ_E`.

`unipolar.m` drives the bias continuation loop (using each converged solution as the initial guess for the next bias point), and after convergence at each `V_a` it runs a battery of **residual/consistency diagnostics**: current conservation across the mesh, Poisson residual, energy-equation residual, a global energy-balance identity (`energy_balance_diagnostics.m`), a drift-velocity vs. saturation-velocity validity check, and a temperature positivity check.

**Key result:** for this device and bias range, the non-isothermal correction produces only a modest departure from the isothermal baseline — the relative current deviation stays below ≈0.5% and the maximum electron-temperature rise is ≈3.7% above the lattice temperature at `V_a = 0.5 V` — while still providing a principled, quantitative way to check whether the isothermal approximation is adequate for a given device/bias regime (see the report, §9, for the full discussion of model validity).

**Run:** open `Project3/MatlabScriptsAndData/unipolar.m` in MATLAB and run it (`data.m` is called automatically). The two `.mat` baseline files must be present in the same folder — they already are — since `unipolar.m` uses the isothermal solution both as the initial guess for the continuation loop and as the comparison baseline; it will raise an error if either file is missing. The script reproduces the bias-sweep figures, the I-V characteristic, and the mesh-refinement comparison saved in `Plot/`.

## Authorship & attribution

Not all files in this repository were written from scratch by the author. Files below were provided as course/lab material and reused as-is (or with minor adaptation); everything else in each project was implemented independently.

- **Project 1** — `MATLAB Codes/supportfunctions/` (`H1_error.m`, `L2_error.m` by A. Tonini; `setfonts.m` by R. Sacco; `convergence_rate.m`) is course-provided material, used to compute errors/convergence rates and to format plots. `mainProblem.m`, `gfem_primal_dir_hom.m`, and `gfem_mixed_dir_hom.m` are original work.
- **Project 3** — the Bernoulli function implementations (`bern.m`, `bern_mike.m`, `BernoulliF.m`) are course-provided material. `data.m`, `unipolar.m`, `continuity_n.m`, `electron_current_sg.m`, `energy_theta.m`, `log_mean.m`, and `energy_balance_diagnostics.m` are original work.

## Requirements

MATLAB, no additional toolboxes (only core functions: `spdiags`, sparse linear algebra, `arrayfun`/`integral` for quadrature).

## Reports

Each project folder includes the corresponding write-up with full derivations, figures, and discussion:
- [`Project1/Report Project 1.pdf`](Project1/Report%20Project%201.pdf)
- [`Project3/Samuele_Ferraro_project3_Report.pdf`](Project3/Samuele_Ferraro_project3_Report.pdf)

## Author

Samuele Ferraro
