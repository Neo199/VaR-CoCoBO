# VaR-CoCoBO: Variational Bayesian Optimisation for Constrained Combinatorial Problems

Repository for the simulations and methods presented in:

> **Variational Bayesian Optimisation for Constrained Combinatorial Problems**
> Niyati Seth, Michael Fop — September 2025

VaR-CoCoBO (Variational pRobabilistic **Co**nstrained **Co**mbinatorial **B**ayesian **O**ptimisation) extends the [VaR-CBO framework](https://github.com/Neo199/VaR-cBO) to settings where candidate solutions must satisfy hard feasibility constraints. Rather than handling constraints indirectly via penalty terms, VaR-CoCoBO integrates them explicitly into both the surrogate model and the acquisition function, ensuring all evaluated configurations are feasible by construction.

---

## Overview

Many real-world combinatorial optimisation problems — in facility location, scheduling, resource allocation, and molecular design — require solutions to satisfy structural, logical, or budget constraints. Ignoring these constraints wastes expensive objective evaluations on infeasible configurations.

VaR-CoCoBO addresses this by:

- **Inheriting the sparse variational Bayes surrogate** from VaR-CBO (spike-and-slab prior with mean-field variational inference), providing fast and scalable surrogate updates within the BO loop.
- **Extending probabilistic reparameterisation** to constrained domains, relaxing binary decisions to continuous Bernoulli probabilities optimised via L-BFGS-B.
- **Introducing a Gumbel-based constraint handling mechanism** that perturbs the relaxed probabilities with Gumbel noise to produce a stochastic ranking, then greedily constructs feasible binary solutions that satisfy all hard constraints. Feasibility is enforced during discretisation, not via penalty terms in the relaxed objective.

---

## Method

The core acquisition step proceeds as follows:

1. Fit the sparse VB surrogate to current observations.
2. Optimise the probabilistically reparameterised acquisition function over `θ ∈ [0, 1]^p` using L-BFGS-B.
3. Draw independent Gumbel noise `gⱼ ~ Gumbel(0,1)` for each variable.
4. Compute scores `sⱼ = θ̂ⱼ + gⱼ + λϕⱼ`, where `ϕⱼ` encodes problem-specific structural information (e.g. marginal coverage, contamination levels, item profits).
5. Rank variables by score and greedily select a feasible binary configuration satisfying all constraints `cₖ(x) ≤ 0`.
6. Evaluate the true objective and update the dataset.

This is summarised in **Algorithm 1** of the paper. The main implementation lives in `varcoco-bo.R`.

---

## Benchmarks and Applications

Methods are evaluated on:

- **0/1 Knapsack Problem** (`p = 24` items, instance P08): constrained maximisation under a weight budget, with feasibility guaranteed by construction across all runs.
- **Contamination Control** (`p = 25` stages): minimising prevention costs subject to probabilistic safety constraints on bacterial prevalence at each supply chain stage.
- **Facility Location — San Francisco** (`p = 16` candidate sites, 205 demand points):
  - *Location Set Covering Problem (LSCP)*: minimise number of facilities to achieve full demand coverage within 5 km.
  - *Maximal Covering Location Problem (MCLP)*: maximise weighted demand coverage with a budget of 4 facilities.

On the LSCP, VaR-CoCoBO recovers the global optimum (8 facilities, 205/205 demand points covered) with only 255 objective evaluations, compared to 3,637 for a genetic algorithm and 1.5 seconds for an exact GLPK solver. On the MCLP, VaR-CoCoBO achieves a mean coverage of 90.06% (SD 1.54%) across 20 runs, reaching the GLPK optimum of 91.63% in 7 of 20 instances, while requiring far fewer objective evaluations than the GA.

---

## Installation

### Required: `sparsevb` (GitHub version)

A modified version of the `sparsevb` R package is required. Install directly from GitHub:

```r
# install.packages("remotes")
remotes::install_github("<your-github-org>/sparsevb")
```

> **Note:** The CRAN version of `sparsevb` will not work. Changes were made to integrate the package into the VaR-CBO/VaR-CoCoBO framework.

### Required: Stan

Stan is required to run BOCS and PRCBO comparison methods (via `rstanarm`). Follow the installation instructions at [mc-stan.org](https://mc-stan.org/users/interfaces/rstan).

### Required: GLPK (for facility location comparisons)

The `Rglpk` package is used to compute exact solutions for the LSCP and MCLP benchmarks:

```r
install.packages("Rglpk")
```

This requires the GLPK library to be installed on your system. See [gnu.org/software/glpk](https://www.gnu.org/software/glpk/) for instructions.

### All other dependencies

```r
install.packages(c("GA", "optimx", "rstanarm", "maxcovr"))
```

---

## Repository Structure

```
.
├── minimisation/        # VaR-CoCoBO and comparison methods for minimisation problems
├── maximisation/        # Variants adapted for maximisation problems
├── facility-location/   # LSCP and MCLP experiments (San Francisco dataset)
├── varcoco-bo.R         # Main VaR-CoCoBO implementation (primary contribution)
└── README.md
```

Minimisation and maximisation variants are separated into distinct folders for ease of use. The facility location experiments use data sourced from the `maxcovr` R package and the associated [GitHub repository](https://github.com/njtierney/maxcovr), with coverage matrices constructed using ArcGIS Network Analyst.

---

## Usage

To run the benchmark examples, navigate to the relevant folder and source the corresponding script. Problem setups follow the experimental configurations described in the paper.

To apply VaR-CoCoBO to a custom constrained problem, the main function in `varcoco-bo.R` requires:

- A black-box objective function `f` over binary inputs `x ∈ {0, 1}^p`
- A set of constraint functions `cₖ(x) ≤ 0` for `k = 1, ..., K`
- A problem-specific auxiliary score `ϕⱼ` for each variable (encodes structural information for the Gumbel ranking step)
- An initial dataset of evaluated points and a total evaluation budget `N_max`

The function returns the best feasible solution observed within the budget.

---

## Relationship to VaR-CBO

This repository extends [VaR-CBO](https://github.com/<your-github-org>/VaR-cBO) (Seth and Fop, 2023). The surrogate model and probabilistic reparameterisation are inherited directly from that work; the contribution here is the Gumbel-based constraint handling mechanism and its application to constrained combinatorial problems. Both repositories share the modified `sparsevb` dependency.

---

## Citation

If you use this code in your work, please cite both papers:

```bibtex
@article{seth2025varcocbo,
  title     = {Variational {B}ayesian Optimisation for Constrained Combinatorial Problems},
  author    = {Seth, Niyati and Fop, Michael},
  year      = {2025}
}

@article{seth2023varcbo,
  title     = {Scalable Variational {B}ayesian Optimisation for Combinatorial Problems},
  author    = {Seth, Niyati and Fop, Michael},
  year      = {2023}
}
```

---

## Contact

For questions or issues, please open a GitHub issue or contact the authors.
