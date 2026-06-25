# Experiments

This repository contains the numerical experiments and supporting code for **Block-Filtering Dürr–Hoyer (BFDH)**, an analysis of the trade-off between the number of elements marked by the oracle and the resulting threshold drop in the [Dürr–Hoyer](https://arxiv.org/abs/quant-ph/9607014) quantum minimum-finding algorithm, which itself builds on [Grover's algorithm](https://arxiv.org/abs/quant-ph/9605043) and the iteration bounds derived by [Boyer, Brassard, Høyer, and Tapp](https://arxiv.org/abs/quant-ph/9605034). The accompanying paper is written and maintained separately; this repository exists to make the numerical claims in that paper independently reproducible and verifiable.

The central object of study is a [filtering fraction](https://en.wikipedia.org/wiki/Amplitude_amplification) $\alpha$, which controls how aggressively the oracle restricts the search space relative to the current threshold, and a [step budget](https://en.wikipedia.org/wiki/Big_O_notation) $S$, the fixed number of iterations over which $\alpha$ is held constant. The code in this repository searches for the value of $\alpha$ that maximizes the efficiency ratio $R_S(\alpha)$ for a given $S$, and verifies that this optimum shifts predictably as $S$ grows. The motivation, full derivation, and proofs are presented in the paper; this repository contains only the code, raw output, and figures that support those derivations numerically.

## What the code computes

The script evaluates a closed-form expression for $R_S(\alpha)$, the ratio of BFDH's efficiency to that of standard Dürr–Hoyer over a fixed number of steps $S$. For a single step, this reduces to the classical efficiency-versus-cost trade-off familiar from [amplitude amplification](https://en.wikipedia.org/wiki/Amplitude_amplification): a narrower filtering fraction produces a larger threshold drop per step but requires more Grover iterations, and vice versa. The script searches numerically for the optimal $\alpha^*(S)$ at each step budget, using a dense parameter sweep cross-checked against a [Newton–Raphson](https://en.wikipedia.org/wiki/Newton%27s_method) root search on the derivative, and traces how this optimum moves as the step budget increases.

Three figures are generated. The first overlays the efficiency curve $R_S(\alpha)$ for several step budgets on a single set of axes, marking the optimum on each curve. The second plots the optimal $\alpha^*(S)$ against the step budget $S$ on a logarithmic axis, showing convergence toward a structural upper bound. The third illustrates what happens when an $\alpha$ tuned for one step budget is extended beyond that budget, contrasting the threshold trajectory and cumulative cost against the standard algorithm.

## Repository structure

| Path | Contents |
|---|---|
| `data/alpha_star_results.csv` | Raw numerical output: $`\alpha^*(S)`$ and $`R_S(\alpha^*)`$ for each tested step budget, plus the comparison data behind the third figure. |
| `figures/fig1_RS_curves.pdf` | $`R_S(\alpha)`$ overlaid for several step budgets, with the optimum marked on each curve. |
| `figures/fig2_alpha_star_vs_S.pdf` | $`\alpha^*(S)`$ plotted against the step budget $`S`$ on a logarithmic axis, alongside the structural bound. |
| `figures/fig3_step_budget_comparison.pdf` | Threshold trajectory and cumulative cost when a step-budget-tuned $`\alpha`$ is extended past its intended budget. |


## Running the experiments

The script depends on [NumPy](https://numpy.org/), [Matplotlib](https://matplotlib.org/), and [SciPy](https://scipy.org/). With these installed, running the script from the repository root regenerates every figure and the results table from scratch:

```bash
pip install -r requirements.txt
python3 bfdh_experiments.py
```

Output is written to `figures/` and `data/`, overwriting the committed copies. This is intentional: the committed files represent one verified run, and re-running the script should reproduce them exactly, since every computation here is deterministic.

## Verifying the results independently

Each value in `data/alpha_star_results.csv` can be checked by hand by substituting the corresponding $\alpha^*(S)$ back into the closed-form expression for $R_S(\alpha)$ inside `bfdh_experiments.py` and confirming it sits at the peak of the curve for that step budget. The script itself performs this cross-check by comparing the result of the dense parameter sweep against the Newton–Raphson solution at each step budget, and the two are reported alongside each other in the console output when the script runs.
