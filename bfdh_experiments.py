"""
Numerical experiments supporting the BFDH paper.
Generates three figures (vector PDF, ready for LaTeX \\includegraphics)
and a raw CSV of results so every number can be checked independently.

Run: python3 bfdh_experiments.py
Output: figures/fig1_RS_curves.pdf
        figures/fig2_alpha_star_vs_S.pdf
        figures/fig3_step_budget_comparison.pdf
        data/alpha_star_results.csv
"""

import numpy as np
import matplotlib.pyplot as plt
import csv
import os

os.makedirs("figures", exist_ok=True)
os.makedirs("data", exist_ok=True)

# Plot style: grayscale, serif font, suitable for an academic paper.
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 9,
    "axes.labelsize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": ":",
})


# ---------------------------------------------------------------------------
# Core function: R_S(alpha), the closed-form efficiency ratio from the paper.
# ---------------------------------------------------------------------------
def R_S(alpha, S):
    """Efficiency ratio of BFDH versus standard DH over a fixed step budget S.
    alpha may be a scalar or a numpy array; S must be a scalar (int)."""
    alpha = np.asarray(alpha, dtype=np.float64)
    k = np.arange(S)
    drop_ratio = (1 - (alpha / 2) ** S) / (1 - (1 / 2) ** S)

    sum_baseline = np.sum(np.sqrt(2.0) ** k)
    r = np.sqrt(2.0 / alpha)
    if alpha.ndim == 0:
        sum_bfdh = alpha ** (-0.5) * np.sum(r ** k)
    else:
        sum_bfdh = alpha ** (-0.5) * np.sum(r[:, None] ** k[None, :], axis=1)

    return drop_ratio * (sum_baseline / sum_bfdh)


def find_optimal_alpha(S, lo=1e-3, hi=1.999999, n=20_000):
    """Vectorized numerical sweep to locate alpha*(S)."""
    alphas = np.linspace(lo, hi, n)
    Rs = R_S(alphas, S)
    idx = np.argmax(Rs)
    return alphas[idx], Rs[idx]


# ---------------------------------------------------------------------------
# FIGURE 1: R_S(alpha) curves for several values of S, with each peak marked.
# Demonstrates visually that each S has a single maximum, and that the
# location of that maximum shifts to the right as S grows.
# ---------------------------------------------------------------------------
def make_figure1():
    fig, ax = plt.subplots(figsize=(3.4, 2.6))

    S_list = [1, 2, 3, 5, 10]
    alphas_plot = np.linspace(0.05, 1.999, 2000)
    linestyles = ["-", "--", "-.", ":", (0, (3, 1, 1, 1))]
    grays = ["0.05", "0.25", "0.40", "0.55", "0.70"]

    for S, ls, gray in zip(S_list, linestyles, grays):
        Rs = R_S(alphas_plot, S)
        ax.plot(alphas_plot, Rs, label=f"$S={S}$", linewidth=1.2,
                linestyle=ls, color=gray)

        a_star, R_star = find_optimal_alpha(S)
        ax.plot(a_star, R_star, "o", color="black", markersize=3, zorder=5)

    ax.axhline(1.0, color="gray", linewidth=0.6, linestyle="--")
    ax.set_xlabel(r"Filtering fraction $\alpha$")
    ax.set_ylabel(r"$R_S(\alpha)$")
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 3)
    ax.legend(loc="upper left", ncol=2, frameon=False, fontsize=7)
    fig.tight_layout(pad=0.4)
    fig.savefig("figures/fig1_RS_curves.pdf")
    plt.close(fig)
    print("Figure 1 done: figures/fig1_RS_curves.pdf")


# ---------------------------------------------------------------------------
# FIGURE 2: alpha*(S) vs S, compared against the alpha=2 limiting line.
# Demonstrates visually the monotonically increasing trend toward the
# structural bound.
# ---------------------------------------------------------------------------
def make_figure2():
    S_values = [1, 2, 3, 4, 5, 6, 7, 8, 10, 15, 20, 30, 50, 70, 100, 150, 200, 300, 500]
    alpha_stars = []
    for S in S_values:
        a_star, _ = find_optimal_alpha(S)
        alpha_stars.append(a_star)

    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    ax.plot(S_values, alpha_stars, "o-", color="black", markersize=3, linewidth=1.0,
            label=r"$\alpha^*(S)$ (numerical)")
    ax.axhline(2.0, color="gray", linewidth=0.8, linestyle="--", label=r"Bound $\alpha=2$")
    ax.axhline(2/3, color="gray", linewidth=0.6, linestyle=":", label=r"$\alpha^*(1)=2/3$")

    ax.set_xscale("log")
    ax.set_xlabel(r"Step budget $S$ (log scale)")
    ax.set_ylabel(r"$\alpha^*(S)$")
    ax.set_ylim(0.5, 2.05)
    ax.legend(loc="lower right", frameon=False, fontsize=7)
    fig.tight_layout(pad=0.4)
    fig.savefig("figures/fig2_alpha_star_vs_S.pdf")
    plt.close(fig)
    print("Figure 2 done: figures/fig2_alpha_star_vs_S.pdf")
    return S_values, alpha_stars


# ---------------------------------------------------------------------------
# FIGURE 3: behavior of a step-budget-tuned alpha when extended past its
# intended budget. Case study: T0=2^50, S=50.
# ---------------------------------------------------------------------------
def simulate(alpha, S, T0, N):
    T = float(T0)
    total_cost = 0.0
    history_T = [T]
    for k in range(S):
        M = alpha * T
        cost = (np.pi / 4) * np.sqrt(N / M)
        T = M / 2
        total_cost += cost
        history_T.append(T)
    return T, total_cost, history_T


def continue_to_target(T_start, target, N, alpha_continue=1.0):
    T = T_start
    extra_cost = 0.0
    extra_steps = 0
    while T > target and extra_steps < 5000:
        M = alpha_continue * T
        cost = (np.pi / 4) * np.sqrt(N / M)
        T = M / 2
        extra_cost += cost
        extra_steps += 1
    return T, extra_cost, extra_steps


def make_figure3():
    T0 = 2 ** 50
    N = 2 * T0
    S = 50

    T_std, C_std, hist_std = simulate(1.0, S, T0, N)
    T_b, C_b, hist_b = simulate(1.93, S, T0, N)

    # Continue alpha=1.93 past the original S=50 budget until T actually
    # reaches 1, switching to alpha=1.0 for the continuation.
    T_final_b, extra_cost_b, extra_steps_b = continue_to_target(T_b, 1.0, N)
    total_cost_b_completed = C_b + extra_cost_b

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6))

    # Left panel: T_k vs k, fixed S=50.
    ax = axes[0]
    ax.plot(range(len(hist_std)), np.array(hist_std) / T0, "-", color="black",
             linewidth=1.1, label=r"$\alpha=1.0$")
    ax.plot(range(len(hist_b)), np.array(hist_b) / T0, "--", color="black",
             linewidth=1.1, label=r"$\alpha=1.93$")
    ax.set_yscale("log")
    ax.set_xlabel(r"Step $k$")
    ax.set_ylabel(r"$T_k/T_0$ (log scale)")
    ax.set_title(r"(a) Fixed $S=50$", fontsize=8)
    ax.legend(frameon=False, fontsize=7)

    # Right panel: cumulative cost, fixed S=50 vs continuation past budget.
    ax = axes[1]
    labels = [r"$\alpha=1.0$" + "\n($S=50$)",
              r"$\alpha=1.93$" + "\n($S=50$)",
              r"$\alpha=1.93$" + "\n(continued past"+"\nbudget to completion)"]
    values = [C_std, C_b, total_cost_b_completed]
    bars = ax.bar(labels, values, color=["0.3", "0.6", "0.6"],
                   edgecolor="black", linewidth=0.6,
                   hatch=["", "", "//"])
    ax.set_yscale("log")
    ax.set_ylabel("Total Grover iterations (log scale)")
    ax.set_title(r"(b) Total cost", fontsize=8)
    for bar, v in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, v*1.3, f"{v:.2e}",
                 ha="center", va="bottom", fontsize=6.5, rotation=0)

    fig.tight_layout(pad=0.5)
    fig.savefig("figures/fig3_step_budget_comparison.pdf")
    plt.close(fig)
    print("Figure 3 done: figures/fig3_step_budget_comparison.pdf")

    return {
        "T0": T0, "N": N, "S": S,
        "T_std_final": T_std, "C_std": C_std,
        "T_b_final": T_b, "C_b": C_b,
        "extra_steps_b": extra_steps_b, "extra_cost_b": extra_cost_b,
        "total_cost_b_completed": total_cost_b_completed,
        "R_50": (T0 - T_b) / C_b / ((T0 - T_std) / C_std),
    }


# ---------------------------------------------------------------------------
# Save all raw data to CSV for independent verification.
# ---------------------------------------------------------------------------
def save_csv(S_values, alpha_stars, figure3_results):
    with open("data/alpha_star_results.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["=== alpha*(S) table (Figure 2) ==="])
        writer.writerow(["S", "alpha_star", "R_S(alpha_star)"])
        for S, a in zip(S_values, alpha_stars):
            writer.writerow([S, f"{a:.6f}", f"{R_S(a, S):.6f}"])

        writer.writerow([])
        writer.writerow(["=== Experiment 3: alpha extended past its step budget (Figure 3) ==="])
        for k, v in figure3_results.items():
            writer.writerow([k, v])
    print("CSV saved: data/alpha_star_results.csv")


if __name__ == "__main__":
    make_figure1()
    S_values, alpha_stars = make_figure2()
    figure3_results = make_figure3()
    save_csv(S_values, alpha_stars, figure3_results)

    print("\n--- Quick summary for manual verification ---")
    for S in [1, 2, 3, 5, 10, 50]:
        a, R = find_optimal_alpha(S)
        print(f"S={S:4d}  alpha*={a:.4f}  R_S(alpha*)={R:.4f}")
