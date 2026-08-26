"""
Generates results/fit_plot.png : overlays the fitted curve on top of the
raw data points, for visual sanity-checking.

Run:
    python src/plot_fit.py
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from verify_fit import load_result, curve_xy

DATA_PATH = "data/xy_data.csv"
OUT_PATH = "results/fit_plot.png"


def main():
    res = load_result()
    theta, M, X = res["theta_rad"], res["M"], res["X"]

    df = pd.read_csv(DATA_PATH)
    x, y = df["x"].to_numpy(), df["y"].to_numpy()

    t_grid = np.linspace(6, 60, 3000)
    x_fit, y_fit = curve_xy(theta, M, X, t_grid)

    plt.figure(figsize=(8, 6))
    plt.scatter(x, y, s=6, alpha=0.35, color="#1f77b4", label="data points (xy_data.csv)")
    plt.plot(x_fit, y_fit, color="#d62728", linewidth=2, label="fitted curve")
    plt.title(f"Fitted curve  |  theta={res['theta_deg']:.2f} deg, M={M:.4f}, X={X:.2f}")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.gca().set_aspect("equal", adjustable="datalim")
    plt.tight_layout()
    plt.savefig(OUT_PATH, dpi=150)
    print(f"Saved plot to {OUT_PATH}")


if __name__ == "__main__":
    main()
