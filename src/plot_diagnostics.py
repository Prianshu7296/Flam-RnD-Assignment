"""Generate the two R&D diagnostics missing from the original repo."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from verify_fit import load_result, uniform_grid_validation

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "data" / "xy_data.csv"


def main():
    res = load_result()
    df = pd.read_csv(DATA_PATH)
    x, y = df["x"].to_numpy(float), df["y"].to_numpy(float)

    theta, M, X = res["theta_rad"], res["M"], res["X"]
    c, s = np.cos(theta), np.sin(theta)

    t_est = (x - X) * c + (y - 42.0) * s
    B_est = -(x - X) * s + (y - 42.0) * c

    order = np.argsort(t_est)
    t = t_est[order]
    B = B_est[order]
    B_model = np.exp(M * np.abs(t)) * np.sin(0.3 * t)

    # Diagnostic 1: the de-rotated signal should follow the one-dimensional model.
    plt.figure(figsize=(9, 5))
    plt.scatter(t, B, s=6, alpha=0.35, label="de-rotated data")
    plt.plot(t, B_model, linewidth=2, label="model")
    plt.xlabel("recovered t")
    plt.ylabel("B(t)")
    plt.title("De-rotation check: recovered oscillation vs model")
    plt.legend()
    plt.tight_layout()
    plt.savefig(ROOT / "results" / "derotation_diagnostic.png", dpi=180)
    plt.close()

    # Diagnostic 2: exact metric residual on the same uniform t-grid.
    val = uniform_grid_validation(theta, M, X, x, y, n_grid=5000)
    plt.figure(figsize=(9, 5))
    plt.plot(val["t_grid"], val["l1"], linewidth=1.5)
    plt.axhline(val["mean_l1"], linestyle="--",
                label=f"mean = {val['mean_l1']:.2e}")
    plt.xlabel("t")
    plt.ylabel("L1 error")
    plt.title("Uniform-grid L1 residual")
    plt.legend()
    plt.tight_layout()
    plt.savefig(ROOT / "results" / "residual_vs_t.png", dpi=180)
    plt.close()

    print("Saved diagnostic plots.")


if __name__ == "__main__":
    main()
