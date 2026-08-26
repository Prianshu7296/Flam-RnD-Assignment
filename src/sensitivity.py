"""Small local sensitivity check around the recovered solution."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from verify_fit import load_result, uniform_grid_validation

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "data" / "xy_data.csv"


def metric(p, x, y):
    return uniform_grid_validation(*p, x, y, n_grid=3000)["mean_l1"]


def main():
    res = load_result()
    x, y = (
        pd.read_csv(DATA_PATH)["x"].to_numpy(float),
        pd.read_csv(DATA_PATH)["y"].to_numpy(float),
    )

    base = np.array([res["theta_rad"], res["M"], res["X"]])
    labels = [
        "theta -1°", "theta +1°",
        "M -0.001", "M +0.001",
        "X -1", "X +1",
    ]
    perturbations = [
        [-np.deg2rad(1), 0, 0], [np.deg2rad(1), 0, 0],
        [0, -0.001, 0], [0, 0.001, 0],
        [0, 0, -1], [0, 0, 1],
    ]
    values = [metric(base + np.array(d), x, y) for d in perturbations]

    plt.figure(figsize=(9, 5))
    plt.bar(labels, values)
    plt.yscale("log")
    plt.ylabel("mean uniform-grid L1")
    plt.title("Local parameter sensitivity")
    plt.xticks(rotation=25, ha="right")
    plt.tight_layout()
    plt.savefig(ROOT / "results" / "sensitivity.png", dpi=180)
    plt.close()

    with open(ROOT / "results" / "sensitivity.txt", "w") as f:
        for label, value in zip(labels, values):
            f.write(f"{label}={value:.10g}\n")

    print("Saved sensitivity results.")


if __name__ == "__main__":
    main()
