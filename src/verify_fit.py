"""
Standalone verification script.

Loads the fitted parameters from results/fit_result.txt, re-derives the
curve, and reports the L1 distance metric described in the assignment's
assessment criteria: for uniformly sampled points along the fitted curve,
compute the (nearest-neighbour) L1 distance to the real data points.

Run:
    python src/verify_fit.py
"""
import numpy as np
import pandas as pd

RESULT_PATH = "results/fit_result.txt"
DATA_PATH = "data/xy_data.csv"


def load_result(path=RESULT_PATH):
    vals = {}
    with open(path) as f:
        for line in f:
            k, v = line.strip().split("=")
            vals[k] = float(v)
    return vals


def curve_xy(theta, M, X, t):
    B = np.exp(M * np.abs(t)) * np.sin(0.3 * t)
    x = t * np.cos(theta) - B * np.sin(theta) + X
    y = 42 + t * np.sin(theta) + B * np.cos(theta)
    return x, y


def main():
    res = load_result()
    theta, M, X = res["theta_rad"], res["M"], res["X"]

    df = pd.read_csv(DATA_PATH)
    x, y = df["x"].to_numpy(), df["y"].to_numpy()

    t_grid = np.linspace(6, 60, 5000)
    x_fit, y_fit = curve_xy(theta, M, X, t_grid)

    l1s = np.array([
        np.min(np.abs(x_fit - xi) + np.abs(y_fit - yi))
        for xi, yi in zip(x, y)
    ])

    print(f"theta = {res['theta_deg']:.4f} deg")
    print(f"M     = {M:.6f}")
    print(f"X     = {X:.6f}")
    print(f"Mean L1 distance (data -> fitted curve): {l1s.mean():.6f}")
    print(f"Max  L1 distance (data -> fitted curve): {l1s.max():.6f}")
    print(f"95th pct L1 distance                   : {np.percentile(l1s, 95):.6f}")


if __name__ == "__main__":
    main()
