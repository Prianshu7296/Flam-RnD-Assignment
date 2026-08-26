"""Independent validation using the recovered t ordering.

The assignment describes an L1 distance between uniformly sampled points.
Because the supplied points do not include t, we first recover t_est using
the fitted parameters, sort by t_est, interpolate the observed x(t), y(t)
onto a uniform t-grid, and compare that grid directly with the analytical
curve at the same t values.

This is deliberately separate from the fitting objective.
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = ROOT / "results" / "fit_result.txt"
DATA_PATH = ROOT / "data" / "xy_data.csv"


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
    y = 42.0 + t * np.sin(theta) + B * np.cos(theta)
    return x, y


def uniform_grid_validation(theta, M, X, x, y, n_grid=5000):
    c, s = np.cos(theta), np.sin(theta)

    # Recover the hidden parameter ordering from the fitted geometry.
    t_est = (x - X) * c + (y - 42.0) * s
    order = np.argsort(t_est)
    t_sorted = t_est[order]
    x_sorted = x[order]
    y_sorted = y[order]

    # Only compare over the actually observed domain.
    t_grid = np.linspace(
        max(6.0, t_sorted.min()),
        min(60.0, t_sorted.max()),
        n_grid,
    )

    # Reconstruct the observed curve at the same uniform t values.
    x_obs = np.interp(t_grid, t_sorted, x_sorted)
    y_obs = np.interp(t_grid, t_sorted, y_sorted)

    # Analytical prediction at those exact t values.
    x_pred, y_pred = curve_xy(theta, M, X, t_grid)

    point_l1 = np.abs(x_pred - x_obs) + np.abs(y_pred - y_obs)

    return {
        "t_grid": t_grid,
        "x_obs": x_obs,
        "y_obs": y_obs,
        "x_pred": x_pred,
        "y_pred": y_pred,
        "l1": point_l1,
        "mean_l1": float(point_l1.mean()),
        "max_l1": float(point_l1.max()),
        "p95_l1": float(np.percentile(point_l1, 95)),
        "t_min": float(t_sorted.min()),
        "t_max": float(t_sorted.max()),
    }


def main():
    res = load_result()
    df = pd.read_csv(DATA_PATH)
    x = df["x"].to_numpy(float)
    y = df["y"].to_numpy(float)

    out = uniform_grid_validation(
        res["theta_rad"], res["M"], res["X"], x, y
    )

    print(f"theta = {res['theta_deg']:.6f} deg")
    print(f"M     = {res['M']:.8f}")
    print(f"X     = {res['X']:.8f}")
    print(f"Recovered t range = [{out['t_min']:.6f}, {out['t_max']:.6f}]")
    print(f"Mean uniform-grid L1 = {out['mean_l1']:.8f}")
    print(f"Max  uniform-grid L1 = {out['max_l1']:.8f}")
    print(f"P95  uniform-grid L1 = {out['p95_l1']:.8f}")

    with open(ROOT / "results" / "validation.txt", "w") as f:
        f.write(f"mean_uniform_l1={out['mean_l1']:.10f}\n")
        f.write(f"max_uniform_l1={out['max_l1']:.10f}\n")
        f.write(f"p95_uniform_l1={out['p95_l1']:.10f}\n")
        f.write(f"t_min={out['t_min']:.10f}\n")
        f.write(f"t_max={out['t_max']:.10f}\n")


if __name__ == "__main__":
    main()
