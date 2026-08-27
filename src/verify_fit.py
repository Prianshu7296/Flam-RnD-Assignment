"""
Standalone verification script.

Loads the fitted parameters from results/fit_result.txt and independently
validates the final reported curve against the supplied unordered data.

The final assignment answer is reported to four significant decimal places
for theta/M/X as:

    theta = 30 degrees
    M     = 0.03
    X     = 55

Because the supplied points are unordered, an estimated t value is first
recovered for each observation using the inverse rotation. The observations
are then sorted by the recovered t value and linearly interpolated onto a
common uniform t-grid.

The fitted parametric curve is evaluated on the same grid and the
coordinate-wise L1 error is calculated as:

    |x_pred - x_ref| + |y_pred - y_ref|

Run:
    python src/verify_fit.py
"""

import numpy as np
import pandas as pd

RESULT_PATH = "results/fit_result.txt"
DATA_PATH = "data/xy_data.csv"
N_GRID = 5000


def load_result(path=RESULT_PATH):
    """Load the fitted parameters saved by fit_curve.py."""
    vals = {}

    with open(path) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            key, value = line.split("=", 1)
            vals[key] = float(value)

    required = {"theta_rad", "theta_deg", "M", "X"}
    missing = required - vals.keys()

    if missing:
        raise ValueError(
            f"Missing required values in {path}: {sorted(missing)}"
        )

    return vals


def curve_xy(theta, M, X, t):
    """Evaluate the original parametric curve."""
    B = np.exp(M * np.abs(t)) * np.sin(0.3 * t)

    x = t * np.cos(theta) - B * np.sin(theta) + X
    y = 42.0 + t * np.sin(theta) + B * np.cos(theta)

    return x, y


def recover_t(x, y, theta, X):
    """Recover the latent t coordinate using the inverse rotation."""
    return (
        (x - X) * np.cos(theta)
        + (y - 42.0) * np.sin(theta)
    )


def main():
    # Load the saved optimizer result.
    res = load_result()

    # Validate the rounded parameters used as the final assignment answer.
    # Rounding here makes the verification consistent with the parameters
    # reported in the README and Desmos equation.
    theta_deg = round(res["theta_deg"])
    theta = np.deg2rad(theta_deg)
    M = round(res["M"], 2)
    X = round(res["X"])

    # Load supplied data.
    df = pd.read_csv(DATA_PATH)

    if not {"x", "y"}.issubset(df.columns):
        raise ValueError("CSV must contain 'x' and 'y' columns.")

    x = df["x"].to_numpy(dtype=float)
    y = df["y"].to_numpy(dtype=float)

    if len(x) != 1500:
        raise ValueError(
            f"Expected 1500 observations, found {len(x)}."
        )

    # Recover latent t values from the unordered observations.
    t_est = recover_t(x, y, theta, X)

    # Sort observations according to recovered t.
    order = np.argsort(t_est)

    t_sorted = t_est[order]
    x_sorted = x[order]
    y_sorted = y[order]

    # Use the observed support inside the assignment range.
    t_min = float(t_sorted.min())
    t_max = float(t_sorted.max())

    if not (6.0 < t_min < t_max < 60.0):
        raise ValueError(
            f"Recovered t-range [{t_min}, {t_max}] "
            "is outside the assignment range 6 < t < 60."
        )

    # Uniform sampling grid over the observed support.
    t_grid = np.linspace(t_min, t_max, N_GRID)

    # Reconstruct the supplied/reference curve on the uniform grid.
    x_ref = np.interp(t_grid, t_sorted, x_sorted)
    y_ref = np.interp(t_grid, t_sorted, y_sorted)

    # Evaluate the final fitted curve on the same grid.
    x_pred, y_pred = curve_xy(theta, M, X, t_grid)

    # Coordinate-wise L1 error between corresponding points.
    l1 = (
        np.abs(x_pred - x_ref)
        + np.abs(y_pred - y_ref)
    )

    mean_l1 = float(np.mean(l1))
    max_l1 = float(np.max(l1))
    p95_l1 = float(np.percentile(l1, 95))

    print(f"theta = {theta_deg:.8f} deg")
    print(f"M     = {M:.10f}")
    print(f"X     = {X:.10f}")
    print(f"t_min = {t_min:.10f}")
    print(f"t_max = {t_max:.10f}")
    print(f"uniform grid points = {N_GRID}")
    print(f"Mean uniform-grid L1: {mean_l1:.10f}")
    print(f"Max  uniform-grid L1: {max_l1:.10f}")
    print(f"95th pct L1         : {p95_l1:.10f}")


if __name__ == "__main__":
    main()
