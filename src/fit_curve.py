"""
Fit unknown parameters (theta, M, X) of the parametric curve:

    x(t) = t*cos(theta) - e^(M*|t|) * sin(0.3t) * sin(theta) + X
    y(t) = 42 + t*sin(theta) + e^(M*|t|) * sin(0.3t) * cos(theta)

given a cloud of (x, y) points sampled at unknown t in (6, 60).

KEY INSIGHT (de-rotation trick)
--------------------------------
Group the "t" part and the "wiggle" part:

    A(t) = t
    B(t) = e^(M*|t|) * sin(0.3t)

Then the equations are exactly a 2D rotation by theta plus a translation:

    [ x - X ]   [ cos(theta)  -sin(theta) ] [ A ]
    [ y - 42] = [ sin(theta)   cos(theta) ] [ B ]

So for ANY candidate (theta, X), we can invert the rotation for every data
point to recover an ESTIMATE of t and of B, without needing correspondence
or ordering of points:

    t_est = (x - X)*cos(theta) + (y - 42)*sin(theta)
    B_est = -(x - X)*sin(theta) + (y - 42)*cos(theta)

If (theta, M, X) are correct, then for every point:

    B_est  ==  exp(M * |t_est|) * sin(0.3 * t_est)

This turns a hard "unknown correspondence" curve-fitting problem into a
plain nonlinear least-squares problem in 3 unknowns (theta, M, X), solved
with scipy.optimize.least_squares.
"""
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

DATA_PATH = "data/xy_data.csv"


def load_data(path=DATA_PATH):
    df = pd.read_csv(path)
    return df["x"].to_numpy(), df["y"].to_numpy()


def residuals(params, x, y):
    theta, M, X = params
    ct, st = np.cos(theta), np.sin(theta)

    t_est = (x - X) * ct + (y - 42.0) * st
    B_est = -(x - X) * st + (y - 42.0) * ct

    B_model = np.exp(M * np.abs(t_est)) * np.sin(0.3 * t_est)

    r_wiggle = B_est - B_model

    lo, hi = 6.0, 60.0
    pen_lo = np.clip(lo - t_est, 0, None)
    pen_hi = np.clip(t_est - hi, 0, None)

    return np.concatenate([r_wiggle, 5.0 * pen_lo, 5.0 * pen_hi])


def global_search(x, y, n_starts=60, seed=42):
    rng = np.random.default_rng(seed)
    theta_bounds = (np.deg2rad(0.01), np.deg2rad(49.99))
    M_bounds = (-0.0499, 0.0499)
    X_bounds = (0.01, 99.99)

    best = None
    for _ in range(n_starts):
        p0 = np.array([
            rng.uniform(*theta_bounds),
            rng.uniform(*M_bounds),
            rng.uniform(*X_bounds),
        ])
        try:
            sol = least_squares(
                residuals, p0, args=(x, y),
                bounds=([theta_bounds[0], M_bounds[0], X_bounds[0]],
                        [theta_bounds[1], M_bounds[1], X_bounds[1]]),
                method="trf", max_nfev=5000,
            )
        except Exception:
            continue
        if best is None or sol.cost < best.cost:
            best = sol
    return best


def main():
    x, y = load_data()
    print(f"Loaded {len(x)} points")

    best = global_search(x, y, n_starts=60, seed=42)
    theta, M, X = best.x
    theta_deg = np.rad2deg(theta)

    print("\n=== Best fit ===")
    print(f"theta = {theta:.6f} rad = {theta_deg:.4f} deg")
    print(f"M     = {M:.6f}")
    print(f"X     = {X:.6f}")
    print(f"final cost (sum sq resid) = {best.cost:.6f}")

    t_grid = np.linspace(6, 60, 2000)
    ct, st = np.cos(theta), np.sin(theta)
    B = np.exp(M * np.abs(t_grid)) * np.sin(0.3 * t_grid)
    x_fit = t_grid * ct - B * st + X
    y_fit = 42 + t_grid * st + B * ct

    l1s = []
    for xi, yi in zip(x, y):
        d = np.abs(x_fit - xi) + np.abs(y_fit - yi)
        l1s.append(d.min())
    l1s = np.array(l1s)
    print(f"\nMean L1 point-to-curve distance : {l1s.mean():.6f}")
    print(f"Max  L1 point-to-curve distance : {l1s.max():.6f}")

    with open("results/fit_result.txt", "w") as f:
        f.write(f"theta_rad={theta:.8f}\n")
        f.write(f"theta_deg={theta_deg:.6f}\n")
        f.write(f"M={M:.8f}\n")
        f.write(f"X={X:.8f}\n")
        f.write(f"mean_L1={l1s.mean():.8f}\n")
        f.write(f"max_L1={l1s.max():.8f}\n")

    latex = (
        f"\\left(t*\\cos({theta:.4f})-e^{{{M:.4f}\\left|t\\right|}}\\cdot"
        f"\\sin(0.3t)\\sin({theta:.4f})+{X:.4f},"
        f"42+t*\\sin({theta:.4f})+e^{{{M:.4f}\\left|t\\right|}}\\cdot"
        f"\\sin(0.3t)\\cos({theta:.4f})\\right)"
    )
    print("\nDesmos / LaTeX string:")
    print(latex)
    with open("results/desmos_equation.txt", "w") as f:
        f.write(latex + "\n")


if __name__ == "__main__":
    main()
