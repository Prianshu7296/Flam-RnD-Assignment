"""
Recover theta, M, X for the assigned parametric curve.

Core idea:
    after undoing the rotation/translation,
        t_est = (x-X) cos(theta) + (y-42) sin(theta)
        B_est = -(x-X) sin(theta) + (y-42) cos(theta)

    and the model requires
        B_est = exp(M |t_est|) sin(0.3 t_est)

This keeps the original Prianshu approach: a 3-parameter nonlinear
least-squares problem. The additions are only better initialization and
an independent global-search cross-check.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "data" / "xy_data.csv"
RESULTS_DIR = ROOT / "results"

T_LO, T_HI = 6.0, 60.0
THETA_BOUNDS = (np.deg2rad(0.001), np.deg2rad(49.999))
M_BOUNDS = (-0.049999, 0.049999)
X_BOUNDS = (0.001, 99.999)


def load_data(path=DATA_PATH):
    df = pd.read_csv(path)
    return df["x"].to_numpy(float), df["y"].to_numpy(float)


def de_rotate(theta, M, X, x, y):
    c, s = np.cos(theta), np.sin(theta)
    t_est = (x - X) * c + (y - 42.0) * s
    B_est = -(x - X) * s + (y - 42.0) * c
    B_model = np.exp(M * np.abs(t_est)) * np.sin(0.3 * t_est)
    return t_est, B_est, B_model


def residuals(params, x, y):
    theta, M, X = params
    t_est, B_est, B_model = de_rotate(theta, M, X, x, y)

    r_wiggle = B_est - B_model

    # Keep the optimizer inside the assignment's t-domain.
    pen_lo = np.clip(T_LO - t_est, 0.0, None)
    pen_hi = np.clip(t_est - T_HI, 0.0, None)

    return np.concatenate([r_wiggle, 5.0 * pen_lo, 5.0 * pen_hi])


def estimate_x(theta, x, y):
    """
    Since the projection along the backbone is
        q = t + X cos(theta),
    and t is restricted to (6, 60), use the midpoint of the observed
    projection range to estimate X.
    """
    c = np.cos(theta)
    q = x * c + (y - 42.0) * np.sin(theta)
    x0 = (q.min() + q.max() - (T_LO + T_HI)) / (2.0 * c)
    return float(np.clip(x0, *X_BOUNDS))


def estimate_m(theta, X, x, y):
    """
    From B = exp(M t) sin(0.3t),
        log(|B/sin(0.3t)|) = M t.
    Ignore points close to sine zeros.
    """
    t, B, _ = de_rotate(theta, 0.0, X, x, y)
    sine = np.sin(0.3 * t)
    mask = (
        (t > T_LO)
        & (t < T_HI)
        & (np.abs(sine) > 0.20)
        & (np.abs(B) > 1e-12)
    )
    if mask.sum() < 20:
        return 0.0

    slope = np.polyfit(
        t[mask],
        np.log(np.abs(B[mask] / sine[mask])),
        1,
    )[0]
    return float(np.clip(slope, *M_BOUNDS))


def informed_starts(x, y):
    """
    Coarse theta sweep gives a small set of geometry-informed starts.
    This is deterministic and supplements, rather than replaces,
    Prianshu's multi-start least-squares strategy.
    """
    candidates = []
    for theta_deg in np.linspace(0.5, 49.5, 99):
        theta = np.deg2rad(theta_deg)
        X = estimate_x(theta, x, y)
        M = estimate_m(theta, X, x, y)
        score = np.sum(residuals([theta, M, X], x, y) ** 2)
        candidates.append((score, theta, M, X))

    candidates.sort(key=lambda z: z[0])
    return [np.array([theta, M, X]) for _, theta, M, X in candidates[:5]]


def global_search(x, y, n_starts=60, seed=42):
    rng = np.random.default_rng(seed)

    starts = informed_starts(x, y)

    # Preserve the original multi-start idea: 60 total LS runs.
    while len(starts) < n_starts:
        starts.append(
            np.array([
                rng.uniform(*THETA_BOUNDS),
                rng.uniform(*M_BOUNDS),
                rng.uniform(*X_BOUNDS),
            ])
        )

    lower = [THETA_BOUNDS[0], M_BOUNDS[0], X_BOUNDS[0]]
    upper = [THETA_BOUNDS[1], M_BOUNDS[1], X_BOUNDS[1]]

    best = None
    costs = []

    for p0 in starts[:n_starts]:
        try:
            sol = least_squares(
                residuals,
                p0,
                args=(x, y),
                bounds=(lower, upper),
                method="trf",
                max_nfev=5000,
            )
        except Exception:
            continue

        costs.append(sol.cost)
        if best is None or sol.cost < best.cost:
            best = sol

    if best is None:
        raise RuntimeError("All least-squares starts failed.")

    return best, np.array(costs)


def scalar_objective(params, x, y):
    r = residuals(params, x, y)
    return float(np.dot(r, r))


def differential_evolution_check(x, y, seed=42):
    """
    Independent global-search cross-check. It is not used to choose the
    final answer; it is used to test whether the LS solution is reproducible
    with a different optimization strategy.
    """
    result = differential_evolution(
        lambda p: scalar_objective(p, x, y),
        bounds=[THETA_BOUNDS, M_BOUNDS, X_BOUNDS],
        seed=seed,
        popsize=8,
        maxiter=120,
        tol=1e-9,
        polish=False,
        updating="immediate",
    )
    return result


def save_result(best, ls_costs, de_result):
    RESULTS_DIR.mkdir(exist_ok=True)

    theta, M, X = best.x
    theta_deg = np.rad2deg(theta)

    with open(RESULTS_DIR / "fit_result.txt", "w") as f:
        f.write(f"theta_rad={theta:.10f}\n")
        f.write(f"theta_deg={theta_deg:.8f}\n")
        f.write(f"M={M:.10f}\n")
        f.write(f"X={X:.10f}\n")
        f.write(f"least_squares_cost={best.cost:.12g}\n")
        f.write(f"ls_runs={len(ls_costs)}\n")
        f.write(f"ls_runs_with_cost_below_1e-4={int(np.sum(ls_costs < 1e-4))}\n")
        f.write(f"de_theta_deg={np.rad2deg(de_result.x[0]):.8f}\n")
        f.write(f"de_M={de_result.x[1]:.10f}\n")
        f.write(f"de_X={de_result.x[2]:.10f}\n")
        f.write(f"de_objective={de_result.fun:.12g}\n")


def main():
    x, y = load_data()
    print(f"Loaded {len(x)} points")

    best, ls_costs = global_search(x, y, n_starts=60, seed=42)
    de_result = differential_evolution_check(x, y, seed=42)

    theta, M, X = best.x
    print("\n=== Final fit ===")
    print(f"theta = {np.rad2deg(theta):.8f} deg")
    print(f"M     = {M:.10f}")
    print(f"X     = {X:.10f}")
    print(f"LS cost = {best.cost:.12g}")

    print("\n=== Independent DE cross-check ===")
    print(f"theta = {np.rad2deg(de_result.x[0]):.8f} deg")
    print(f"M     = {de_result.x[1]:.10f}")
    print(f"X     = {de_result.x[2]:.10f}")
    print(f"objective = {de_result.fun:.12g}")

    save_result(best, ls_costs, de_result)


if __name__ == "__main__":
    main()
