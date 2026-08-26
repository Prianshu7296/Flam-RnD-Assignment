"""Generate the main data-vs-model overlay."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from verify_fit import load_result, curve_xy

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "data" / "xy_data.csv"
OUT_PATH = ROOT / "results" / "fit_plot.png"


def main():
    res = load_result()
    df = pd.read_csv(DATA_PATH)
    x, y = df["x"].to_numpy(float), df["y"].to_numpy(float)

    t = np.linspace(6, 60, 5000)
    xf, yf = curve_xy(res["theta_rad"], res["M"], res["X"], t)

    plt.figure(figsize=(8, 6))
    plt.scatter(x, y, s=7, alpha=0.35, label="data")
    plt.plot(xf, yf, linewidth=2, label="fitted curve")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(
        f"Curve fit: theta={res['theta_deg']:.3f}°, "
        f"M={res['M']:.5f}, X={res['X']:.3f}"
    )
    plt.legend()
    plt.gca().set_aspect("equal", adjustable="datalim")
    plt.tight_layout()
    plt.savefig(OUT_PATH, dpi=180)
    plt.close()
    print(f"Saved {OUT_PATH}")


if __name__ == "__main__":
    main()
