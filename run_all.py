"""One-command reproducibility entry point."""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def run(script):
    subprocess.run([sys.executable, str(ROOT / "src" / script)], check=True)


if __name__ == "__main__":
    run("fit_curve.py")
    run("verify_fit.py")
    run("plot_fit.py")
    run("plot_diagnostics.py")
    run("sensitivity.py")
