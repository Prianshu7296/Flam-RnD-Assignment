import numpy as np
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fit_curve import de_rotate
from verify_fit import curve_xy


def test_curve_and_derotation_are_inverse():
    t = np.linspace(6.2, 59.8, 100)
    theta, M, X = np.deg2rad(30.0), 0.03, 55.0
    x, y = curve_xy(theta, M, X, t)
    recovered_t, B, _ = de_rotate(theta, M, X, x, y)
    expected_B = np.exp(M * np.abs(t)) * np.sin(0.3 * t)

    assert np.max(np.abs(recovered_t - t)) < 1e-10
    assert np.max(np.abs(B - expected_B)) < 1e-10


def test_parameter_bounds():
    assert 0 < np.deg2rad(30) < np.deg2rad(50)
    assert -0.05 < 0.03 < 0.05
    assert 0 < 55 < 100
