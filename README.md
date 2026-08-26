# Parametric Curve Fitting — R&D / AI Assignment

## Final result

The recovered parameters are:

| Parameter | Result | Given range |
|---|---:|---:|
| θ | **30.0000°** | 0° < θ < 50° |
| M | **0.030000** | −0.05 < M < 0.05 |
| X | **55.0000** | 0 < X < 100 |

## 1. Mathematical reduction

The curve can be written as

```text
A(t) = t
B(t) = exp(M|t|) sin(0.3t)

[x-X]   [ cosθ  -sinθ ] [A]
[y-42] = [ sinθ   cosθ ] [B]
```

So for any candidate `(θ, X)` the rotation can be inverted analytically:

```text
t_est = (x-X)cosθ + (y-42)sinθ
B_est = -(x-X)sinθ + (y-42)cosθ
```

At the correct parameters,

```text
B_est = exp(M|t_est|) sin(0.3t_est)
```

This removes the unknown point-to-t correspondence and leaves only three parameters to optimize.

## 2. Optimization

The original multi-start nonlinear least-squares approach is retained.

Three additions make the search more defensible:

1. A coarse θ sweep gives geometry-informed initial values.
2. X is initialized from the known t-domain `(6,60)`.
3. M is initialized from
   `log(|B/sin(0.3t)|) ≈ M t`, away from sine zeros.

The final solution is still selected by bounded `scipy.optimize.least_squares` over 60 starts.

As an independent check, the same problem is solved with Differential Evolution. The two methods converge to the same parameters to numerical precision.

## 3. Validation

The supplied points do not contain their original `t` values. Therefore the validation first recovers `t_est` using the fitted geometry, sorts the observations by `t_est`, interpolates the observed curve onto a uniform t-grid, and compares the analytical prediction and observed interpolation at the same t values.

The reported metric is:

```text
mean(|x_pred-x_obs| + |y_pred-y_obs|)
```

on 5,000 uniformly spaced points over the observed t-domain.

This validation is kept separate from the fitting objective.

## 4. Diagnostics

`results/fit_plot.png`

Raw data and fitted curve overlay.

`results/derotation_diagnostic.png`

Shows that after undoing the rotation and translation, the recovered one-dimensional signal follows `exp(Mt)sin(0.3t)`.

`results/residual_vs_t.png`

Shows the L1 error at each point of the uniform validation grid.

`results/sensitivity.png`

Shows how the validation error changes under small perturbations of θ, M and X.

## 5. Reproduce

```bash
pip install -r requirements.txt
python run_all.py
```

Or run individual steps:

```bash
python src/fit_curve.py
python src/verify_fit.py
python src/plot_fit.py
python src/plot_diagnostics.py
python src/sensitivity.py
pytest -q
```

## 6. Why the solution is robust

The fit is not accepted merely because the final parameters look plausible.

There are three independent checks:

- 60 bounded nonlinear least-squares starts.
- A geometry-informed initialization based on the curve structure.
- An independent Differential Evolution optimization.

The final result is then evaluated with a separate uniform-grid L1 procedure and visual/geometric diagnostics.

## Repository structure

```text
├── README.md
├── requirements.txt
├── run_all.py
├── data/
│   └── xy_data.csv
├── src/
│   ├── fit_curve.py
│   ├── verify_fit.py
│   ├── plot_fit.py
│   ├── plot_diagnostics.py
│   └── sensitivity.py
├── tests/
│   └── test_geometry.py
└── results/
    ├── fit_result.txt
    ├── validation.txt
    ├── sensitivity.txt
    ├── fit_plot.png
    ├── derotation_diagnostic.png
    ├── residual_vs_t.png
    └── sensitivity.png
```

## Notes

The assignment only asks for θ, M and X. The recovered `t` values are an internal fitting/validation device and are not part of the submitted answer.
