# Parametric Curve Fitting — R&D Assignment

Recover the unknown parameters `θ`, `M`, and `X` from an unordered set of `(x, y)` points.

## Problem

The curve is:

```text
x(t) = t*cos(θ) - exp(M*|t|)*sin(0.3*t)*sin(θ) + X

y(t) = 42 + t*sin(θ) + exp(M*|t|)*sin(0.3*t)*cos(θ)
```

with:

```text
6 < t < 60
```

The task is to recover:

```text
θ, M, X
```

from the supplied 1,500 unordered points.

## Final Answer

| Parameter | Recovered Value |      Allowed Range |
| --------- | --------------: | -----------------: |
| `θ`       |  **29.999973°** |     `0° < θ < 50°` |
| `M`       | **0.029999997** | `-0.05 < M < 0.05` |
| `X`       |   **54.999998** |      `0 < X < 100` |

Rounded to the required precision:

```text
θ = 30°
M = 0.03
X = 55
```

Equivalent value in radians:

```text
θ = 0.523599 rad
```

## Desmos

[Open the fitted curve in Desmos](https://www.desmos.com/calculator/jjcufejdax)

Copy-paste version:

```text
\left(t*\cos(0.523599)-e^{0.03\left|t\right|}\cdot\sin(0.3t)\sin(0.523599)+55,
42+t*\sin(0.523599)+e^{0.03\left|t\right|}\sin(0.3t)\cos(0.523599)\right)
```

with:

```text
6 < t < 60
```

# 1. Main Insight

The points are unordered, so we don't know the `t` value for each point.

One option would be to introduce a separate `t_i` for every observation. With 1,500 points this becomes a large optimization problem.

There is a simpler way because the equation has a useful geometric structure.

Let:

```text
A(t) = t

B(t) = exp(M*|t|) * sin(0.3*t)
```

Then the curve can be written as:

```text
[x - X]   [ cos(θ)  -sin(θ) ] [ A(t) ]
[y - 42] = [ sin(θ)   cos(θ) ] [ B(t) ]
```

So basically the original curve is:

```text
(t, exp(M*|t|)*sin(0.3*t))
```

rotated by `θ` and shifted by `(X, 42)`.

This means we can undo the rotation and recover `t` directly.

# 2. Analytical De-Rotation

For a candidate `(θ, X)`, rotate each point back:

```text
t_est =
(x - X)*cos(θ)
+
(y - 42)*sin(θ)
```

and:

```text
B_est =
-(x - X)*sin(θ)
+
(y - 42)*cos(θ)
```

For the correct parameters, we should have:

```text
B_est ≈ exp(M*|t_est|) * sin(0.3*t_est)
```

So instead of finding 1,500 unknown `t` values, we only need to fit:

```text
θ, M, X
```

This is the main trick used in the solution.

The inverse rotation also means we don't need to first order the points or do point-to-curve matching during the main optimization.

# 3. Optimization

For every candidate parameter set `(θ, M, X)`:

1. De-rotate all the observed points.
2. Get `t_est` and `B_est`.
3. Compute the model value:

```text
B_model = exp(M*|t_est|) * sin(0.3*t_est)
```

4. Use the residual:

```text
r = B_est - B_model
```

5. Add a soft penalty if `t_est` goes outside `6 < t < 60`.
6. Minimize the residual using bounded nonlinear least squares.

The implementation uses SciPy's Trust Region Reflective (`trf`) solver.

I used 60 different initial points with a fixed random seed (`42`) to check that the result is not dependent on a lucky initialization.

Parameter bounds used in the code:

```text
θ: 0.01° to 49.99°
M: -0.0499 to 0.0499
X: 0.01 to 99.99
```

The small margins are just to avoid starting exactly on the boundary.

# 4. Multi-Start Stability

A single nonlinear optimization run can sometimes end up in a different local solution depending on where it starts.

So the fit is run 60 times with different starting values.

All 60 runs converged to the same solution basin.

Best values:

```text
θ = 29.9999729°
M = 0.0299999969
X = 54.9999982
```

Final least-squares cost:

```text
9.11 × 10^-9
```

The values were consistent across the runs to the displayed precision.

# 5. Validation

The fitting objective and the final validation are kept separate.

During fitting, the points are transformed into `(t_est, B_est)` space because this gives a simple 3-variable optimization problem.

After fitting, the recovered parameters are put back into the original parametric equation and checked in `(x, y)` space.

Since the input points are unordered, I first recover their `t` values and sort them. Then I linearly interpolate the observations onto a common uniform `t` grid.

The fitted curve is evaluated on the same grid using 5,000 points.

The observed points cover approximately:

```text
6.0494 < t < 59.9952
```

which is inside the allowed range.

For each grid point, the L1 error is:

```text
|x_pred - x_ref| + |y_pred - y_ref|
```

The following values are reported:

```text
Mean L1
95th percentile L1
Maximum L1
```

A separate point-to-curve L1 check is also kept as an extra geometric check, but it is not the main validation metric.

# 6. Validation Results

Using the final rounded parameters:

```text
θ = 30°
M = 0.03
X = 55
```

the uniform-grid validation gives:

| Metric                  |           Result |
| ----------------------- | ---------------: |
| Mean uniform-grid L1    | **0.0001744963** |
| Maximum uniform-grid L1 | **0.0090396116** |
| 95th percentile L1      | **0.0005913611** |

These values are calculated by comparing the predicted and interpolated reference curve at the same uniform `t` values.

The full results are saved in:

```text
results/validation.txt
```

# 7. Why the Optimization Cost and L1 Error Are Different

The optimizer works in the de-rotated coordinate system:

```text
(x, y) → (t_est, B_est)
```

and minimizes:

```text
B_est - exp(M*|t_est|)*sin(0.3*t_est)
```

The final validation instead reconstructs the curve in `(x, y)` space and compares corresponding points on a uniform `t` grid.

So the optimizer cost and the final L1 error are measuring different things, and they don't have to be numerically the same.

In simple terms:

```text
optimization
    ↓
recover θ, M, X

validation
    ↓
check the reconstructed curve
```

# 8. Visual Results

## Fitted Curve

The supplied points and fitted curve are plotted together:

![Observed points and fitted parametric curve](results/fit_plot.png)

This gives a quick visual check that the fitted curve follows the data and oscillation pattern.

## L1 Residual

![L1 residual versus t](results/residual_vs_t.png)

This shows how the reconstruction error changes along `t`.

## Local Parameter Sensitivity

![Local parameter sensitivity](results/sensitivity.png)

This shows how the curve changes when the recovered parameters are perturbed around the fitted values.

# 9. Reproducibility

Clone the repository:

```bash
git clone https://github.com/Prianshu7296/Flam-RnD-Assignment.git
cd Flam-RnD-Assignment
```

Install dependencies:

```bash
pip install -r requirements.txt
```

Run the fitting:

```bash
python src/fit_curve.py
```

Run the validation:

```bash
python src/verify_fit.py
```

Regenerate the plots:

```bash
python src/plot_fit.py
```

The fitting script writes:

```text
results/fit_result.txt
```

and the validation script writes:

```text
results/validation.txt
```

# 10. Repository Structure

```text
Flam-RnD-Assignment/
├── README.md
├── requirements.txt
├── .gitignore
│
├── data/
│   └── xy_data.csv
│
├── src/
│   ├── fit_curve.py
│   ├── verify_fit.py
│   └── plot_fit.py
│
└── results/
    ├── fit_result.txt
    ├── validation.txt
    ├── desmos_equation.txt
    ├── fit_plot.png
    ├── residual_vs_t.png
    └── sensitivity.png
```

# 11. Design Decisions

### Analytical reduction

The inverse rotation lets us recover `t` without optimizing a separate `t` for every point.

### Bounded optimization

The search is kept inside the parameter ranges given in the assignment.

### Multi-start

60 starting points are used to check sensitivity to initialization.

### Fixed seed

The same initialization can be reproduced using seed `42`.

### Separate validation

The final curve is checked separately from the fitting objective.

### Visual checks

Plots are included to make it easier to inspect the fit and residuals.

# 12. Final Recovered Curve

Using:

```text
θ = 30°
M = 0.03
X = 55
```

the final curve is:

```text
x(t) =
t*cos(30°)
- exp(0.03*|t|)*sin(0.3*t)*sin(30°)
+ 55

y(t) =
42
+ t*sin(30°)
+ exp(0.03*|t|)*sin(0.3*t)*cos(30°)
```

for:

```text
6 < t < 60
```

# Conclusion

The recovered parameters are:

```text
θ = 30°
M = 0.03
X = 55
```

The main idea was to use the fact that the given curve is just a rotated and translated version of a simpler curve.

By undoing that rotation, the unordered point problem becomes a 3-parameter optimization problem instead of trying to estimate a `t` value separately for every observation.

The solution was tested from 60 different initializations, and the final curve was also checked independently using a uniform `t` grid.

Mean uniform-grid L1 error:

```text
1.744963 × 10^-4
```




