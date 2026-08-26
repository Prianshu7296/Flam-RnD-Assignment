# Parametric Curve Fitting — R&D Assignment

Recover the unknown parameters `θ`, `M`, and `X` from an unordered set of `(x, y)` points.

The curve is defined by:

```text
x(t) = t*cos(θ) - exp(M*|t|)*sin(0.3*t)*sin(θ) + X

y(t) = 42 + t*sin(θ) + exp(M*|t|)*sin(0.3*t)*cos(θ)
```

with:

```text
6 < t < 60
```

The objective is to recover:

```text
θ, M, X
```

from the supplied 1,500 unordered points.

---

## Final Answer

| Parameter | Recovered Value |      Allowed Range |
| --------- | --------------: | -----------------: |
| `θ`       |  **29.999973°** |     `0° < θ < 50°` |
| `M`       | **0.029999997** | `-0.05 < M < 0.05` |
| `X`       |   **54.999998** |      `0 < X < 100` |

Rounded to the precision required by the assignment:

```text
θ = 30°
M = 0.03
X = 55
```

Equivalent angle in radians:

```text
θ = 0.523599 rad
```

---

## Desmos

[Open the fitted curve in Desmos](https://www.desmos.com/calculator/jjcufejdax)

Copy-paste ready:

```text
\left(t*\cos(0.523599)-e^{0.03\left|t\right|}\cdot\sin(0.3t)\sin(0.523599)+55,42+t*\sin(0.523599)+e^{0.03\left|t\right|}\cdot\sin(0.3t)\cos(0.523599)\right)
```

---

# 1. Main Insight

The supplied points are unordered, so their corresponding parameter values `t` are unknown.

A direct formulation could introduce a separate unknown `t_i` for every observed point. With 1,500 observations, that would create a large nonlinear optimization problem with thousands of variables.

Instead, the structure of the equation can be exploited analytically.

Define:

```text
A(t) = t

B(t) = exp(M*|t|) * sin(0.3*t)
```

Then the curve becomes:

```text
[x - X]   [ cos(θ)  -sin(θ) ] [ A(t) ]
[y - 42] = [ sin(θ)   cos(θ) ] [ B(t) ]
```

So the observed curve is the base curve

```text
(t, exp(M*|t|)*sin(0.3*t))
```

after a rotation by `θ` and a translation by `(X, 42)`.

This geometric structure lets us undo the rotation analytically.

---

# 2. Analytical De-Rotation

For any candidate `(θ, X)`, transform each observed point back into the original curve coordinate system:

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

For the correct parameters, the transformed point must satisfy:

```text
B_est ≈ exp(M*|t_est|) * sin(0.3*t_est)
```

Therefore, instead of estimating an unknown `t_i` for every point, the entire dataset can be explained using only three unknowns:

```text
θ, M, X
```

This is the key reduction used in the solution.

### Why this matters

The transformation is a rigid rotation followed by a translation. The inverse rotation recovers the coordinate system in which the first coordinate is directly `t`.

That removes the need for:

```text
- point ordering
- explicit point-to-curve correspondence
- thousands of independent t variables
- nearest-neighbour matching during the main optimization
```

The problem is reduced to a low-dimensional nonlinear fit.

---

# 3. Optimization

The fitting objective is constructed in the de-rotated coordinate system.

For every candidate parameter vector:

```text
(θ, M, X)
```

the algorithm:

1. Inverts the rotation for all observed points.
2. Recovers `t_est` and `B_est`.
3. Evaluates the model:

```text
B_model = exp(M*|t_est|) * sin(0.3*t_est)
```

4. Computes the transformed residual:

```text
r = B_est - B_model
```

5. Adds soft penalties when `t_est` falls outside the permitted interval.
6. Minimizes the combined residual using bounded nonlinear least squares.

The implementation uses SciPy's Trust Region Reflective (`trf`) solver with:

```text
60 independent initializations
fixed random seed = 42
```

The parameter bounds used by the implementation are:

```text
θ: 0.01° to 49.99°
M: -0.0499 to 0.0499
X: 0.01 to 99.99
```

The small interior margins avoid starting exactly on an optimization boundary.

---

# 4. Multi-Start Stability

A single nonlinear optimization run can depend on its initial point.

To test this, 60 independent starting points are sampled uniformly inside the allowed parameter ranges.

All 60 runs converged successfully to the same solution basin.

Best recovered values:

```text
θ = 29.9999729°
M = 0.0299999969
X = 54.9999982
```

The final least-squares cost was approximately:

```text
9.11 × 10^-9
```

Across all 60 runs, the optimization converged to the same numerical solution to the displayed precision.

This is useful as a practical robustness check: the result is not dependent on a single lucky initialization.

---

# 5. Validation

The fitting objective and final curve validation are treated separately.

The optimizer minimizes the transformed-coordinate model residual because that gives a compact three-variable optimization problem.

After fitting, the recovered parameters are inserted back into the original parametric equations and the reconstructed curve is evaluated independently.

The validation script samples the fitted curve uniformly over:

```text
6 <= t <= 60
```

using 5,000 points.

For every observed point `(x_i, y_i)`, it computes the minimum L1 distance to the sampled fitted curve:

```text
d_i =
min(
    |x(t) - x_i| + |y(t) - y_i|
)
```

The final statistics are calculated over all 1,500 supplied points.

This provides a dense numerical approximation of the point-to-curve L1 reconstruction error.

---

# 6. Validation Results

Measured directly from the supplied 1,500-point dataset using the repository's validation procedure:

| Metric              |       Result |
| ------------------- | -----------: |
| Mean L1 distance    | **0.004074** |
| 95th percentile L1  | **0.008108** |
| Maximum L1 distance | **0.012357** |

The mean L1 error is approximately:

```text
4.07 × 10^-3
```

The validation is performed independently of the optimizer's residual calculation by reconstructing the curve from the saved fitted parameters.

---

# 7. Why the Validation Error Is Larger Than the Optimizer Residual

The optimization residual and the final L1 validation measure different quantities.

The optimizer operates after transforming every point into the curve's intrinsic coordinate system:

```text
(x, y) → (t_est, B_est)
```

and asks whether:

```text
B_est ≈ exp(M*|t_est|)*sin(0.3*t_est)
```

The final validation instead works in the original `(x, y)` space and asks how close each observed point is to the discretely sampled reconstructed curve under L1 distance.

Therefore, the two numbers are not expected to be identical.

This distinction is intentional:

```text
optimization objective
        ↓
parameter recovery

independent geometric validation
        ↓
curve reconstruction quality
```

---

# 8. Visual Results

## Fitted Curve

The supplied data points and reconstructed parametric curve are shown together.

![Observed points and fitted parametric curve](results/fit_plot.png)

The visual overlay provides a qualitative check that the recovered parameters reproduce the overall geometry and oscillation pattern.

---

## L1 Residual

The residual plot shows how the reconstruction error varies across the curve.

![L1 residual versus t](results/residual_vs_t.png)

This is useful for identifying whether the model has localized regions of larger error rather than relying only on a single aggregate metric.

---

## Local Parameter Sensitivity

The sensitivity plot shows the effect of perturbing the recovered parameters around the fitted solution.

![Local parameter sensitivity](results/sensitivity.png)

This provides an additional diagnostic of the local behaviour of the fitted solution.

---

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

Recover the parameters:

```bash
python src/fit_curve.py
```

Run independent curve validation:

```bash
python src/verify_fit.py
```

Regenerate the figures:

```bash
python src/plot_fit.py
```

The fitting script writes the recovered parameters to:

```text
results/fit_result.txt
```

and the validation script reads those saved values for its reconstruction check.

---

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

---

# 11. Design Decisions

### Analytical reduction instead of brute force

The inverse rotation eliminates the need to optimize one unknown `t` value per data point.

### Bounded optimization

The search is constrained to the parameter ranges specified by the assignment.

### Multi-start fitting

60 independent starting points are used to reduce sensitivity to initialization.

### Fixed random seed

The initialization process is reproducible.

### Separate validation

The final reconstructed curve is evaluated separately from the internal optimizer residual.

### Visual diagnostics

The repository includes:

```text
fit_plot.png
residual_vs_t.png
sensitivity.png
```

so that the numerical result can also be inspected visually.

---

# 12. Final Recovered Curve

Using the rounded recovered parameters:

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

---

# Conclusion

The final recovered parameters are:

```text
θ = 30°
M = 0.03
X = 55
```

The main idea is to exploit the fact that the given parametric curve is a rotated and translated version of a simpler curve.

By analytically undoing that rigid transformation, the unordered point-cloud problem is reduced to a bounded optimization over only three parameters:

```text
θ, M, X
```

The solution is then tested from 60 independent initializations, with all runs converging to the same parameter basin.

The fitted parameters reproduce the supplied curve closely, while the separate geometric validation provides an independent measure of the final reconstruction quality.




