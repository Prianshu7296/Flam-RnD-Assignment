# Parametric Curve Fitting — R&D Assignment

Recover the unknown parameters `theta`, `M`, and `X` from an unordered set of `(x, y)` points.

The curve is defined by the parametric equations:

```text
x(t) = t*cos(theta) - exp(M*|t|)*sin(0.3*t)*sin(theta) + X

y(t) = 42 + t*sin(theta) + exp(M*|t|)*sin(0.3*t)*cos(theta)
```

with domain:

```text
6 < t < 60
```

## Search Bounds

| Parameter | Range              |
| --------- | ------------------ |
| `theta`   | `0° < theta < 50°` |
| `M`       | `-0.05 < M < 0.05` |
| `X`       | `0 < X < 100`      |

---

# Final Answer

| Parameter | Recovered Value |
| --------- | --------------: |
| `theta`   |    **30.0000°** |
| `M`       |      **0.0300** |
| `X`       |     **55.0000** |

Equivalent angle in radians:

```text
theta = 0.523599 rad
```

The recovered parameters are well inside the allowed search bounds.

## Desmos

[Open fitted curve in Desmos](https://www.desmos.com/calculator/jjcufejdax)

Copy-paste ready equation:

```text
\left(t*\cos(0.523599)-e^{0.03\left|t\right|}\cdot\sin(0.3t)\sin(0.523599)+55,42+t*\sin(0.523599)+e^{0.03\left|t\right|}\cdot\sin(0.3t)\cos(0.523599)\right)
```

---

# 1. Key Mathematical Insight

The main challenge is that the supplied `(x, y)` points are unordered and their corresponding values of `t` are unknown.

A naive approach would introduce one unknown `t_i` for every observed point. With 1,500 points, this would create a very high-dimensional optimization problem.

Instead, the structure of the curve can be exploited analytically.

Define:

```text
A(t) = t

B(t) = exp(M*|t|) * sin(0.3*t)
```

The curve can then be written as a rotation followed by a translation:

```text
[x - X]   [ cos(theta)  -sin(theta) ] [ A(t) ]
[y - 42] =[ sin(theta)   cos(theta) ] [ B(t) ]
```

Because the rotation matrix is orthogonal, its inverse is simply its transpose.

Therefore, for any candidate `theta` and `X`:

```text
t_est =
(x - X)*cos(theta)
+
(y - 42)*sin(theta)
```

and

```text
B_est =
-(x - X)*sin(theta)
+
(y - 42)*cos(theta)
```

For the correct parameters:

```text
B_est ≈ exp(M*|t_est|) * sin(0.3*t_est)
```

This is the key simplification.

Instead of optimizing thousands of unknown `t_i` values, the problem is reduced to only three unknown parameters:

```text
theta, M, X
```

This makes the fitting problem substantially lower-dimensional and easier to optimize robustly.

---

# 2. Why the Inverse Rotation Works

The original curve is a simple base curve:

```text
(t, exp(M*|t|) * sin(0.3*t))
```

which is then:

1. Rotated by `theta`
2. Translated by `(X, 42)`

Undoing the translation and rotation recovers the coordinates in the original curve frame.

In that frame, the first coordinate is directly equal to `t`.

Therefore:

```text
t_est = transformed x-coordinate
```

and the second transformed coordinate can be compared directly with the exponential-sine model.

This avoids explicit point-to-curve correspondence search during optimization.

---

# 3. Optimization Strategy

The optimization is performed in only three variables:

```text
(theta, M, X)
```

For each candidate parameter set, the implementation:

1. Applies the inverse rotation to every observed point.
2. Recovers `t_est` and `B_est`.
3. Evaluates the model

```text
B_model = exp(M*|t_est|) * sin(0.3*t_est)
```

4. Computes the transformed residual

```text
residual = B_est - B_model
```

5. Applies soft penalties when recovered `t` values move outside the allowed interval.
6. Optimizes the residual using bounded nonlinear least squares.

The implementation uses 60 independent random initializations with a fixed random seed.

Each initialization uses SciPy's Trust Region Reflective (`trf`) least-squares solver.

The purpose of the multi-start procedure is to reduce sensitivity to the initial guess and check whether different starting points consistently converge to the same parameter region.

All successful runs converge to the same solution basin:

```text
theta ≈ 30°
M     ≈ 0.03
X     ≈ 55
```

---

# 4. Optimization Results

The final recovered parameters are:

```text
theta = 30.0000°
M     = 0.0300
X     = 55.0000
```

The repeated convergence across many initializations provides an additional consistency check on the solution.

---

# 5. Validation

The optimization residual and the final validation metric are treated separately.

The optimizer works in the transformed coordinate system because the inverse rotation provides a simple residual for estimating the parameters.

The final validation is performed directly against the reconstructed parametric curve.

## Validation procedure

The fitted curve is sampled densely over the allowed domain:

```text
6 <= t <= 60
```

using 5,000 uniformly spaced values of `t`.

For every observed point `(x_i, y_i)`, the implementation computes the minimum L1 distance to the sampled fitted curve:

```text
d_i =
min(
    |x(t) - x_i| + |y(t) - y_i|
)
```

over all sampled curve points.

The final reported metrics are computed from these point-to-curve distances.

This provides a dense numerical approximation of the reconstruction error using the L1 distance specified by the assignment.

---

# 6. Validation Results

Results on the supplied 1,500-point dataset:

| Metric              |            Value |
| ------------------- | ---------------: |
| Mean L1 distance    | **0.0001744963** |
| 95th percentile L1  | **0.0005913611** |
| Maximum L1 distance | **0.0090396116** |

The mean error is approximately:

```text
1.745 × 10^-4
```

The small mean error indicates that the recovered parameters reproduce the supplied curve closely.

The maximum error is larger than the mean because the validation compares a finite set of observed points against a discretely sampled representation of the fitted curve.

---

# 7. Visual Validation

## Fitted Curve

The supplied data points are shown together with the reconstructed parametric curve.

![Observed points and fitted parametric curve](results/fit_plot.png)

The fitted curve closely follows the observed oscillatory trajectory across the full domain.

---

## L1 Residual versus t

The residual plot shows how reconstruction error varies across the curve.

![L1 residual versus t](results/residual_vs_t.png)

This helps identify regions where the fitted model differs most from the observed data.

---

## Local Parameter Sensitivity

The sensitivity plot shows how the fitted solution changes under local parameter perturbations.

![Local parameter sensitivity](results/sensitivity.png)

This provides an additional check on the stability of the recovered solution.

---

# 8. Reproducibility

Clone the repository:

```bash
git clone https://github.com/Prianshu7296/Flam-RnD-Assignment.git
cd Flam-RnD-Assignment
```

Install dependencies:

```bash
pip install -r requirements.txt
```

Run the fitting procedure:

```bash
python src/fit_curve.py
```

Run validation:

```bash
python src/verify_fit.py
```

Regenerate the figures:

```bash
python src/plot_fit.py
```

---

# 9. Repository Structure

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

# 10. Reproducibility and Design Choices

A few design choices are intentional:

### Analytical parameter reduction

The inverse-rotation transformation avoids introducing one separate `t_i` variable for every data point.

### Multi-start optimisation

Multiple initializations reduce dependence on a single starting point and make the optimization more robust to local minima.

### Fixed random seed

The initialization procedure is reproducible.

### Separate validation

The optimization residual is not used as the only measure of final fit quality. The fitted curve is evaluated separately using a dense point-to-curve L1 calculation.

---

# 11. Final Recovered Curve

The final recovered curve is:

```text
x(t) = t*cos(30°)
       - exp(0.03*|t|)*sin(0.3*t)*sin(30°)
       + 55

y(t) = 42
       + t*sin(30°)
       + exp(0.03*|t|)*sin(0.3*t)*cos(30°)
```

for:

```text
6 < t < 60
```

---

# Conclusion

The recovered parameters are:

```text
theta = 30°
M     = 0.03
X     = 55
```

The main contribution of the approach is the analytical inverse-rotation step.

It transforms the unordered point-cloud problem into a three-parameter nonlinear optimization problem without explicitly optimizing the unknown `t` value for every data point.

The final fitted curve closely reproduces the supplied 1,500-point dataset, with a mean L1 reconstruction error of approximately:

```text
1.745 × 10^-4
```

This combination of analytical reduction, multi-start optimization, and independent curve-level validation provides a compact and reproducible solution to the assignment.



