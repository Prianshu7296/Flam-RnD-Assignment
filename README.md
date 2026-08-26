# Parametric Curve Fitting — R&D Assignment

Recover the unknown parameters \(\theta\), \(M\), and \(X\) from an unordered set of \((x,y)\) points.

The curve is defined by

$$
\begin{aligned}
x(t) &= t\cos\theta
      - e^{M|t|}\sin(0.3t)\sin\theta + X \\
y(t) &= 42 + t\sin\theta
      + e^{M|t|}\sin(0.3t)\cos\theta
\end{aligned}
$$

with

$$
6 < t < 60.
$$

## Search Bounds

| Parameter  | Range                           |
| ---------- | ------------------------------- |
| \(\theta\) | \(0^\circ < \theta < 50^\circ\) |
| \(M\)      | \(-0.05 < M < 0.05\)            |
| \(X\)      | \(0 < X < 100\)                 |

---

# Final Answer

| Parameter  | Recovered Value |
| ---------- | --------------: |
| \(\theta\) |    **30.0000°** |
| \(M\)      |      **0.0300** |
| \(X\)      |     **55.0000** |

Equivalent angle in radians:

$$
\theta \approx 0.523599\text{ rad}
$$

The recovered parameters lie comfortably inside the allowed search bounds.

## Desmos

[Open fitted curve in Desmos](https://www.desmos.com/calculator/jjcufejdax)

Copy-paste ready equation:

```text
\left(t*\cos(0.523599)-e^{0.03\left|t\right|}\cdot\sin(0.3t)\sin(0.523599)+55,42+t*\sin(0.523599)+e^{0.03\left|t\right|}\cdot\sin(0.3t)\cos(0.523599)\right)
```

---

# 1. Key Mathematical Insight

The main difficulty is that the supplied \((x,y)\) points are unordered and the corresponding values of \(t\) are unknown.

A direct approach would introduce one additional unknown \(t_i\) for every observed point. That would turn the problem into a very high-dimensional nonlinear optimization problem.

Instead, the structure of the curve can be exploited analytically.

Define

$$
A(t)=t
$$

and

$$
B(t)=e^{M|t|}\sin(0.3t).
$$

Then the curve becomes

$$
\begin{pmatrix}
x-X\\
y-42
\end{pmatrix}
=
\begin{pmatrix}
\cos\theta & -\sin\theta\\
\sin\theta & \cos\theta
\end{pmatrix}
\begin{pmatrix}
A(t)\\
B(t)
\end{pmatrix}.
$$

Therefore, the observed curve is simply the base curve

$$
\left(t,\ e^{M|t|}\sin(0.3t)\right)
$$

after a rotation by \(\theta\) and a translation by \((X,42)\).

The important consequence is that the rotation can be inverted analytically.

---

# 2. Inverse Rotation

Because a rotation matrix is orthogonal, its inverse is its transpose.

For any candidate \((\theta,X)\),

$$
t_{\text{est}}
=
(x-X)\cos\theta+(y-42)\sin\theta
$$

and

$$
B_{\text{est}}
=
-(x-X)\sin\theta+(y-42)\cos\theta.
$$

For the correct parameters, the transformed point should satisfy

$$
B_{\text{est}}
\approx
e^{M|t_{\text{est}}|}
\sin(0.3t_{\text{est}}).
$$

This removes the need to explicitly search for point-to-curve correspondences.

The original unordered 2-D fitting problem is therefore reduced to a bounded optimization problem involving only three unknowns:

$$
\boxed{(\theta,M,X)}
$$

This reduction is the main mathematical idea used in the solution.

---

# 3. Optimization Strategy

For each candidate parameter vector

$$
(\theta,M,X)
$$

the implementation performs the following steps:

1. Apply the inverse rotation to every observed point.
2. Recover an estimated \(t\) value and transformed oscillation value \(B\).
3. Evaluate the model

$$
B_{\text{model}}
=
e^{M|t_{\text{est}}|}
\sin(0.3t_{\text{est}}).
$$

4. Compute the transformed-coordinate residual

$$
r =
B_{\text{est}}-B_{\text{model}}.
$$

5. Penalize recovered \(t\) values outside the allowed domain \(6<t<60\).
6. Optimize the resulting residual vector using bounded nonlinear least squares.

The implementation uses **60 independent random initializations** with a fixed random seed. Each initialization is optimized using SciPy's Trust Region Reflective (`trf`) least-squares solver.

Using multiple starting points reduces dependence on the initial guess and provides a practical check that the optimization is repeatedly reaching the same basin.

The best solution converges to:

```text
theta = 30.0000 degrees
M     = 0.0300
X     = 55.0000
```

---

# 4. Why the Optimization Works

The sinusoidal component is nonlinear in \(t\), but the inverse rotation gives a direct estimate of \(t\) for every candidate \((\theta,X)\).

That means the implementation does **not** optimize thousands of independent \(t_i\) values.

Instead, every observed point is transformed analytically, and the only variables being optimized are:

```text
theta
M
X
```

This reduces the dimensionality dramatically and makes the problem much easier to solve robustly.

The fitted solution is also checked across many different initial parameter guesses rather than relying on a single optimizer run.

---

# 5. Validation

The fitting objective and the final validation are intentionally separated.

The optimizer works in the transformed coordinate system because this provides a simple differentiable residual for recovering the parameters.

The final validation instead measures the geometric reconstruction error between the observed data and the fitted parametric curve.

## Validation procedure

The fitted curve is uniformly sampled over

$$
6 \leq t \leq 60
$$

using 5,000 points.

For every observed point \((x_i,y_i)\), the implementation computes its L1 distance to the closest sampled point on the fitted curve:

$$
d_i
=
\min_t
\left(
|x(t)-x_i|
+
|y(t)-y_i|
\right).
$$

The reported metrics are then computed over all supplied data points.

This gives a dense numerical approximation to the point-to-curve L1 reconstruction error.

The validation is implemented independently in:

```text
src/verify_fit.py
```

and therefore does not depend on the optimizer's internal residual.

---

# 6. Validation Results

Results on the supplied **1,500-point dataset**:

| Metric              |            Value |
| ------------------- | ---------------: |
| Mean L1 distance    | **0.0001744963** |
| 95th percentile L1  | **0.0005913611** |
| Maximum L1 distance | **0.0090396116** |

The mean reconstruction error is approximately

$$
1.745\times10^{-4}.
$$

The fitted curve therefore stays very close to the supplied point cloud across the full parameter range.

The maximum error is larger than the mean, which is expected when evaluating a dense oscillatory curve against a finite set of sampled points.

---

# 7. Visual Validation

## Fitted Curve

The supplied data points are shown together with the reconstructed parametric curve.

![Observed points and fitted parametric curve](results/fit_plot.png)

The fitted curve follows the observed oscillatory trajectory across the full domain.

---

## L1 Residual

The residual plot shows how reconstruction error varies across the curve.

![L1 residual versus t](results/residual_vs_t.png)

This provides a more useful diagnostic than the fitted-curve plot alone because it shows where the model deviates most strongly from the supplied data.

---

## Parameter Sensitivity

The sensitivity plot shows the local effect of perturbing the recovered parameters.

![Parameter sensitivity](results/sensitivity.png)

This provides an additional check that the reported solution is not simply an arbitrary point in parameter space.

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

This performs the multi-start nonlinear least-squares search and writes the recovered parameters to:

```text
results/fit_result.txt
```

Run the independent validation:

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

# 10. Final Result

The recovered parameters are:

```text
theta = 30°
M     = 0.03
X     = 55
```

Therefore the recovered curve is

$$
\boxed{
\begin{aligned}
x(t) &= t\cos(30^\circ)
- e^{0.03|t|}\sin(0.3t)\sin(30^\circ)
+55\\
y(t) &= 42+t\sin(30^\circ)
+e^{0.03|t|}\sin(0.3t)\cos(30^\circ)
\end{aligned}}
$$

for

$$
6<t<60.
$$

---

# Conclusion

The solution uses the geometric structure of the parametric equation rather than treating the problem as a generic black-box curve fit.

The key step is the analytical inverse rotation, which removes the unknown point-to-\(t\) correspondence and reduces the problem to a three-parameter nonlinear optimization.

The final parameters

$$
\boxed{\theta=30^\circ,\quad M=0.03,\quad X=55}
$$

produce a close reconstruction of the supplied 1,500-point dataset, with a mean L1 point-to-curve error of approximately

$$
\boxed{1.745\times10^{-4}}.
$$


