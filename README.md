# Parametric Curve Fitting — R&D Assignment

Recover the unknown parameters \(\theta\), \(M\), and \(X\) from an unordered set of \((x,y)\) points.

The curve is defined by the parametric equations

$$
\begin{align*}
x(t) &= t\cos\theta - e^{M|t|}\sin(0.3t)\sin\theta + X \\
y(t) &= 42 + t\sin\theta + e^{M|t|}\sin(0.3t)\cos\theta
\end{align*}
$$

with domain

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

## Final Answer

| Parameter  |                                   Value | Constraint                      |
| ---------- | --------------------------------------: | ------------------------------- |
| \(\theta\) | **30.0000°** (\(\approx 0.523599\) rad) | \(0^\circ < \theta < 50^\circ\) |
| \(M\)      |                              **0.0300** | \(-0.05 < M < 0.05\)            |
| \(X\)      |                             **55.0000** | \(0 < X < 100\)                 |

### Desmos

Copy-paste ready:

```text
\left(t*\cos(0.523599)-e^{0.03\left|t\right|}\cdot\sin(0.3t)\sin(0.523599)+55,42+t*\sin(0.523599)+e^{0.03\left|t\right|}\cdot\sin(0.3t)\cos(0.523599)\right)
```

[Open fitted curve in Desmos](https://www.desmos.com/calculator/jjcufejdax)

---

# Mathematical Approach

## 1. Rewrite as a Rigid Motion

Define

$$
A(t)=t,
\qquad
B(t)=e^{M|t|}\sin(0.3t).
$$

The parametric equations can then be written as a 2-D rotation followed by a translation:

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

Because a rotation matrix is orthogonal, its inverse is simply its transpose:

$$
\begin{pmatrix}
u\\
v
\end{pmatrix}
=
\begin{pmatrix}
\cos\theta & \sin\theta\\
-\sin\theta & \cos\theta
\end{pmatrix}
\begin{pmatrix}
x-X\\
y-42
\end{pmatrix}.
$$

Therefore,

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

For the correct parameters,

$$
B_{\text{est}}
\approx
e^{M|t_{\text{est}}|}\sin(0.3t_{\text{est}}).
$$

This is the key simplification in the solution: the unknown correspondence between the unordered \((x,y)\) points and their parameter values \(t\) does not need to be searched explicitly. For any candidate \((\theta,M,X)\), the corresponding \(t\) can be recovered directly by the inverse rotation.

Because the transformation is a rigid rotation and translation, the residual in the rotated coordinate system preserves Euclidean distance to the corresponding point on the curve.

---

## 2. Optimisation

The problem is reduced to a bounded nonlinear least-squares problem in only three unknowns:

$$
(\theta,M,X).
$$

For each candidate parameter set:

1. Recover \(t_{\text{est}}\) and \(B_{\text{est}}\) using the inverse rotation.
2. Evaluate the model

$$
B_{\text{model}}
=
e^{M|t_{\text{est}}|}\sin(0.3t_{\text{est}}).
$$

3. Form the residual

$$
r =
B_{\text{est}}-B_{\text{model}}.
$$

4. Add soft quadratic penalties when recovered \(t\) values move outside the allowed interval \((6,60)\).
5. Run multiple bounded nonlinear least-squares optimisations from random initialisations.
6. Cross-check the best solution using a global Differential Evolution search.

### Multi-start optimisation

The implementation performs **60 independent bounded nonlinear least-squares runs** using Trust-Region Reflective optimisation, with starting points sampled uniformly within the allowed parameter ranges.

The purpose of the multi-start approach is to reduce the chance of accepting a poor local minimum.

The best solution is then independently cross-checked using Differential Evolution.

All runs converge to the same parameter basin:

```text
theta ≈ 30°
M     ≈ 0.03
X     ≈ 55
```

---

# Validation

The assignment evaluates the quality of the recovered curve using an L1 distance between uniformly sampled points on the expected and predicted curves.

For the validation step, the observed points are first ordered using their recovered \(t_{\text{est}}\) values and interpolated onto a dense uniform \(t\)-grid.

The predicted curve is evaluated analytically on the same grid.

The point-wise L1 error is then

$$
L_1(t)
=
|x_{\text{pred}}(t)-x_{\text{obs}}(t)|
+
|y_{\text{pred}}(t)-y_{\text{obs}}(t)|.
$$

## Validation Results

Results on the supplied **1,500-point dataset**:

| Metric               |                  Value |
| -------------------- | ---------------------: |
| Mean uniform-grid L1 | \(1.745\times10^{-4}\) |
| 95th-percentile L1   | \(5.914\times10^{-4}\) |
| Maximum L1           | \(9.040\times10^{-3}\) |

These errors are very small relative to the scale of the curve and indicate that the recovered parameters reproduce the supplied data closely across the full domain.

---

# Results

## Fitted Curve

The recovered parameters are compared visually against the supplied data points.

![Fitted Curve](results/fit_plot.png)

---

## Uniform-grid L1 Residual versus \(t\)

The following plot shows how the L1 reconstruction error varies across the parameter domain.

![Uniform-grid L1 Residual versus t](results/residual_vs_t.png)

---

## Local Parameter Sensitivity

This plot shows how the fitted solution changes when the recovered parameters are perturbed locally.

![Local Parameter Sensitivity](results/sensitivity.png)

---

# Reproduction

Clone the repository and install the required dependencies:

```bash
git clone https://github.com/Prianshu7296/Flam-RnD-Assignment.git
cd Flam-RnD-Assignment

pip install -r requirements.txt
```

Run the fitting pipeline:

```bash
python src/fit_curve.py
```

This recovers \(\theta\), \(M\), and \(X\).

Run validation:

```bash
python src/verify_fit.py
```

This evaluates the assignment-style L1 reconstruction error.

Regenerate the plots:

```bash
python src/plot_fit.py
```

---

# Repository Structure

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
│   ├── fit_curve.py      # Main parameter optimisation
│   ├── verify_fit.py     # L1 validation
│   └── plot_fit.py       # Visualisation
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

# Notes

* The supplied points are unordered.
* The original \(t\) values are not required as input.
* The inverse-rotation step eliminates the need for a point-to-curve correspondence search during optimisation.
* The fitting problem is therefore reduced to only three unknown parameters: \(\theta\), \(M\), and \(X\).
* Multiple local optimisation runs are used to reduce sensitivity to initialisation.
* Differential Evolution provides an independent global-optimisation cross-check.
* The final parameters are:

```text
theta = 30°
M     = 0.03
X     = 55
```

These values reproduce the supplied curve with a very small reconstruction error.

---

# Conclusion

The final recovered parametric curve is

$$
\boxed{\theta=30^\circ,\quad M=0.03,\quad X=55}
$$

The main strength of the approach is the analytical inverse-rotation step, which converts an unordered point-cloud fitting problem into a low-dimensional optimisation problem while avoiding explicit point-to-curve correspondence search.


