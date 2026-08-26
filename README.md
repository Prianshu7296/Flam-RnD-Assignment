# Parametric Curve Fitting — R&D / AI Assignment

Recover the unknown parameters `θ`, `M`, and `X` from an unordered set of `(x, y)` points.

The curve is:

```text
x(t) = t·cos(θ) − e^(M|t|)·sin(0.3t)·sin(θ) + X

y(t) = 42 + t·sin(θ) + e^(M|t|)·sin(0.3t)·cos(θ)

6 < t < 60
```

## Final Answer

| Parameter |                       Value | Given constraint |
| --------- | --------------------------: | ---------------- |
| `θ`       | **30.0000°** (0.523599 rad) | 0° < θ < 50°     |
| `M`       |                  **0.0300** | −0.05 < M < 0.05 |
| `X`       |                 **55.0000** | 0 < X < 100      |

These are the values I get from the fit.

### Desmos

The fitted curve is available here:

https://www.desmos.com/calculator/jjcufejdax

The same equation can also be pasted into Desmos manually:

```text
\left(
t\cos(0.5236)-e^{0.0300|t|}\sin(0.3t)\sin(0.5236)+55,
42+t\sin(0.5236)+e^{0.0300|t|}\sin(0.3t)\cos(0.5236)
\right)
```

with the domain:

```text
6 ≤ t ≤ 60
```

The fitted curve overlaps the supplied points very closely.

## How I approached it

### 1. Rewrite the equation

At first I looked at the two equations together rather than fitting `x(t)` and `y(t)` separately.

Let

```text
A(t) = t

B(t) = e^(M|t|) · sin(0.3t)
```

Then the equations become

```text
x - X = A cos(θ) - B sin(θ)

y - 42 = A sin(θ) + B cos(θ)
```

This is just a rotation of the `(A, B)` coordinates by `θ`, followed by a translation.

In matrix form:

```text
[ x - X ]   [ cos(θ)  -sin(θ) ] [ A ]
[ y - 42 ] = [ sin(θ)   cos(θ) ] [ B ]
```

The useful part here is that a rotation can be inverted exactly.

### 2. Recover t from each point

For a candidate `θ` and `X`, I can rotate each point back:

```text
t_est = (x - X) cos(θ) + (y - 42) sin(θ)

B_est = -(x - X) sin(θ) + (y - 42) cos(θ)
```

So I don't need to know which `t` belongs to each point beforehand.

For the correct values of `θ` and `X`, `t_est` should be the original parameter value, and `B_est` should follow:

```text
B_est ≈ exp(M|t_est|) · sin(0.3t_est)
```

This reduces the problem to estimating only three parameters: `θ`, `M`, and `X`.

### 3. Fit the three parameters

I used bounded nonlinear least squares for the fit.

For each candidate `(θ, M, X)` I:

1. Rotate all the points back.
2. Calculate `t_est` and `B_est`.
3. Compare `B_est` with `exp(M|t_est|) sin(0.3t_est)`.
4. Penalize points whose estimated `t` goes outside the allowed range.
5. Minimize the resulting residual.

Because the objective can have more than one local minimum, I ran the optimizer from multiple random starting points instead of relying on one initialization. I used 60 starts and kept the best result.

The runs mostly converged to the same parameter values, which gave me some confidence that the solution was not just one lucky local minimum.

### 4. Check the result

The final parameters are:

```text
θ = 30°
M = 0.03
X = 55
```

I then reconstructed the curve using these values and checked it against the supplied data.

The fit is very close. Using the evaluation in the repository gives:

```text
Mean L1 distance: 0.0041
Max L1 distance:  0.0124
95th percentile:  0.0081
```

I also generated the overlay plot in `results/fit_plot.png` as a visual check.

The important thing for me was that both the numerical fit and the plotted curve gave the same conclusion: `30°, 0.03, 55` is the parameter set that fits the supplied points.

## Why this works

The main simplification is the rotation.

Instead of treating this as a completely unknown parametric curve, I can rotate the points back into the coordinate system where one coordinate is simply `t`.

That means the correspondence problem becomes much easier. Once `θ` and `X` are close to the correct values, the recovered `t` values have the expected structure and the remaining parameter `M` controls the amplitude of the oscillation.

## Repository structure

```text
curve-param-fitting/

├── README.md
├── requirements.txt
├── .gitignore

├── data/
│   └── xy_data.csv

├── src/
│   ├── fit_curve.py
│   ├── verify_fit.py
│   └── plot_fit.py

└── results/
    ├── fit_result.txt
    ├── desmos_equation.txt
    └── fit_plot.png
```

## How to reproduce

```bash
git clone https://github.com/Prianshu7296/Flam-RnD-Assignment.git
cd Flam-RnD-Assignment

pip install -r requirements.txt

python src/fit_curve.py
python src/verify_fit.py
python src/plot_fit.py
```

The first script fits the parameters, the second checks the saved result, and the third generates the plot.

## Notes

The `t` values are not part of the required output, so I only use them internally during fitting.

The final answer required by the assignment is:

```text
θ = 30°
M = 0.03
X = 55
```

The Desmos version of the fitted curve is here:

https://www.desmos.com/calculator/jjcufejdax
