# Classify an equilibrium point of an ODE system

Computes the Jacobian matrix of a one- or two-dimensional ODE system at
a known equilibrium point using forward finite differences, then
classifies the equilibrium using eigenvalue analysis. Returns a tidy
data frame row that can be combined with results from multiple
equilibria using [`rbind()`](https://rdrr.io/r/base/cbind.html).

## Usage

``` r
classify_equilibrium(
  deriv,
  equilibrium,
  system = c("two.dim", "one.dim"),
  parameters = NULL,
  h = 1e-06,
  tol = 1e-08
)
```

## Arguments

- deriv:

  A function describing the ODE system, in Convention A or B.

- equilibrium:

  Numeric vector of length 1 (1D) or 2 (2D) giving the equilibrium
  coordinates. Typically obtained from
  [`find_equilibrium()`](https://jmgraham30.github.io/ggphasr/reference/find_equilibrium.md).

- system:

  Character: `"two.dim"` (default) or `"one.dim"`.

- parameters:

  Parameter vector or list passed to `deriv`.

- h:

  Numeric. Step size for finite-difference Jacobian computation.
  Default: `1e-6`.

- tol:

  Numeric. Tolerance for treating eigenvalue real/imaginary parts as
  zero. Default: `1e-8`.

## Value

A [`data.frame()`](https://rdrr.io/r/base/data.frame.html) with one row
and the following columns:

- `x`:

  x-coordinate of the equilibrium (2D) or `NA` (1D).

- `y`:

  y-coordinate of the equilibrium (both 1D and 2D).

- `classification`:

  Character. One of: `"Stable node"`, `"Unstable node"`,
  `"Stable spiral"`, `"Unstable spiral"`, `"Center"`, `"Saddle"`,
  `"Non-isolated equilibrium"` (2D); or `"Stable"`, `"Unstable"`,
  `"Inconclusive (df/dy = 0)"` (1D).

- `tr`:

  Trace of the Jacobian (2D only; `NA` for 1D).

- `det`:

  Determinant of the Jacobian (2D only; `NA` for 1D).

- `jacobian_11`:

  `J[1,1]` — always present.

- `jacobian_12`:

  `J[1,2]` (2D only; `NA` for 1D).

- `jacobian_21`:

  `J[2,1]` (2D only; `NA` for 1D).

- `jacobian_22`:

  `J[2,2]` (2D only; `NA` for 1D).

- `lambda_1_re`:

  Real part of eigenvalue 1 (2D only; `NA` for 1D).

- `lambda_1_im`:

  Imaginary part of eigenvalue 1 (2D only).

- `lambda_2_re`:

  Real part of eigenvalue 2 (2D only).

- `lambda_2_im`:

  Imaginary part of eigenvalue 2 (2D only).

## Details

Results from multiple equilibria can be combined into a summary table:

    eq_list <- find_equilibrium(ode_example_11, y0 = NULL,
                                 xlim = c(0,4), ylim = c(0,4))
    results <- do.call(rbind, lapply(eq_list, function(eq) {
      classify_equilibrium(ode_example_11, equilibrium = eq)
    }))

## See also

[`find_equilibrium()`](https://jmgraham30.github.io/ggphasr/reference/find_equilibrium.md)

## Examples

``` r
# Classify the interior equilibrium of Lotka-Volterra (expect: Center)
classify_equilibrium(
  ode_lotka_volterra,
  equilibrium = c(2, 2),
  parameters  = c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1)
)
#>   x y classification tr det jacobian_11 jacobian_12 jacobian_21 jacobian_22
#> 1 2 2         Center  0   1           0          -1           1           0
#>   lambda_1_re lambda_1_im lambda_2_re lambda_2_im
#> 1           0           1           0          -1

# Classify all equilibria of example 11
eq_list <- find_equilibrium(ode_example_11, y0 = NULL,
                             xlim = c(0, 4), ylim = c(0, 4))
do.call(rbind, lapply(eq_list, function(eq) {
  classify_equilibrium(ode_example_11, equilibrium = eq)
}))
#>               x y classification        tr       det jacobian_11   jacobian_12
#> 1 -9.524185e-10 2    Stable node -3.000002  2.000003   -1.000001  1.904837e-09
#> 2  0.000000e+00 0  Unstable node  4.999998  5.999995    2.999999  0.000000e+00
#> 3  1.000000e+00 1         Saddle -2.000002 -0.999998   -1.000001 -2.000000e+00
#> 4  3.000000e+00 0    Stable node -4.000002  3.000004   -3.000001 -6.000000e+00
#>   jacobian_21 jacobian_22 lambda_1_re lambda_1_im lambda_2_re lambda_2_im
#> 1          -2   -2.000001   -2.000001           0  -1.0000010           0
#> 2           0    1.999999    2.999999           0   1.9999990           0
#> 3          -1   -1.000001   -2.414215           0   0.4142126           0
#> 4           0   -1.000001   -3.000001           0  -1.0000010           0
```
