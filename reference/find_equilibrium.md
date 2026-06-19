# Find equilibria of an ODE system numerically

Locates one or more equilibrium points of a one- or two-dimensional
autonomous ODE system using Newton-Raphson root-finding via
[`rootSolve::multiroot()`](https://rdrr.io/pkg/rootSolve/man/multiroot.html).
Accepts either a single initial guess or a grid of starting points to
search for multiple equilibria.

## Usage

``` r
find_equilibrium(
  deriv,
  y0 = NULL,
  system = c("two.dim", "one.dim"),
  parameters = NULL,
  xlim = NULL,
  ylim = NULL,
  n_grid = 10L,
  tol = 1e-08,
  dedup_tol = 1e-04
)
```

## Arguments

- deriv:

  A function describing the ODE system, in Convention A or B. See
  [ggphasr-package](https://jmgraham30.github.io/ggphasr/reference/ggphasr-package.md)
  for details.

- y0:

  Initial guess(es) for the root-finder. One of:

  - A numeric vector of length 1 (1D) or 2 (2D): a single starting
    point.

  - A numeric matrix with one row per starting point.

  - A list of numeric vectors.

  - `NULL`: triggers automatic grid search over `xlim` × `ylim`
    (requires `xlim` and `ylim` to be supplied).

- system:

  Character: `"two.dim"` (default) or `"one.dim"`.

- parameters:

  Parameter vector or list passed to `deriv`.

- xlim:

  Numeric vector of length 2. x-axis search range. Required when
  `y0 = NULL`. Default: `NULL`.

- ylim:

  Numeric vector of length 2. y-axis search range. Required when
  `y0 = NULL`. Default: `NULL`.

- n_grid:

  Integer. Grid resolution per axis for automatic search when
  `y0 = NULL`. Default: `10` (giving up to 100 starting points for 2D
  systems).

- tol:

  Numeric. Convergence tolerance passed to
  [`rootSolve::multiroot()`](https://rdrr.io/pkg/rootSolve/man/multiroot.html).
  Default: `1e-8`.

- dedup_tol:

  Numeric. Distance tolerance for removing duplicate roots. Two
  equilibria closer than this are treated as the same point. Default:
  `1e-4`.

## Value

A list of numeric vectors, each giving the coordinates of one
equilibrium point. If only one equilibrium is found, a list of length 1
is returned (not an unwrapped vector) for consistency.

## Details

Newton-Raphson is sensitive to the initial guess — it may fail to
converge or converge to a non-equilibrium point for badly chosen
starting values. When `y0 = NULL`, the function tries `n_grid^2`
starting points on a regular grid, collects all converged roots, and
deduplicates them. This is more reliable than a single guess but still
not exhaustive.

Failed convergences are silently discarded; only roots where the
residual is below `tol * 100` are retained.

## See also

[`classify_equilibrium()`](https://jmgraham30.github.io/ggphasr/reference/classify_equilibrium.md)

## Examples

``` r
# Single guess: find the interior equilibrium of Lotka-Volterra
find_equilibrium(
  ode_lotka_volterra,
  y0         = c(1.5, 1.5),
  parameters = c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1)
)
#> [[1]]
#> [1] 2 2
#> 

# Automatic grid search: find all equilibria of example 11
find_equilibrium(
  ode_example_11,
  y0    = NULL,
  xlim  = c(0, 4),
  ylim  = c(0, 4)
)
#> [[1]]
#> [1] -9.524216e-10  2.000000e+00
#> 
#> [[2]]
#> [1] 0 0
#> 
#> [[3]]
#> [1] 1 1
#> 
#> [[4]]
#> [1] 3 0
#> 

# 1D system: find equilibria of the logistic equation
find_equilibrium(
  ode_logistic,
  y0         = NULL,
  system     = "one.dim",
  ylim       = c(-1, 12),
  parameters = c(r = 1, K = 10)
)
#> [[1]]
#> [1] -1.709727e-16
#> 
#> [[2]]
#> [1] 10
#> 
```
