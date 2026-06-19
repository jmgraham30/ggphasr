# Logistic growth model

The ODE for logistic (Verhulst) population growth:

## Usage

``` r
ode_logistic(t, y, parameters = c(r = 1, K = 10))
```

## Arguments

- t:

  Numeric scalar. Time (not used directly, included for deSolve
  compatibility).

- y:

  Numeric vector of length 1. Current state: \\y\[1\]\\ is the
  population size.

- parameters:

  Named numeric vector with elements:

  `r`

  :   Intrinsic growth rate. Default: `1`.

  `K`

  :   Carrying capacity. Default: `10`.

## Value

A list with one element: a numeric vector of length 1 containing
\\dy/dt\\, as required by
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html).

## Details

\$\$\frac{dy}{dt} = r \\ y \left(1 - \frac{y}{K}\right)\$\$

where \\y(t)\\ is the population size, \\r \> 0\\ is the intrinsic
growth rate, and \\K \> 0\\ is the carrying capacity. The model has two
equilibria: an unstable equilibrium at \\y = 0\\ and a stable
equilibrium at \\y = K\\.

The analytic solution is \\y(t) = K / (1 + ((K - y_0)/y_0) \\
e^{-rt})\\.

## See also

[`ode_exponential()`](https://jmgraham30.github.io/ggphasr/reference/ode_exponential.md),
[`ode_monomolecular()`](https://jmgraham30.github.io/ggphasr/reference/ode_monomolecular.md)

## Examples

``` r
# Evaluate at y = 5 with default parameters (r = 1, K = 10)
# dy/dt = 1 * 5 * (1 - 5/10) = 2.5
ode_logistic(t = 0, y = c(5), parameters = c(r = 1, K = 10))
#> [[1]]
#> [1] 2.5
#> 

if (FALSE) { # \dontrun{
gg_phase_portrait(ode_logistic, ylim = c(-2, 12))
} # }
```
