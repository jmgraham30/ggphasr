# Monomolecular growth model

The ODE for monomolecular (saturating, or Mitscherlich) growth:

## Usage

``` r
ode_monomolecular(t, y, parameters = c(r = 1, K = 10))
```

## Arguments

- t:

  Numeric scalar. Time (not used directly, included for deSolve
  compatibility).

- y:

  Numeric vector of length 1. Current state: \\y\[1\]\\ is the current
  value.

- parameters:

  Named numeric vector with elements:

  `r`

  :   Rate constant. Default: `1`.

  `K`

  :   Asymptote (maximum value). Default: `10`.

## Value

A list with one element: a numeric vector of length 1 containing
\\dy/dt\\, as required by
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html).

## Details

\$\$\frac{dy}{dt} = r \left(K - y\right)\$\$

where \\y(t)\\ is the quantity of interest (e.g., biomass), \\r \> 0\\
is the rate constant, and \\K \> 0\\ is the asymptote (maximum
attainable value). Growth decelerates monotonically as \\y\\ approaches
\\K\\. There is a single stable equilibrium at \\y = K\\.

Unlike logistic growth, there is no inflection point — growth rate is
greatest at \\y = 0\\ and decreases linearly to zero at \\y = K\\.

The analytic solution is \\y(t) = K(1 - e^{-rt}) + y_0 \\ e^{-rt}\\.

## See also

[`ode_logistic()`](https://jmgraham30.github.io/ggphasr/reference/ode_logistic.md),
[`ode_von_bertalanffy()`](https://jmgraham30.github.io/ggphasr/reference/ode_von_bertalanffy.md)

## Examples

``` r
# Evaluate at y = 4 with default parameters (r = 1, K = 10)
# dy/dt = 1 * (10 - 4) = 6
ode_monomolecular(t = 0, y = c(4), parameters = c(r = 1, K = 10))
#> [[1]]
#> [1] 6
#> 

if (FALSE) { # \dontrun{
gg_phase_portrait(ode_monomolecular, ylim = c(0, 12))
} # }
```
