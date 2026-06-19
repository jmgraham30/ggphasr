# Exponential growth model

The ODE for exponential (Malthusian) population growth:

## Usage

``` r
ode_exponential(t, y, parameters = c(r = 0.5))
```

## Arguments

- t:

  Numeric scalar. Time (not used directly, included for deSolve
  compatibility).

- y:

  Numeric vector of length 1. Current state: \\y\[1\]\\ is the
  population size.

- parameters:

  Named numeric vector with element:

  `r`

  :   Intrinsic growth rate. Default: `0.5`.

## Value

A list with one element: a numeric vector of length 1 containing
\\dy/dt\\, as required by
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html).

## Details

\$\$\frac{dy}{dt} = r \\ y\$\$

where \\y(t)\\ is the population size (or any quantity growing
proportionally to itself) and \\r\\ is the intrinsic growth rate. When
\\r \> 0\\ the population grows without bound; when \\r \< 0\\ it decays
to zero; when \\r = 0\\ it is constant.

The analytic solution is \\y(t) = y_0 \\ e^{rt}\\.

## See also

[`ode_logistic()`](https://jmgraham30.github.io/ggphasr/reference/ode_logistic.md),
[`ode_monomolecular()`](https://jmgraham30.github.io/ggphasr/reference/ode_monomolecular.md)

## Examples

``` r
# Evaluate at y = 2 with default parameters (r = 0.5)
# dy/dt = 0.5 * 2 = 1
ode_exponential(t = 0, y = c(2), parameters = c(r = 0.5))
#> [[1]]
#> [1] 1
#> 

# Use with gg_flow_field() for a 1D phase portrait
if (FALSE) { # \dontrun{
gg_phase_portrait(ode_exponential, ylim = c(-3, 3))
} # }
```
