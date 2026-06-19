# Example ODE system 1 (1D)

The derivative function of example one-dimensional ODE system 1 from
Grayling (2014):

## Usage

``` r
ode_example_01(t, y, parameters = NULL)
```

## Arguments

- t:

  Numeric scalar. Time (autonomous; included for deSolve compatibility).

- y:

  Numeric vector of length 1. Current state.

- parameters:

  Not used. Accepted for deSolve compatibility.

## Value

A list with one element: a numeric vector of length 1 containing
\\dy/dt\\, as required by
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html).

## Details

\$\$\frac{dy}{dt} = 4 - y^2\$\$

Equilibria at \\y^\* = \pm 2\\: \\y^\* = 2\\ is stable, \\y^\* = -2\\ is
unstable.

## References

Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
[doi:10.32614/RJ-2014-023](https://doi.org/10.32614/RJ-2014-023)

## Examples

``` r
ode_example_01(t = 0, y = c(1), parameters = NULL)
#> [[1]]
#> [1] 3
#> 

if (FALSE) { # \dontrun{
gg_phase_portrait(ode_example_01, ylim = c(-4, 4))
} # }
```
