# Example ODE system 2 (1D)

The derivative function of example one-dimensional ODE system 2 from
Grayling (2014):

## Usage

``` r
ode_example_02(t, y, parameters = NULL)
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

\$\$\frac{dy}{dt} = y(1 - y)(2 - y)\$\$

Three equilibria: \\y^\* = 0\\ (unstable), \\y^\* = 1\\ (stable), \\y^\*
= 2\\ (unstable).

## References

Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
[doi:10.32614/RJ-2014-023](https://doi.org/10.32614/RJ-2014-023)

## Examples

``` r
ode_example_02(t = 0, y = c(0.5), parameters = NULL)
#> [[1]]
#> [1] 0.375
#> 

if (FALSE) { # \dontrun{
gg_phase_portrait(ode_example_02, ylim = c(-0.5, 2.5))
} # }
```
