# Example ODE system 9 (2D)

The derivative function of example two-dimensional ODE system 9 from
Grayling (2014):

## Usage

``` r
ode_example_09(t, y, parameters = NULL)
```

## Arguments

- t:

  Numeric scalar. Time (autonomous; included for deSolve compatibility).

- y:

  Numeric vector of length 2. State vector: \\y\[1\] = x\\, \\y\[2\] =
  y\\.

- parameters:

  Not used. Accepted for deSolve compatibility.

## Value

A list with one element: a numeric vector of length 2 containing
\\(dx/dt, dy/dt)\\, as required by
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html).

## Details

\$\$\frac{dx}{dt} = -x + y, \qquad \frac{dy}{dt} = -x - y\$\$

A stable spiral at the origin (eigenvalues \\-1 \pm i\\). Trajectories
spiral inward to \\(0, 0)\\.

## References

Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
[doi:10.32614/RJ-2014-023](https://doi.org/10.32614/RJ-2014-023)

## Examples

``` r
ode_example_09(t = 0, y = c(2, 0), parameters = NULL)
#> [[1]]
#> [1] -2 -2
#> 

if (FALSE) { # \dontrun{
gg_flow_field(ode_example_09, xlim = c(-3, 3), ylim = c(-3, 3))
} # }
```
