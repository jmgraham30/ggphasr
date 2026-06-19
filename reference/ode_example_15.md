# Example ODE system 15 (2D)

The derivative function of example two-dimensional ODE system 15 from
Grayling (2014):

## Usage

``` r
ode_example_15(t, y, parameters = NULL)
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

\$\$\frac{dx}{dt} = x - y - x(x^2+y^2), \qquad \frac{dy}{dt} = x + y -
y(x^2+y^2)\$\$

In polar coordinates this becomes \\dr/dt = r(1-r^2)\\, \\d\theta/dt =
1\\, revealing a stable limit cycle at \\r = 1\\ (the unit circle) and
an unstable spiral at the origin. A classic example for illustrating
limit cycles.

## References

Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
[doi:10.32614/RJ-2014-023](https://doi.org/10.32614/RJ-2014-023)

## Examples

``` r
ode_example_15(t = 0, y = c(0.5, 0), parameters = NULL)
#> [[1]]
#> [1] 0.375 0.500
#> 

if (FALSE) { # \dontrun{
gg_flow_field(ode_example_15, xlim = c(-2, 2), ylim = c(-2, 2))
} # }
```
