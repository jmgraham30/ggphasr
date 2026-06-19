# Example ODE system 12 (2D)

The derivative function of example two-dimensional ODE system 12 from
Grayling (2014):

## Usage

``` r
ode_example_12(t, y, parameters = NULL)
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

\$\$\frac{dx}{dt} = y + x(x^2 + y^2 - 1), \qquad \frac{dy}{dt} = -x +
y(x^2 + y^2 - 1)\$\$

An unstable spiral at \\(0,0)\\ and a second equilibrium at \\(1,1)\\, a
saddle point. Useful for illustrating unstable equilibria and limit
cycle-like behavior near the unit circle.

## References

Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
[doi:10.32614/RJ-2014-023](https://doi.org/10.32614/RJ-2014-023)

## Examples

``` r
ode_example_12(t = 0, y = c(2, 2), parameters = NULL)
#> [[1]]
#> [1] 16 12
#> 

if (FALSE) { # \dontrun{
gg_flow_field(ode_example_12, xlim = c(-4, 4), ylim = c(-4, 4))
} # }
```
