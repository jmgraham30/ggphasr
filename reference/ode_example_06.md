# Example ODE system 6 (2D)

The derivative function of example two-dimensional ODE system 6 from
Grayling (2014):

## Usage

``` r
ode_example_06(t, y, parameters = NULL)
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

\$\$\frac{dx}{dt} = y, \qquad \frac{dy}{dt} = -x\$\$

The harmonic oscillator. All trajectories are closed ellipses (circles
when plotted on equal-aspect axes) centered on the origin, which is a
neutrally stable center.

## References

Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
[doi:10.32614/RJ-2014-023](https://doi.org/10.32614/RJ-2014-023)

## Examples

``` r
ode_example_06(t = 0, y = c(1, 0), parameters = NULL)
#> [[1]]
#> [1]  0 -1
#> 

if (FALSE) { # \dontrun{
gg_flow_field(ode_example_06, xlim = c(-3, 3), ylim = c(-3, 3))
} # }
```
