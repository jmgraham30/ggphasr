# Example ODE system 13 (2D)

The derivative function of example two-dimensional ODE system 13 from
Grayling (2014):

## Usage

``` r
ode_example_13(t, y, parameters = NULL)
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

\$\$\frac{dx}{dt} = x(1-x) - xy, \qquad \frac{dy}{dt} =
y\\\left(\frac{x}{x+0.5} - 0.5\right)\$\$

A predator-prey system with a saturating functional response (Holling
type II). Equilibria include the trivial \\(0,0)\\, the prey-only state
\\(1,0)\\, and a coexistence equilibrium that may be stable or unstable
depending on parameters, giving rise to a limit cycle via Hopf
bifurcation.

## References

Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
[doi:10.32614/RJ-2014-023](https://doi.org/10.32614/RJ-2014-023)

## Examples

``` r
ode_example_13(t = 0, y = c(0.5, 0.5), parameters = NULL)
#> [[1]]
#> [1] 0 0
#> 

if (FALSE) { # \dontrun{
gg_flow_field(ode_example_13, xlim = c(0, 1.5), ylim = c(0, 1))
} # }
```
