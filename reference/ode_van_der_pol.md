# Van der Pol oscillator

The Van der Pol oscillator, a nonlinear oscillator with self-sustaining
oscillations (a limit cycle), written as a 2D first-order system:

## Usage

``` r
ode_van_der_pol(t, y, parameters = c(mu = 1))
```

## Arguments

- t:

  Numeric scalar. Time (autonomous; included for deSolve compatibility).

- y:

  Numeric vector of length 2. State vector: \\y\[1\]\\ = position \\x\\,
  \\y\[2\]\\ = velocity.

- parameters:

  Named numeric vector with element:

  `mu`

  :   Nonlinearity parameter. Default: `1`. Must be \\\geq 0\\.

## Value

A list with one element: a numeric vector of length 2 containing
\\(dx/dt, dy/dt)\\, as required by
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html).

## Details

\$\$\frac{dx}{dt} = y\$\$ \$\$\frac{dy}{dt} = \mu (1 - x^2) y - x\$\$

where \\x(t)\\ is the position (displacement), \\y(t)\\ is the velocity,
and \\\mu \geq 0\\ is the nonlinearity/damping parameter.

When \\\mu = 0\\ the system reduces to a harmonic oscillator with
circular orbits. For \\\mu \> 0\\ a stable limit cycle exists; the shape
becomes increasingly relaxation-oscillator-like as \\\mu\\ increases.
The origin \\(0, 0)\\ is an unstable spiral for \\\mu \> 0\\.

## Examples

``` r
# At (x=2, y=0) with mu=1: dx/dt=0, dy/dt = 1*(1-4)*0 - 2 = -2
ode_van_der_pol(t = 0, y = c(2, 0), parameters = c(mu = 1))
#> [[1]]
#> [1]  0 -2
#> 

if (FALSE) { # \dontrun{
gg_flow_field(ode_van_der_pol, xlim = c(-4, 4), ylim = c(-4, 4))
} # }
```
