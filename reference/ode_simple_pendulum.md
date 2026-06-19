# Simple pendulum model

The simple pendulum written as a 2D first-order system:

## Usage

``` r
ode_simple_pendulum(t, y, parameters = c(gL = 1, b = 0))
```

## Arguments

- t:

  Numeric scalar. Time (autonomous; included for deSolve compatibility).

- y:

  Numeric vector of length 2. State vector: \\y\[1\]\\ = angle
  \\\theta\\ (radians), \\y\[2\]\\ = angular velocity \\\omega\\.

- parameters:

  Named numeric vector with elements:

  `gL`

  :   Ratio \\g/L\\ (gravitational acceleration divided by pendulum
      length). Default: `1`.

  `b`

  :   Damping coefficient. Default: `0` (undamped).

## Value

A list with one element: a numeric vector of length 2 containing
\\(d\theta/dt, d\omega/dt)\\, as required by
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html).

## Details

\$\$\frac{d\theta}{dt} = \omega\$\$ \$\$\frac{d\omega}{dt} =
-\frac{g}{L} \sin(\theta) - b\\\omega\$\$

where \\\theta\\ is the angle from vertical (radians), \\\omega\\ is the
angular velocity, \\g/L\\ is the ratio of gravitational acceleration to
pendulum length (combined into a single parameter `gL` for convenience),
and \\b \geq 0\\ is a damping coefficient.

Setting \\b = 0\\ gives the undamped (conservative) pendulum with
heteroclinic orbits connecting the unstable equilibria at \\\theta =
\pm\pi\\. Adding damping (\\b \> 0\\) makes the stable equilibrium at
\\\theta = 0\\ a stable spiral (for small \\b\\) or stable node (for
large \\b\\).

## Examples

``` r
# Undamped pendulum at theta = pi/4, omega = 0
ode_simple_pendulum(t = 0, y = c(pi/4, 0), parameters = c(gL = 1, b = 0))
#> [[1]]
#> [1]  0.0000000 -0.7071068
#> 

# Damped pendulum
ode_simple_pendulum(t = 0, y = c(pi/4, 0), parameters = c(gL = 1, b = 0.5))
#> [[1]]
#> [1]  0.0000000 -0.7071068
#> 

if (FALSE) { # \dontrun{
gg_flow_field(ode_simple_pendulum,
              xlim = c(-pi, pi), ylim = c(-3, 3))
} # }
```
