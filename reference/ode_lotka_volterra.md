# Lotka-Volterra predator-prey model

The classic two-species predator-prey system:

## Usage

``` r
ode_lotka_volterra(
  t,
  y,
  parameters = c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1)
)
```

## Arguments

- t:

  Numeric scalar. Time (autonomous system; included for deSolve
  compatibility).

- y:

  Numeric vector of length 2. State vector: \\y\[1\]\\ = prey (\\x\\),
  \\y\[2\]\\ = predator.

- parameters:

  Named numeric vector with elements:

  `alpha`

  :   Prey growth rate. Default: `1`.

  `beta`

  :   Predation rate. Default: `0.5`.

  `delta`

  :   Predator growth rate per prey consumed. Default: `0.5`.

  `gamma`

  :   Predator death rate. Default: `1`.

## Value

A list with one element: a numeric vector of length 2 containing
\\(dx/dt, dy/dt)\\, as required by
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html).

## Details

\$\$\frac{dx}{dt} = \alpha x - \beta x y\$\$ \$\$\frac{dy}{dt} = \delta
x y - \gamma y\$\$

where \\x(t)\\ is prey abundance, \\y(t)\\ is predator abundance,
\\\alpha \> 0\\ is the prey growth rate, \\\beta \> 0\\ is the predation
rate, \\\delta \> 0\\ is the predator growth rate per prey consumed, and
\\\gamma \> 0\\ is the predator death rate.

The system has two equilibria: a trivial unstable equilibrium at \\(0,
0)\\ and a neutrally stable center at \\(x^\*, y^\*) = (\gamma/\delta,\\
\alpha/\beta)\\. All trajectories starting in the positive quadrant are
closed orbits around the interior equilibrium.

## See also

[`ode_competition()`](https://jmgraham30.github.io/ggphasr/reference/ode_competition.md)

## Examples

``` r
# Interior equilibrium at (gamma/delta, alpha/beta) = (2, 2)
ode_lotka_volterra(t = 0, y = c(2, 2),
                  parameters = c(alpha = 1, beta = 0.5,
                                 delta = 0.5, gamma = 1))
#> [[1]]
#> [1] 0 0
#> 

if (FALSE) { # \dontrun{
gg_flow_field(ode_lotka_volterra, xlim = c(0, 5), ylim = c(0, 5))
} # }
```
