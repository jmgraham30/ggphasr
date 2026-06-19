# Two-species Lotka-Volterra competition model

The classic two-species interspecific competition model:

## Usage

``` r
ode_competition(
  t,
  y,
  parameters = c(r1 = 1, r2 = 1, K1 = 10, K2 = 10, a12 = 0.5, a21 = 0.5)
)
```

## Arguments

- t:

  Numeric scalar. Time (autonomous; included for deSolve compatibility).

- y:

  Numeric vector of length 2. State vector: \\y\[1\]\\ = \\N_1\\,
  \\y\[2\]\\ = \\N_2\\.

- parameters:

  Named numeric vector with elements:

  `r1`

  :   Growth rate of species 1. Default: `1`.

  `r2`

  :   Growth rate of species 2. Default: `1`.

  `K1`

  :   Carrying capacity of species 1. Default: `10`.

  `K2`

  :   Carrying capacity of species 2. Default: `10`.

  `a12`

  :   Effect of species 2 on species 1 (\\\alpha\_{12}\\). Default:
      `0.5`.

  `a21`

  :   Effect of species 1 on species 2 (\\\alpha\_{21}\\). Default:
      `0.5`.

## Value

A list with one element: a numeric vector of length 2 containing
\\(dN_1/dt, dN_2/dt)\\, as required by
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html).

## Details

\$\$\frac{dN_1}{dt} = r_1 N_1 \left(1 - \frac{N_1 + \alpha\_{12}
N_2}{K_1}\right)\$\$ \$\$\frac{dN_2}{dt} = r_2 N_2 \left(1 - \frac{N_2 +
\alpha\_{21} N_1}{K_2}\right)\$\$

where \\N_1, N_2\\ are species abundances, \\r_1, r_2\\ are intrinsic
growth rates, \\K_1, K_2\\ are carrying capacities, and \\\alpha\_{12},
\alpha\_{21}\\ are interspecific competition coefficients (the effect of
species 2 on species 1, and vice versa).

The four possible outcomes (competitive exclusion of species 1 or 2,
stable coexistence, or an unstable equilibrium with priority effects)
depend on the relative magnitudes of \\K_1\\, \\K_2\\, \\\alpha\_{12}\\,
and \\\alpha\_{21}\\, making this an excellent teaching example for
phase plane analysis.

## See also

[`ode_lotka_volterra()`](https://jmgraham30.github.io/ggphasr/reference/ode_lotka_volterra.md)

## Examples

``` r
ode_competition(t = 0, y = c(5, 5),
               parameters = c(r1 = 1, r2 = 1,
                              K1 = 10, K2 = 10,
                              a12 = 0.5, a21 = 0.5))
#> [[1]]
#> [1] 1.25 1.25
#> 

if (FALSE) { # \dontrun{
gg_flow_field(ode_competition, xlim = c(0, 15), ylim = c(0, 15))
} # }
```
