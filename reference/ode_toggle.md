# Genetic toggle switch model

The Gardner et al. (2000) mutual repressor model of a synthetic genetic
toggle switch:

## Usage

``` r
ode_toggle(t, y, parameters = c(alpha1 = 3, alpha2 = 3, beta = 2, gamma = 2))
```

## Arguments

- t:

  Numeric scalar. Time (autonomous; included for deSolve compatibility).

- y:

  Numeric vector of length 2. State vector: \\y\[1\]\\ = \\u\\,
  \\y\[2\]\\ = \\v\\.

- parameters:

  Named numeric vector with elements:

  `alpha1`

  :   Effective synthesis rate of \\u\\. Default: `3`.

  `alpha2`

  :   Effective synthesis rate of \\v\\. Default: `3`.

  `beta`

  :   Cooperativity coefficient in \\u\\ equation. Default: `2`.

  `gamma`

  :   Cooperativity coefficient in \\v\\ equation. Default: `2`.

## Value

A list with one element: a numeric vector of length 2 containing
\\(du/dt, dv/dt)\\, as required by
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html).

## Details

\$\$\frac{du}{dt} = \frac{\alpha_1}{1 + v^\beta} - u\$\$
\$\$\frac{dv}{dt} = \frac{\alpha_2}{1 + u^\gamma} - v\$\$

where \\u\\ and \\v\\ are the concentrations of two mutually repressing
proteins. \\\alpha_1, \alpha_2 \> 0\\ are the effective synthesis rates
(incorporating promoter strength and degradation), and \\\beta, \gamma
\> 0\\ are the cooperativity (Hill) coefficients.

For sufficiently large \\\alpha_1\\ and \\\alpha_2\\ with cooperativity
\\\> 1\\, the system is bistable: two stable equilibria (each
corresponding to one protein being dominant) separated by an unstable
saddle point. This bistability is the basis of a synthetic biological
memory device.

## References

Gardner TS, Cantor CR, Collins JJ (2000). Construction of a genetic
toggle switch in Escherichia coli. *Nature* 403: 339–342.
[doi:10.1038/35002131](https://doi.org/10.1038/35002131)

## Examples

``` r
ode_toggle(t = 0, y = c(2, 0.5),
           parameters = c(alpha1 = 3, alpha2 = 3, beta = 2, gamma = 2))
#> [[1]]
#> [1] 0.4 0.1
#> 

if (FALSE) { # \dontrun{
gg_flow_field(ode_toggle, xlim = c(0, 4), ylim = c(0, 4))
} # }
```
