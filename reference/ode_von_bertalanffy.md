# Von Bertalanffy growth model

The ODE for von Bertalanffy growth, widely used in fisheries biology and
ecology to model growth of organisms:

## Usage

``` r
ode_von_bertalanffy(t, y, parameters = c(alpha = 1, beta = 0.5))
```

## Arguments

- t:

  Numeric scalar. Time (not used directly, included for deSolve
  compatibility).

- y:

  Numeric vector of length 1. Current state: \\y\[1\]\\ is body mass
  (must be non-negative).

- parameters:

  Named numeric vector with elements:

  `alpha`

  :   Anabolism coefficient. Default: `1`.

  `beta`

  :   Catabolism coefficient. Default: `0.5`.

## Value

A list with one element: a numeric vector of length 1 containing
\\dy/dt\\, as required by
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html).

## Details

\$\$\frac{dy}{dt} = \alpha \\ y^{2/3} - \beta \\ y\$\$

where \\y(t)\\ is body mass (or a proportional measure of body size),
\\\alpha \> 0\\ is the anabolism coefficient (growth rate proportional
to surface area, hence the \\2/3\\ power), and \\\beta \> 0\\ is the
catabolism coefficient (decay proportional to mass).

The model has a single non-trivial stable equilibrium at \\y^\* =
(\alpha / \beta)^3\\, which corresponds to the asymptotic body size
\\W\_\infty\\ in fisheries notation.

Note that \\y^{2/3}\\ is not defined for \\y \< 0\\. Although the
biological interpretation restricts \\y \geq 0\\, no guard is applied
internally — users should restrict phase portrait axes accordingly.

## See also

[`ode_monomolecular()`](https://jmgraham30.github.io/ggphasr/reference/ode_monomolecular.md),
[`ode_logistic()`](https://jmgraham30.github.io/ggphasr/reference/ode_logistic.md)

## Examples

``` r
# Evaluate at y = 8 with default parameters (alpha = 1, beta = 0.5)
# dy/dt = 1 * 8^(2/3) - 0.5 * 8 = 4 - 4 = 0  (equilibrium)
ode_von_bertalanffy(t = 0, y = c(8), parameters = c(alpha = 1, beta = 0.5))
#> [[1]]
#> [1] -4.440892e-16
#> 

if (FALSE) { # \dontrun{
gg_phase_portrait(ode_von_bertalanffy, ylim = c(0, 12))
} # }
```
