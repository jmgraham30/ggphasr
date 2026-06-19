# SIR epidemic model

The classic Kermack-McKendrick SIR (Susceptible-Infected-Recovered)
epidemic model. Supports two formulations via the `scale` parameter:

## Usage

``` r
ode_sir(t, y, parameters = c(beta = 0.5, gamma = 0.1))
```

## Arguments

- t:

  Numeric scalar. Time (autonomous system; included for deSolve
  compatibility).

- y:

  Numeric vector of length 2. State vector: \\y\[1\]\\ = \\S\\
  (susceptible), \\y\[2\]\\ = \\I\\ (infected).

- parameters:

  Named numeric vector with elements:

  `beta`

  :   Transmission rate. Default: `0.5`.

  `gamma`

  :   Recovery rate. Default: `0.1`.

  `N`

  :   Total population size (used only when `scale = "counts"`).
      Default: `1000`.

  `scale`

  :   Character: `"proportions"` (default) or `"counts"`. Can be passed
      as a named element of the `parameters` list, or as a separate
      argument.

## Value

A list with one element: a numeric vector of length 2 containing
\\(dS/dt, dI/dt)\\, as required by
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html).

## Details

**Proportions** (`scale = "proportions"`, default): \$\$\frac{dS}{dt} =
-\beta S I\$\$ \$\$\frac{dI}{dt} = \beta S I - \gamma I\$\$

where \\S + I + R = 1\\ at all times. Here \\\beta\\ is the transmission
rate and \\\gamma\\ is the recovery rate. The basic reproduction number
is \\R_0 = \beta / \gamma\\.

**Counts** (`scale = "counts"`): \$\$\frac{dS}{dt} = -\beta S I / N\$\$
\$\$\frac{dI}{dt} = \beta S I / N - \gamma I\$\$

where \\N\\ is the total (constant) population size. The proportions
formulation is recovered by setting \\N = 1\\.

Note: only \\S\\ and \\I\\ are tracked as state variables. \\R\\ can be
recovered as \\R = N - S - I\\.

The `scale` argument is passed as part of `parameters` to maintain
deSolve compatibility. Because mixing character and numeric values in
[`c()`](https://rdrr.io/r/base/c.html) coerces everything to character,
pass `parameters` as a [`list()`](https://rdrr.io/r/base/list.html) when
including `scale`:
`parameters = list(beta = 0.5, gamma = 0.1, N = 1000, scale = "counts")`.
Passing a plain numeric [`c()`](https://rdrr.io/r/base/c.html) vector
(without `scale`) also works and is the recommended form for the
proportions formulation. When `scale` is absent, `"proportions"` is
assumed.

## Examples

``` r
# Proportions formulation — plain numeric vector is fine
ode_sir(t = 0, y = c(0.99, 0.01),
        parameters = c(beta = 0.5, gamma = 0.1))
#> [[1]]
#> [1] -0.00495  0.00395
#> 

# Counts formulation — use list() to avoid character coercion
ode_sir(t = 0, y = c(990, 10),
        parameters = list(beta = 0.5, gamma = 0.1, N = 1000,
                          scale = "counts"))
#> [[1]]
#> [1] -4.95  3.95
#> 

if (FALSE) { # \dontrun{
gg_flow_field(ode_sir, xlim = c(0, 1), ylim = c(0, 1))
} # }
```
