# Lindemann mechanism (chemical kinetics)

The Lindemann-Christiansen mechanism for unimolecular gas-phase
reactions, tracking the concentrations of reactant \\A\\ and activated
intermediate \\A^\*\\:

## Usage

``` r
ode_lindemann(t, y, parameters = c(k1 = 1, k_1 = 1, k2 = 0.5))
```

## Arguments

- t:

  Numeric scalar. Time (autonomous; included for deSolve compatibility).

- y:

  Numeric vector of length 2. State vector: \\y\[1\]\\ = \\\[A\]\\,
  \\y\[2\]\\ = \\\[A^\*\]\\.

- parameters:

  Named numeric vector with elements:

  `k1`

  :   Activation rate constant. Default: `1`.

  `k_1`

  :   Deactivation rate constant. Default: `1`.

  `k2`

  :   Product formation rate constant. Default: `0.5`.

## Value

A list with one element: a numeric vector of length 2 containing
\\(d\[A\]/dt, d\[A^\*\]/dt)\\, as required by
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html).

## Details

Reaction scheme: \$\$A + A \xrightarrow{k_1} A^\* + A \quad
\text{(activation)}\$\$ \$\$A^\* + A \xrightarrow{k\_{-1}} A + A \quad
\text{(deactivation)}\$\$ \$\$A^\* \xrightarrow{k_2} P \quad
\text{(product formation)}\$\$

The resulting ODE system for \\\[A\]\\ and \\\[A^\*\]\\ is:

\$\$\frac{d\[A\]}{dt} = -k_1 \[A\]^2 + k\_{-1} \[A^\*\]\[A\]\$\$
\$\$\frac{d\[A^\*\]}{dt} = k_1 \[A\]^2 - k\_{-1} \[A^\*\]\[A\] - k_2
\[A^\*\]\$\$

## Examples

``` r
ode_lindemann(t = 0, y = c(2, 0.1),
             parameters = c(k1 = 1, k_1 = 1, k2 = 0.5))
#> [[1]]
#> [1] -3.80  3.75
#> 

if (FALSE) { # \dontrun{
gg_flow_field(ode_lindemann, xlim = c(0, 3), ylim = c(0, 1))
} # }
```
