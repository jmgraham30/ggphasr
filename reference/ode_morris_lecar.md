# Morris-Lecar neuron model

The Morris-Lecar (1981) conductance-based neuron model, a 2D reduction
of Hodgkin-Huxley dynamics for a barnacle muscle fiber:

## Usage

``` r
ode_morris_lecar(
  t,
  y,
  parameters = c(I = 0, C = 20, gCa = 4.4, gK = 8, gL = 2, VCa = 120, VK = -84, VL = -60,
    V1 = -1.2, V2 = 18, V3 = 2, V4 = 30, phi = 0.04)
)
```

## Arguments

- t:

  Numeric scalar. Time (autonomous; included for deSolve compatibility).

- y:

  Numeric vector of length 2. State vector: \\y\[1\]\\ = membrane
  potential \\V\\ (mV), \\y\[2\]\\ = K\\^+\\ channel open probability
  \\N\\ (dimensionless, \\0 \leq N \leq 1\\).

- parameters:

  Named numeric vector with elements:

  `I`

  :   Applied current (\\\mu\\A/cm\\^2\\). Default: `0.`

  `C`

  :   Membrane capacitance (\\\mu\\F/cm\\^2\\). Default: `20`.

  `gCa`

  :   Maximum Ca\\^{2+}\\ conductance (mS/cm\\^2\\). Default: `4.4`.

  `gK`

  :   Maximum K\\^+\\ conductance (mS/cm\\^2\\). Default: `8`.

  `gL`

  :   Leak conductance (mS/cm\\^2\\). Default: `2`.

  `VCa`

  :   Ca\\^{2+}\\ reversal potential (mV). Default: `120`.

  `VK`

  :   K\\^+\\ reversal potential (mV). Default: `-84`.

  `VL`

  :   Leak reversal potential (mV). Default: `-60`.

  `V1`

  :   Voltage at half-activation of Ca\\^{2+}\\ (mV). Default: `-1.2`.

  `V2`

  :   Slope of Ca\\^{2+}\\ activation (mV). Default: `18`.

  `V3`

  :   Voltage at half-activation of K\\^+\\ (mV). Default: `2`.

  `V4`

  :   Slope of K\\^+\\ activation (mV). Default: `30`.

  `phi`

  :   Reference frequency (dimensionless). Default: `0.04`.

## Value

A list with one element: a numeric vector of length 2 containing
\\(dV/dt, dN/dt)\\, as required by
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html).

## Details

\$\$C \frac{dV}{dt} = I - g\_{Ca} M\_\infty(V)(V - V\_{Ca}) - g_K N (V -
V_K) - g_L (V - V_L)\$\$ \$\$\frac{dN}{dt} = \phi \frac{N\_\infty(V) -
N}{\tau_N(V)}\$\$

where the voltage-dependent steady-state functions are: \$\$M\_\infty(V)
= \frac{1}{2}\left(1 + \tanh\\\left(\frac{V -
V_1}{V_2}\right)\right)\$\$ \$\$N\_\infty(V) = \frac{1}{2}\left(1 +
\tanh\\\left(\frac{V - V_3}{V_4}\right)\right)\$\$ \$\$\tau_N(V) =
\left(\cosh\\\left(\frac{V - V_3}{2 V_4}\right)\right)^{-1}\$\$

The state variables are membrane potential \\V\\ (mV) and \\N\\, the
probability that a K\\^+\\ channel is open. Depending on the applied
current \\I\\ and other parameters, the system can exhibit a stable
resting state, a stable limit cycle (repetitive firing), or bistability,
making it a rich teaching example for bifurcation analysis.

## References

Morris C, Lecar H (1981). Voltage oscillations in the barnacle giant
muscle fiber. *Biophysical Journal* 35(1): 193-213.
[doi:10.1016/S0006-3495(81)84782-0](https://doi.org/10.1016/S0006-3495%2881%2984782-0)

## Examples

``` r
# At resting potential with zero applied current
ode_morris_lecar(t = 0, y = c(-60, 0),
                parameters = c(I = 0, C = 20,
                               gCa = 4.4, gK = 8, gL = 2,
                               VCa = 120, VK = -84, VL = -60,
                               V1 = -1.2, V2 = 18, V3 = 2, V4 = 30,
                               phi = 0.04))
#> [[1]]
#> [1] 0.057500749 0.000999041
#> 

if (FALSE) { # \dontrun{
gg_flow_field(ode_morris_lecar,
              xlim = c(-80, 60), ylim = c(0, 0.6))
} # }
```
