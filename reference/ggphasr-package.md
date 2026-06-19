# ggphasr: Phase Plane Analysis of ODE Systems Using ggplot2

`ggphasr` provides tools for qualitative analysis of one- and
two-dimensional autonomous ordinary differential equation (ODE) systems
using phase plane methods. All visualizations are produced with ggplot2
and its extensions, making plots fully customizable with standard
ggplot2 syntax.

## ODE conventions

All functions accept ODE systems in one of two calling conventions:

**Convention A** (deSolve-compatible, primary):

    f <- function(t, y, parameters) {
      list(c(dy1, dy2))
    }

**Convention B** (simplified):

    # 2D:
    f <- function(x, y, parameters = NULL) c(dx, dy)
    # 1D:
    f <- function(y, parameters = NULL) dy

Convention B functions are automatically detected and wrapped into
Convention A internally. All built-in ODE systems use Convention A.

## Main functions

**Plotting (return ggplot objects or layer lists):**

- [`gg_flow_field()`](https://jmgraham30.github.io/ggphasr/reference/gg_flow_field.md):

  Direction/velocity field arrows on a grid

- [`gg_nullclines()`](https://jmgraham30.github.io/ggphasr/reference/gg_nullclines.md):

  Zero-isocline curves (x- and y-nullclines)

- [`gg_trajectory()`](https://jmgraham30.github.io/ggphasr/reference/gg_trajectory.md):

  Numerically integrated solution paths

- [`gg_phase_portrait()`](https://jmgraham30.github.io/ggphasr/reference/gg_phase_portrait.md):

  1D phase line with directional arrows

- [`gg_time_series()`](https://jmgraham30.github.io/ggphasr/reference/gg_time_series.md):

  Time-series plots of state variables

- [`gg_manifolds()`](https://jmgraham30.github.io/ggphasr/reference/gg_manifolds.md):

  Stable/unstable manifolds of saddle points

- [`gg_phase_plane()`](https://jmgraham30.github.io/ggphasr/reference/gg_phase_plane.md):

  All-in-one phase plane analysis wrapper

**Analysis (return data structures, no plots):**

- [`find_equilibrium()`](https://jmgraham30.github.io/ggphasr/reference/find_equilibrium.md):

  Numerical root-finding for equilibria

- [`classify_equilibrium()`](https://jmgraham30.github.io/ggphasr/reference/classify_equilibrium.md):

  Trace-determinant stability classification

**Built-in ODE systems:**

- 1D models:

  [`ode_exponential()`](https://jmgraham30.github.io/ggphasr/reference/ode_exponential.md),
  [`ode_logistic()`](https://jmgraham30.github.io/ggphasr/reference/ode_logistic.md),
  [`ode_monomolecular()`](https://jmgraham30.github.io/ggphasr/reference/ode_monomolecular.md),
  [`ode_von_bertalanffy()`](https://jmgraham30.github.io/ggphasr/reference/ode_von_bertalanffy.md)

- 2D models:

  [`ode_lotka_volterra()`](https://jmgraham30.github.io/ggphasr/reference/ode_lotka_volterra.md),
  [`ode_sir()`](https://jmgraham30.github.io/ggphasr/reference/ode_sir.md),
  [`ode_van_der_pol()`](https://jmgraham30.github.io/ggphasr/reference/ode_van_der_pol.md),
  [`ode_simple_pendulum()`](https://jmgraham30.github.io/ggphasr/reference/ode_simple_pendulum.md),
  [`ode_competition()`](https://jmgraham30.github.io/ggphasr/reference/ode_competition.md),
  [`ode_toggle()`](https://jmgraham30.github.io/ggphasr/reference/ode_toggle.md),
  [`ode_morris_lecar()`](https://jmgraham30.github.io/ggphasr/reference/ode_morris_lecar.md),
  [`ode_lindemann()`](https://jmgraham30.github.io/ggphasr/reference/ode_lindemann.md)

- Textbook examples:

  [`ode_example_01()`](https://jmgraham30.github.io/ggphasr/reference/ode_example_01.md)
  through
  [`ode_example_15()`](https://jmgraham30.github.io/ggphasr/reference/ode_example_15.md)

## Typical workflow

    library(ggphasr)
    library(ggplot2)

    lv_params <- c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1)

    # Quick all-in-one analysis
    result <- gg_phase_plane(ode_lotka_volterra,
                              xlim = c(0, 5), ylim = c(0, 5),
                              parameters = lv_params)
    result$plot
    result$equilibria

    # Composable layer-by-layer approach
    gg_flow_field(ode_lotka_volterra,
                  xlim = c(0, 5), ylim = c(0, 5),
                  parameters = lv_params) +
      gg_nullclines(ode_lotka_volterra,
                    xlim = c(0, 5), ylim = c(0, 5),
                    parameters = lv_params) +
      gg_trajectory(ode_lotka_volterra,
                    y0 = c(1, 1), xlim = c(0, 5), ylim = c(0, 5),
                    parameters = lv_params)

## References

Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
[doi:10.32614/RJ-2014-023](https://doi.org/10.32614/RJ-2014-023)

## See also

Useful links:

- <https://github.com/jmgraham30/ggphasr>

- Report bugs at <https://github.com/jmgraham30/ggphasr/issues>

## Author

**Maintainer**: Jason Graham <prof.jason.m.graham@gmail.com>
([ORCID](https://orcid.org/0000-0003-0047-7178))

Authors:

- Jason Graham <prof.jason.m.graham@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-0047-7178))
