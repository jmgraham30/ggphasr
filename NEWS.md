# ggphasr (development version)

## Initial development release

* Added `gg_flow_field()` for direction/velocity field visualization of 1D
  and 2D autonomous ODE systems, with support for equal-length and
  proportional arrow scaling, magnitude color mapping, and multi-system
  overlay.

* Added `gg_nullclines()` for plotting x- and y-nullclines via contour
  methods, with configurable colors, line types, and legend support.

* Added `gg_trajectory()` for numerically integrated solution paths from
  one or more initial conditions, with forward and backward integration,
  direction arrows, and automatic color scaling across trajectories.

* Added `gg_phase_portrait()` for 1D phase line visualization with
  directional arrows and filled/open circle markers at stable/unstable
  equilibria.

* Added `gg_time_series()` for time-domain plots of ODE solutions, with
  faceted panels for 2D systems and multiple initial condition support.

* Added `gg_manifolds()` for stable and unstable manifold computation and
  visualization at saddle points of 2D systems.

* Added `gg_phase_plane()` as a high-level all-in-one wrapper producing
  flow fields, nullclines, trajectories, and equilibrium annotations in a
  single call. Returns a named list with `$plot` and `$equilibria`
  components.

* Added `add_layer()` for adding ggplot2 layers to a `ggphasr_result`
  object while preserving the `$equilibria` data frame.

* Added `find_equilibrium()` for numerical equilibrium finding via
  Newton-Raphson (rootSolve), with single-guess and grid-search modes.

* Added `classify_equilibrium()` for trace-determinant stability
  classification using a finite-difference Jacobian. Returns a tidy
  data frame row.

* Added `theme_phase_plane()`, a ggplot2 theme designed for phase plane
  plots with white background, light grey grid lines, and centered axes.

* Added built-in ODE systems: four 1D growth models, eight classical 2D
  systems, and fifteen generic textbook examples (ported from phaseR).

* Both phaseR/deSolve-compatible (Convention A) and simplified (Convention B)
  ODE function signatures are supported throughout, with automatic detection.
