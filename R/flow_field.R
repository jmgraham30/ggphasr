# flow_field.R
#
# gg_flow_field() — plot a direction/velocity field for a 1D or 2D
# autonomous ODE system using ggplot2.
#
# Internal helpers defined in this file:
#   .compute_flow_field()  — evaluates the ODE on a grid; returns a data frame
#   .scale_arrows()        — scales arrow vectors for display
#   .build_flow_field_layers() — constructs the ggplot2 layer list


# ---------------------------------------------------------------------------
# .compute_flow_field()
# ---------------------------------------------------------------------------
#
# Evaluates a normalized (Convention A) ODE on a regular grid and returns
# a data frame with one row per grid point.
#
# For 1D systems (`system = "one.dim"`) the grid is over `y` only; `x` is
# set to a fixed value of 1 (used as the time axis placeholder, consistent
# with phaseR's behavior). The returned data frame has columns:
#   x, y, dx, dy
# where for 1D systems `x` is the time-axis position and `dx = 1` (unit
# rightward arrow base), `dy = f(y)`.
#
# For 2D systems the grid is over both `x` and `y`. The returned data frame
# has columns: x, y, dx, dy.
#
# @param deriv       Normalized (Convention A) ODE function.
# @param system      "one.dim" or "two.dim".
# @param xlim        Numeric vector of length 2: x-axis range.
# @param ylim        Numeric vector of length 2: y-axis range.
# @param n_points    Integer: number of grid points along each axis.
# @param parameters  Parameter vector/list passed to `deriv`.
# @return A data frame with columns x, y, dx, dy.
#
.compute_flow_field <- function(deriv,
                                 system,
                                 xlim,
                                 ylim,
                                 n_points,
                                 parameters) {

  if (system == "two.dim") {
    xs <- seq(xlim[[1L]], xlim[[2L]], length.out = n_points)
    ys <- seq(ylim[[1L]], ylim[[2L]], length.out = n_points)
    grid <- expand.grid(x = xs, y = ys)
    derivs <- mapply(
      function(xi, yi) {
        .eval_ode(deriv, t = 0, y = c(xi, yi), parameters = parameters)
      },
      grid$x, grid$y
    )
    grid$dx <- derivs[1L, ]
    grid$dy <- derivs[2L, ]

  } else {
    # 1D: x axis is the independent variable (time-like), y is the state
    xs <- seq(xlim[[1L]], xlim[[2L]], length.out = n_points)
    ys <- seq(ylim[[1L]], ylim[[2L]], length.out = n_points)
    grid <- expand.grid(x = xs, y = ys)
    grid$dx <- 1   # unit rightward arrow (time axis)
    grid$dy <- vapply(grid$y, function(yi) {
      .eval_ode(deriv, t = 0, y = c(yi), parameters = parameters)[[1L]]
    }, numeric(1L))
  }

  grid
}


# ---------------------------------------------------------------------------
# .scale_arrows()
# ---------------------------------------------------------------------------
#
# Scales (dx, dy) vectors for display as arrows on the grid.
#
# Two modes:
#   "equal"        — normalizes all arrows to the same display length
#   "proportional" — scales by magnitude, bounded by `max_magnitude`
#
# The maximum arrow display length is `cell_size * arrow_length_scale`,
# where `cell_size` is the spacing between adjacent grid points (so arrows
# can never overflow into adjacent cells).
#
# @param df              Data frame with columns dx, dy.
# @param arrow_type      "equal" or "proportional".
# @param xlim, ylim      Axis ranges (used to compute cell_size).
# @param n_points        Grid resolution.
# @param arrow_length_scale Fraction of cell size for max arrow length.
#   Default 0.85 leaves a small gap between arrows.
# @param max_magnitude   For "proportional" only: if supplied, magnitudes are
#   scaled relative to this value rather than the grid maximum.
# @return Data frame with additional columns xend, yend (arrow endpoints).
#
.scale_arrows <- function(df,
                           arrow_type,
                           xlim,
                           ylim,
                           n_points,
                           arrow_length_scale = 0.85,
                           max_magnitude      = NULL) {

  # Cell size: the spacing between adjacent grid points
  cell_x <- diff(xlim) / (n_points - 1)
  cell_y <- diff(ylim) / (n_points - 1)
  cell_size <- min(cell_x, cell_y)
  max_len   <- cell_size * arrow_length_scale

  magnitudes <- sqrt(df$dx^2 + df$dy^2)

  if (arrow_type == "equal") {
    # Normalize to unit length then scale to max_len
    # Avoid dividing by zero at equilibrium points
    safe_mag <- ifelse(magnitudes < .Machine$double.eps, 1, magnitudes)
    df$dx_scaled <- df$dx / safe_mag * max_len
    df$dy_scaled <- df$dy / safe_mag * max_len

  } else {
    # Proportional: scale relative to maximum magnitude on the grid
    # (or user-supplied max_magnitude)
    ref_mag <- if (!is.null(max_magnitude)) {
      max_magnitude
    } else {
      max(magnitudes, na.rm = TRUE)
    }
    # Avoid division by zero if all derivatives are zero
    if (ref_mag < .Machine$double.eps) ref_mag <- 1
    df$dx_scaled <- df$dx / ref_mag * max_len
    df$dy_scaled <- df$dy / ref_mag * max_len
  }

  df$xend <- df$x + df$dx_scaled
  df$yend <- df$y + df$dy_scaled

  # Store magnitude for optional color mapping
  df$magnitude <- magnitudes

  df
}


# ---------------------------------------------------------------------------
# gg_flow_field()
# ---------------------------------------------------------------------------

#' Plot a flow field (direction field) for an ODE system
#'
#' Computes and plots the direction or velocity field of a one- or
#' two-dimensional autonomous ODE system on a regular grid. Returns a
#' complete [ggplot2::ggplot()] object to which additional layers (nullclines,
#' trajectories, etc.) can be added with `+`.
#'
#' Multiple ODE systems can be overlaid on a single plot by passing a named
#' list to the `deriv` argument (see Details and Examples).
#'
#' @param deriv A function (or named list of functions) describing the ODE
#'   system, in either Convention A (`f(t, y, parameters)` returning
#'   `list(c(...))`) or Convention B (simplified `f(x, y, parameters)`
#'   returning `c(...)`). See [ggphasr-package] for details on ODE conventions.
#'
#'   To overlay multiple systems, pass a named list:
#'   `deriv = list(system1 = f1, system2 = f2)`.
#'
#' @param xlim Numeric vector of length 2. Range of the x-axis (or the
#'   time axis for 1D systems). Required.
#' @param ylim Numeric vector of length 2. Range of the y-axis (state
#'   variable axis). Required.
#' @param system Character string: `"two.dim"` (default) or `"one.dim"`.
#' @param parameters A numeric vector or list of parameter values passed to
#'   `deriv`. When `deriv` is a list of functions, `parameters` can be a
#'   named list of parameter vectors, one per system; or a single vector
#'   applied to all systems.
#' @param n_points Integer. Number of grid points along each axis.
#'   Default: `21` (matching phaseR's default).
#' @param arrow_type Character string: `"equal"` (default, all arrows the
#'   same length, showing direction only) or `"proportional"` (arrow length
#'   proportional to vector magnitude).
#' @param arrow_color Character. Color of the arrows. Default: `"grey60"`.
#'   Ignored when `color_by_magnitude = TRUE`.
#' @param arrow_size Numeric. Relative size of the arrow heads.
#'   Default: `0.25`.
#' @param arrow_linewidth Numeric. Line width of the arrow shafts.
#'   Default: `0.4`.
#' @param arrow_length_scale Numeric in (0, 1]. Maximum arrow length as a
#'   fraction of the grid cell size. Default: `0.85`.
#' @param max_magnitude Numeric or `NULL`. For `arrow_type = "proportional"`
#'   only: if supplied, arrow lengths are scaled relative to this reference
#'   magnitude rather than the maximum on the grid. Useful for producing
#'   consistent scaling across multiple plots. Default: `NULL` (auto-scale
#'   to grid maximum).
#' @param color_by_magnitude Logical. If `TRUE`, arrows are colored by
#'   vector magnitude using a continuous color scale. Overrides
#'   `arrow_color`. Default: `FALSE`.
#' @param magnitude_palette Character vector of length 2. Low and high colors
#'   for the magnitude color scale, used when `color_by_magnitude = TRUE`.
#'   Default: `c("grey80", "#2c7bb6")`.
#' @param add_origin_lines Logical. If `TRUE` (default), adds thin reference
#'   lines at `x = 0` and `y = 0` (or `y = 0` only for 1D systems) via
#'   [ggplot2::geom_hline()] and [ggplot2::geom_vline()].
#' @param origin_line_color Character. Color of the origin reference lines.
#'   Default: `"grey40"`.
#' @param xlab Character. x-axis label. Default: `"x"` for 2D systems,
#'   `"t"` for 1D systems.
#' @param ylab Character. y-axis label. Default: `"y"`.
#' @param title Character or `NULL`. Plot title. Default: `NULL`.
#'
#' @return A [ggplot2::ggplot()] object. Additional layers, scales, and theme
#'   elements can be added with `+`.
#'
#' @details
#' ## ODE conventions
#' Both phaseR-style (Convention A) and simplified (Convention B) ODE
#' functions are accepted. The calling convention is detected automatically
#' from the function's argument names.
#'
#' ## Multiple systems
#' When `deriv` is a named list, each system is drawn with a different color
#' (taken from a discrete palette) and a legend is added automatically.
#' `parameters` should then be a named list of the same length as `deriv`,
#' with one parameter vector per system. If `parameters` is a single vector,
#' it is applied to all systems.
#'
#' ## 1D systems
#' For `system = "one.dim"`, the x-axis represents time (the independent
#' variable) and the y-axis represents the state variable. Arrows point
#' rightward (increasing time) and up or down according to `dy/dt`.
#'
#' @examples
#' # 2D system: Lotka-Volterra phase plane
#' gg_flow_field(
#'   ode_lotka_volterra,
#'   xlim       = c(0, 5),
#'   ylim       = c(0, 5),
#'   parameters = c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1)
#' )
#'
#' # 1D system: logistic growth phase portrait
#' gg_flow_field(
#'   ode_logistic,
#'   xlim       = c(0, 4),
#'   ylim       = c(-1, 12),
#'   system     = "one.dim",
#'   parameters = c(r = 1, K = 10)
#' )
#'
#' # Proportional arrows colored by magnitude
#' gg_flow_field(
#'   ode_van_der_pol,
#'   xlim                = c(-3, 3),
#'   ylim                = c(-3, 3),
#'   parameters          = c(mu = 1),
#'   arrow_type          = "proportional",
#'   color_by_magnitude  = TRUE
#' )
#'
#' # Compose with nullclines and a custom title
#' gg_flow_field(
#'   ode_lotka_volterra,
#'   xlim       = c(0, 5),
#'   ylim       = c(0, 5),
#'   parameters = c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1),
#'   title      = "Lotka-Volterra"
#' ) +
#'   ggplot2::labs(x = "Prey", y = "Predator")
#'
#' @seealso [ggphasr::gg_nullclines()], [ggphasr::gg_trajectory()],
#'   [ggphasr::gg_phase_portrait()], [ggphasr::theme_phase_plane()]
#' @export
gg_flow_field <- function(deriv,
                           xlim,
                           ylim,
                           system              = c("two.dim", "one.dim"),
                           parameters          = NULL,
                           n_points            = 21L,
                           arrow_type          = c("equal", "proportional"),
                           arrow_color         = "grey60",
                           arrow_size          = 0.25,
                           arrow_linewidth     = 0.4,
                           arrow_length_scale  = 0.85,
                           max_magnitude       = NULL,
                           color_by_magnitude  = FALSE,
                           magnitude_palette   = c("grey80", "#2c7bb6"),
                           add_origin_lines    = TRUE,
                           origin_line_color   = "grey40",
                           xlab                = NULL,
                           ylab                = "y",
                           title               = NULL) {

  # ── Input validation ────────────────────────────────────────────────────
  system     <- match.arg(system)
  arrow_type <- match.arg(arrow_type)

  if (!is.numeric(xlim) || length(xlim) != 2L || xlim[[1L]] >= xlim[[2L]]) {
    rlang::abort("`xlim` must be a numeric vector of length 2 with xlim[1] < xlim[2].")
  }
  if (!is.numeric(ylim) || length(ylim) != 2L || ylim[[1L]] >= ylim[[2L]]) {
    rlang::abort("`ylim` must be a numeric vector of length 2 with ylim[1] < ylim[2].")
  }
  if (!is.numeric(n_points) || n_points < 2L) {
    rlang::abort("`n_points` must be an integer >= 2.")
  }
  n_points <- as.integer(n_points)

  # Default x-axis label
  if (is.null(xlab)) xlab <- if (system == "one.dim") "t" else "x"

  # ── Handle multiple systems ─────────────────────────────────────────────
  multi_system <- is.list(deriv) && !is.function(deriv)

  if (multi_system) {
    # Validate: must be a named list of functions
    if (is.null(names(deriv)) || any(names(deriv) == "")) {
      rlang::abort(
        "When `deriv` is a list, all elements must be named ",
        "(e.g. `list(system1 = f1, system2 = f2)`)."
      )
    }
    if (!all(vapply(deriv, is.function, logical(1L)))) {
      rlang::abort("All elements of `deriv` must be functions.")
    }

    # Broadcast parameters: if a single vector, replicate for each system
    if (!is.list(parameters) || is.null(names(parameters))) {
      params_list <- stats::setNames(
        replicate(length(deriv), parameters, simplify = FALSE),
        names(deriv)
      )
    } else {
      params_list <- parameters
    }

    # Build data for each system
    system_names <- names(deriv)
    all_data <- lapply(system_names, function(nm) {
      fn     <- deriv[[nm]]
      params <- params_list[[nm]]
      norm   <- .normalize_ode(fn, system = system)
      .validate_ode(norm, system = system, parameters = params)
      df     <- .compute_flow_field(norm, system, xlim, ylim,
                                     n_points, params)
      df     <- .scale_arrows(df, arrow_type, xlim, ylim, n_points,
                               arrow_length_scale, max_magnitude)
      df$system_name <- nm
      df
    })
    field_data <- do.call(rbind, all_data)

  } else {
    # Single system
    norm <- .normalize_ode(deriv, system = system)
    .validate_ode(norm, system = system, parameters = parameters)
    field_data <- .compute_flow_field(norm, system, xlim, ylim,
                                       n_points, parameters)
    field_data <- .scale_arrows(field_data, arrow_type, xlim, ylim,
                                 n_points, arrow_length_scale, max_magnitude)
    field_data$system_name <- "single"
  }

  # ── Build the ggplot ────────────────────────────────────────────────────
  arrow_spec <- grid::arrow(length = grid::unit(arrow_size, "lines"),
                             type   = "closed")

  p <- ggplot2::ggplot()

  # ── Arrow layers ─────────────────────────────────────────────────────────
  if (multi_system) {
    # Color by system name (discrete)
    p <- p +
      ggplot2::geom_segment(
        data    = field_data,
        mapping = ggplot2::aes(x      = .data$x,
                               y      = .data$y,
                               xend   = .data$xend,
                               yend   = .data$yend,
                               color  = .data$system_name),
        arrow     = arrow_spec,
        linewidth = arrow_linewidth,
        na.rm     = TRUE
      ) +
      ggplot2::scale_color_brewer(palette = "Set1", name = NULL)

  } else if (color_by_magnitude) {
    # Color by vector magnitude (continuous)
    p <- p +
      ggplot2::geom_segment(
        data    = field_data,
        mapping = ggplot2::aes(x      = .data$x,
                               y      = .data$y,
                               xend   = .data$xend,
                               yend   = .data$yend,
                               color  = .data$magnitude),
        arrow     = arrow_spec,
        linewidth = arrow_linewidth,
        na.rm     = TRUE
      ) +
      ggplot2::scale_color_gradient(
        low  = magnitude_palette[[1L]],
        high = magnitude_palette[[2L]],
        name = "Magnitude"
      )

  } else {
    # Single color
    p <- p +
      ggplot2::geom_segment(
        data    = field_data,
        mapping = ggplot2::aes(x    = .data$x,
                               y    = .data$y,
                               xend = .data$xend,
                               yend = .data$yend),
        color     = arrow_color,
        arrow     = arrow_spec,
        linewidth = arrow_linewidth,
        na.rm     = TRUE
      )
  }

  # ── Origin reference lines ───────────────────────────────────────────────
  if (add_origin_lines) {
    p <- p +
      ggplot2::geom_hline(yintercept = 0,
                          color      = origin_line_color,
                          linewidth  = 0.5,
                          linetype   = "solid")
    if (system == "two.dim") {
      p <- p +
        ggplot2::geom_vline(xintercept = 0,
                            color      = origin_line_color,
                            linewidth  = 0.5,
                            linetype   = "solid")
    }
  }

  # ── Coordinate system and scales ─────────────────────────────────────────
  p <- p +
    ggplot2::coord_cartesian(xlim   = xlim,
                              ylim   = ylim,
                              expand = FALSE) +
    ggplot2::labs(x = xlab, y = ylab)

  # ── Title ────────────────────────────────────────────────────────────────
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }

  # ── Theme ────────────────────────────────────────────────────────────────
  p <- p + theme_phase_plane()

  p
}
