# nullclines.R
#
# gg_nullclines() — add nullcline curves to a phase plane plot.
#
# For a 2D system dx/dt = f(x,y), dy/dt = g(x,y):
#   x-nullcline: the set of points where dx/dt = f(x,y) = 0
#   y-nullcline: the set of points where dy/dt = g(x,y) = 0
#
# For a 1D system dy/dt = f(y):
#   The nullcline is simply the set of y values where f(y) = 0
#   (i.e., the equilibrium points), drawn as vertical lines on the
#   phase portrait.
#
# Nullclines are computed by evaluating the ODE on a fine grid and
# drawing contour lines at zero using ggplot2::geom_contour().
#
# Internal helpers:
#   .compute_nullcline_grid()  — evaluates f and g on a grid
#   .nullcline_layers()        — builds the ggplot2 layer list


# ---------------------------------------------------------------------------
# .compute_nullcline_grid()
# ---------------------------------------------------------------------------
#
# Evaluates the ODE derivatives on a fine rectangular grid and returns
# a data frame suitable for geom_contour().
#
# For 2D systems, returns columns: x, y, f (= dx/dt), g (= dy/dt).
# For 1D systems, returns columns: y, f (= dy/dt) only.
#
# @param deriv       Normalized (Convention A) ODE function.
# @param system      "two.dim" or "one.dim".
# @param xlim        Numeric vector of length 2: x-axis range.
# @param ylim        Numeric vector of length 2: y-axis range.
# @param n_points    Integer: grid resolution (n_points x n_points for 2D).
# @param parameters  Parameter vector/list passed to `deriv`.
# @return A data frame.
#
.compute_nullcline_grid <- function(deriv,
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
    grid$f <- derivs[1L, ]   # dx/dt
    grid$g <- derivs[2L, ]   # dy/dt

  } else {
    # 1D: only y matters
    ys      <- seq(ylim[[1L]], ylim[[2L]], length.out = n_points)
    grid    <- data.frame(y = ys)
    grid$f  <- vapply(ys, function(yi) {
      .eval_ode(deriv, t = 0, y = c(yi), parameters = parameters)[[1L]]
    }, numeric(1L))
  }

  grid
}


# ---------------------------------------------------------------------------
# gg_nullclines()
# ---------------------------------------------------------------------------

#' Add nullclines to a phase plane plot
#'
#' Computes and adds the nullclines of a one- or two-dimensional autonomous
#' ODE system to an existing phase plane plot. Returns a list of
#' [ggplot2::layer()] objects that can be added to a [ggplot2::ggplot()]
#' object with `+`.
#'
#' For a 2D system \eqn{dx/dt = f(x,y)}, \eqn{dy/dt = g(x,y)}:
#' - **x-nullcline**: the curve(s) where \eqn{f(x,y) = 0}
#' - **y-nullcline**: the curve(s) where \eqn{g(x,y) = 0}
#'
#' For a 1D system \eqn{dy/dt = f(y)}, the nullclines are the equilibrium
#' points (where \eqn{f(y) = 0}), drawn as horizontal lines on the phase
#' portrait.
#'
#' @param deriv A function describing the ODE system, in Convention A or B.
#'   See [ggphasr-package] for details.
#' @param xlim Numeric vector of length 2. x-axis range. Should match the
#'   `xlim` passed to [ggphasr::gg_flow_field()].
#' @param ylim Numeric vector of length 2. y-axis range. Should match the
#'   `ylim` passed to [ggphasr::gg_flow_field()].
#' @param system Character: `"two.dim"` (default) or `"one.dim"`.
#' @param parameters Parameter vector or list passed to `deriv`.
#' @param n_points Integer. Grid resolution for nullcline computation.
#'   Default: `250`. Higher values give smoother curves.
#' @param x_color Character. Color of the x-nullcline (where \eqn{dx/dt = 0}).
#'   Default: `"#d73027"` (red).
#' @param y_color Character. Color of the y-nullcline (where \eqn{dy/dt = 0}).
#'   Default: `"#4575b4"` (blue).
#' @param x_linetype Character or integer. Line type for the x-nullcline.
#'   Default: `"solid"`.
#' @param y_linetype Character or integer. Line type for the y-nullcline.
#'   Default: `"dashed"`.
#' @param linewidth Numeric. Line width for both nullclines. Default: `0.75`.
#' @param add_legend Logical. If `TRUE`, adds a legend entry for each
#'   nullcline. Default: `TRUE`.
#' @param legend_position Character string or numeric vector of length 2.
#'   Controls legend placement. One of `"right"` (default), `"left"`,
#'   `"top"`, `"bottom"`, `"none"`, `"inside"` (top-right corner of the
#'   panel, compact styling), or a numeric vector `c(x, y)` with values
#'   in `[0, 1]` for a custom inside position.
#'
#' @return A list of [ggplot2] layer objects. Add to a ggplot with `+`.
#'
#' @examples
#' # Standard workflow: flow field + nullclines
#' gg_flow_field(
#'   ode_lotka_volterra,
#'   xlim       = c(0, 5),
#'   ylim       = c(0, 5),
#'   parameters = c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1)
#' ) +
#'   gg_nullclines(
#'     ode_lotka_volterra,
#'     xlim       = c(0, 5),
#'     ylim       = c(0, 5),
#'     parameters = c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1)
#'   )
#'
#' # Customize nullcline appearance
#' gg_flow_field(
#'   ode_competition,
#'   xlim       = c(0, 15),
#'   ylim       = c(0, 15),
#'   parameters = c(r1=1, r2=1, K1=10, K2=10, a12=0.5, a21=0.5)
#' ) +
#'   gg_nullclines(
#'     ode_competition,
#'     xlim       = c(0, 15),
#'     ylim       = c(0, 15),
#'     parameters = c(r1=1, r2=1, K1=10, K2=10, a12=0.5, a21=0.5),
#'     x_color    = "forestgreen",
#'     y_color    = "darkorange",
#'     linewidth  = 1
#'   )
#'
#' @seealso [ggphasr::gg_flow_field()], [ggphasr::gg_trajectory()]
#' @export
gg_nullclines <- function(deriv,
                           xlim,
                           ylim,
                           system      = c("two.dim", "one.dim"),
                           parameters  = NULL,
                           n_points    = 250L,
                           x_color     = "#d73027",
                           y_color     = "#4575b4",
                           x_linetype  = "solid",
                           y_linetype  = "dashed",
                           linewidth   = 0.75,
                           add_legend       = TRUE,
                           legend_position  = "right") {

  # ── Input validation ─────────────────────────────────────────────────────
  system   <- match.arg(system)
  n_points <- as.integer(n_points)

  if (!is.numeric(xlim) || length(xlim) != 2L || xlim[[1L]] >= xlim[[2L]]) {
    rlang::abort("`xlim` must be a numeric vector of length 2 with xlim[1] < xlim[2].")
  }
  if (!is.numeric(ylim) || length(ylim) != 2L || ylim[[1L]] >= ylim[[2L]]) {
    rlang::abort("`ylim` must be a numeric vector of length 2 with ylim[1] < ylim[2].")
  }

  # ── Normalize and validate ODE ───────────────────────────────────────────
  norm <- .normalize_ode(deriv, system = system)
  .validate_ode(norm, system = system, parameters = parameters)

  # ── Compute grid ─────────────────────────────────────────────────────────
  grid <- .compute_nullcline_grid(norm, system, xlim, ylim,
                                   n_points, parameters)

  # ── Build layers ─────────────────────────────────────────────────────────
  layers <- list()

  if (system == "two.dim") {

    if (add_legend) {
      # Use geom_contour with a color aesthetic mapped to a factor so that
      # a legend is generated automatically
      nc_long <- rbind(
        data.frame(x = grid$x, y = grid$y, z = grid$f,
                   nullcline = "x-nullcline (dx/dt = 0)"),
        data.frame(x = grid$x, y = grid$y, z = grid$g,
                   nullcline = "y-nullcline (dy/dt = 0)")
      )

      # Define colors and linetypes as named vectors for scale_*_manual
      nc_colors    <- stats::setNames(c(x_color, y_color),
                                       c("x-nullcline (dx/dt = 0)",
                                         "y-nullcline (dy/dt = 0)"))
      nc_linetypes <- stats::setNames(c(x_linetype, y_linetype),
                                       c("x-nullcline (dx/dt = 0)",
                                         "y-nullcline (dy/dt = 0)"))

      layers <- list(
        ggplot2::geom_contour(
          data    = nc_long,
          mapping = ggplot2::aes(x        = .data$x,
                                 y        = .data$y,
                                 z        = .data$z,
                                 color    = .data$nullcline,
                                 linetype = .data$nullcline),
          breaks    = 0,
          linewidth = linewidth,
          na.rm     = TRUE
        ),
        ggplot2::scale_color_manual(
          name   = NULL,
          values = nc_colors
        ),
        ggplot2::scale_linetype_manual(
          name   = NULL,
          values = nc_linetypes
        ),
        .legend_theme(legend_position)
      )

    } else {
      # No legend: draw x- and y-nullclines as separate layers
      layers <- list(
        ggplot2::geom_contour(
          data    = grid,
          mapping = ggplot2::aes(x = .data$x, y = .data$y, z = .data$f),
          breaks    = 0,
          color     = x_color,
          linetype  = x_linetype,
          linewidth = linewidth,
          na.rm     = TRUE
        ),
        ggplot2::geom_contour(
          data    = grid,
          mapping = ggplot2::aes(x = .data$x, y = .data$y, z = .data$g),
          breaks    = 0,
          color     = y_color,
          linetype  = y_linetype,
          linewidth = linewidth,
          na.rm     = TRUE
        )
      )
    }

  } else {
    # 1D: draw horizontal lines at equilibrium y-values (where f = 0).
    # Use geom_contour on a 2D grid with a dummy x dimension.
    xs       <- seq(xlim[[1L]], xlim[[2L]], length.out = 10L)
    grid_2d  <- expand.grid(x = xs, y = grid$y)
    grid_2d$f <- rep(grid$f, each = 10L)

    layers <- list(
      ggplot2::geom_contour(
        data    = grid_2d,
        mapping = ggplot2::aes(x = .data$x, y = .data$y, z = .data$f),
        breaks    = 0,
        color     = x_color,
        linetype  = x_linetype,
        linewidth = linewidth,
        na.rm     = TRUE
      )
    )
  }

  layers
}
