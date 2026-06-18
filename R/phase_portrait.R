# phase_portrait.R
#
# gg_phase_portrait() — add a 1D phase line to an existing phase portrait.
#
# The phase line is drawn as a vertical line at a fixed x position
# (defaulting to the left edge of the plot), with:
#   (1) A vertical reference line (the phase line itself)
#   (2) Upward/downward arrows showing the sign of dy/dt
#   (3) Filled circles at stable equilibria
#   (4) Open circles at unstable equilibria
#   (5) Diamond symbols at semi-stable equilibria
#
# Drawing the phase line vertically at the plot edge cleanly separates
# it from the flow field arrows, which occupy the interior of the plot.


# ---------------------------------------------------------------------------
# .find_equilibria_1d() and .classify_equilibrium_1d()
# (unchanged from original — see inline below)
# ---------------------------------------------------------------------------

#' @keywords internal
.find_equilibria_1d <- function(deriv, ylim, parameters, n_points = 500L) {

  ys <- seq(ylim[[1L]], ylim[[2L]], length.out = n_points)
  fs <- vapply(ys, function(yi) {
    .eval_ode(deriv, t = 0, y = c(yi), parameters = parameters)[[1L]]
  }, numeric(1L))

  eq_rows <- list()

  for (i in seq_len(length(fs) - 1L)) {
    f_left  <- fs[[i]]
    f_right <- fs[[i + 1L]]

    if (f_left == 0) {
      stability <- .classify_equilibrium_1d(fs, i, n_points)
      eq_rows[[length(eq_rows) + 1L]] <- data.frame(
        y = ys[[i]], stability = stability
      )
    } else if (sign(f_left) != sign(f_right)) {
      y_eq <- ys[[i]] - f_left * (ys[[i + 1L]] - ys[[i]]) /
        (f_right - f_left)
      stability <- .classify_equilibrium_1d(fs, i, n_points)
      eq_rows[[length(eq_rows) + 1L]] <- data.frame(
        y = y_eq, stability = stability
      )
    }
  }

  if (length(eq_rows) == 0L) {
    return(data.frame(y = numeric(0L), stability = character(0L)))
  }

  do.call(rbind, eq_rows)
}

#' @keywords internal
.classify_equilibrium_1d <- function(fs, i, n_points) {
  left_sign  <- if (i > 1L)              sign(fs[[i - 1L]]) else sign(fs[[i]])
  right_sign <- if (i < n_points - 1L)   sign(fs[[i + 1L]]) else sign(fs[[i + 1L]])

  if      (left_sign > 0 && right_sign < 0) "stable"
  else if (left_sign < 0 && right_sign > 0) "unstable"
  else                                       "semi-stable"
}


# ---------------------------------------------------------------------------
# gg_phase_portrait()
# ---------------------------------------------------------------------------

#' Add a 1D phase line to a phase portrait
#'
#' Computes and adds a one-dimensional phase line to an existing
#' [ggplot2::ggplot()] object. The phase line is drawn as a **vertical
#' line at the left or right edge** of the plot (controlled by
#' `line_x_position`), cleanly separated from the flow field arrows in the
#' plot interior. It consists of:
#' \itemize{
#'   \item A vertical reference line at `line_x_position`
#'   \item Upward/downward arrows showing the sign of \eqn{dy/dt}
#'   \item Filled circles at stable (attracting) equilibria
#'   \item Open circles at unstable (repelling) equilibria
#'   \item Diamond symbols at semi-stable equilibria
#' }
#'
#' @param deriv A function describing the 1D ODE system, in Convention A
#'   or B. See [ggphasr] for details.
#' @param ylim Numeric vector of length 2. Range of the y-axis (state
#'   variable). Should match the `ylim` of the parent plot.
#' @param xlim Numeric vector of length 2. Range of the x-axis. Should
#'   match the `xlim` of the parent plot. Default: `c(0, 1)`.
#' @param parameters Parameter vector or list passed to `deriv`.
#' @param line_x_position Numeric or `NULL`. x-coordinate at which the
#'   vertical phase line is drawn. If `NULL` (default), placed at the
#'   left edge: `xlim[1] + 0.04 * diff(xlim)`.
#' @param n_arrows Integer. Number of directional arrows along the phase
#'   line. Default: `15`.
#' @param n_search Integer. Grid resolution for equilibrium detection.
#'   Default: `500`.
#' @param arrow_color Character. Color of the directional arrows.
#'   Default: `"grey30"`.
#' @param arrow_size Numeric. Size of arrow heads in lines. Default: `0.4`.
#' @param arrow_length_scale Numeric in (0, 1]. Arrow length as a fraction
#'   of the y-range divided by `n_arrows`. Default: `0.7`.
#' @param stable_color Character. Fill color for stable equilibria.
#'   Default: `"black"`.
#' @param unstable_fill Character. Fill color for unstable equilibria.
#'   Default: `"white"`.
#' @param eq_size Numeric. Size of equilibrium points. Default: `4`.
#' @param eq_stroke Numeric. Border width of equilibrium points.
#'   Default: `1`.
#' @param line_color Character. Color of the vertical phase line.
#'   Default: `"grey30"`.
#' @param line_linewidth Numeric. Width of the phase line. Default: `0.6`.
#'
#' @return A list of [ggplot2] layer objects. Add to a ggplot with `+`.
#'
#' @examples
#' # Logistic growth: stable equilibrium at K = 10, unstable at y = 0
#' gg_flow_field(
#'   ode_logistic,
#'   xlim = c(0, 6), ylim = c(-1, 12),
#'   system = "one.dim", parameters = c(r = 1, K = 10)
#' ) +
#'   gg_phase_portrait(
#'     ode_logistic,
#'     ylim = c(-1, 12), xlim = c(0, 6),
#'     parameters = c(r = 1, K = 10)
#'   )
#'
#' # Place phase line on the right edge instead
#' gg_flow_field(
#'   ode_example_01,
#'   xlim = c(0, 4), ylim = c(-4, 4),
#'   system = "one.dim"
#' ) +
#'   gg_phase_portrait(
#'     ode_example_01,
#'     ylim = c(-4, 4), xlim = c(0, 4),
#'     line_x_position = 3.85
#'   )
#'
#' @seealso [ggphasr::gg_flow_field()], [ggphasr::gg_nullclines()],
#'   [ggphasr::gg_time_series()]
#' @export
gg_phase_portrait <- function(deriv,
                               ylim,
                               xlim              = c(0, 1),
                               parameters        = NULL,
                               line_x_position   = NULL,
                               n_arrows          = 15L,
                               n_search          = 500L,
                               arrow_color       = "grey30",
                               arrow_size        = 0.4,
                               arrow_length_scale = 0.7,
                               stable_color      = "black",
                               unstable_fill     = "white",
                               eq_size           = 4,
                               eq_stroke         = 1,
                               line_color        = "grey30",
                               line_linewidth     = 0.6) {

  # ── Input validation ─────────────────────────────────────────────────────
  if (!is.numeric(ylim) || length(ylim) != 2L || ylim[[1L]] >= ylim[[2L]]) {
    rlang::abort("`ylim` must be a numeric vector of length 2 with ylim[1] < ylim[2].")
  }
  if (!is.numeric(xlim) || length(xlim) != 2L || xlim[[1L]] >= xlim[[2L]]) {
    rlang::abort("`xlim` must be a numeric vector of length 2 with xlim[1] < xlim[2].")
  }
  n_arrows <- as.integer(n_arrows)
  n_search <- as.integer(n_search)

  # Default: left edge with a small inset
  if (is.null(line_x_position)) {
    line_x_position <- xlim[[1L]] + 0.04 * diff(xlim)
  }

  # ── Normalize and validate ODE ───────────────────────────────────────────
  norm <- .normalize_ode(deriv, system = "one.dim")
  .validate_ode(norm, system = "one.dim", parameters = parameters)

  # ── Find equilibria ──────────────────────────────────────────────────────
  eq_data <- .find_equilibria_1d(norm, ylim, parameters, n_search)

  # ── Build directional arrow data ─────────────────────────────────────────
  # Evenly spaced y positions, trimming endpoints
  arrow_ys <- seq(ylim[[1L]], ylim[[2L]], length.out = n_arrows + 2L)
  arrow_ys <- arrow_ys[-c(1L, length(arrow_ys))]

  # Remove positions within 2% of ylim range of any equilibrium
  tol <- diff(ylim) * 0.02
  if (nrow(eq_data) > 0L) {
    too_close <- vapply(arrow_ys, function(ay) {
      any(abs(ay - eq_data$y) < tol)
    }, logical(1L))
    arrow_ys <- arrow_ys[!too_close]
  }

  # Evaluate dy/dt at each arrow y position
  arrow_fs <- vapply(arrow_ys, function(yi) {
    .eval_ode(norm, t = 0, y = c(yi), parameters = parameters)[[1L]]
  }, numeric(1L))

  # Arrow length: fraction of cell spacing along y
  cell_size  <- diff(ylim) / n_arrows
  arrow_len  <- cell_size * arrow_length_scale

  # Vertical arrows: x fixed at line_x_position, y moves up or down
  arrow_data <- data.frame(
    x    = line_x_position,
    y    = arrow_ys,
    xend = line_x_position,
    yend = arrow_ys + arrow_len * sign(arrow_fs)
  )
  # Remove zero-derivative arrows (at equilibria)
  arrow_data <- arrow_data[abs(arrow_fs) > .Machine$double.eps, ]

  # ── Build layers ─────────────────────────────────────────────────────────
  layers <- list()

  # Vertical phase line
  layers[[1L]] <- ggplot2::geom_vline(
    xintercept = line_x_position,
    color      = line_color,
    linewidth  = line_linewidth,
    linetype   = "solid"
  )

  # Directional arrows (vertical, up = dy/dt > 0, down = dy/dt < 0)
  if (nrow(arrow_data) > 0L) {
    layers[[2L]] <- ggplot2::geom_segment(
      data    = arrow_data,
      mapping = ggplot2::aes(x    = .data$x,
                             y    = .data$y,
                             xend = .data$xend,
                             yend = .data$yend),
      color     = arrow_color,
      linewidth = 0.5,
      arrow     = grid::arrow(length = grid::unit(arrow_size, "lines"),
                              type   = "closed")
    )
  }

  # Equilibrium points on the phase line
  if (nrow(eq_data) > 0L) {

    stable <- eq_data[eq_data$stability == "stable", ]
    if (nrow(stable) > 0L) {
      layers[[length(layers) + 1L]] <- ggplot2::geom_point(
        data    = data.frame(x = line_x_position, y = stable$y),
        mapping = ggplot2::aes(x = .data$x, y = .data$y),
        shape  = 21,
        fill   = stable_color,
        color  = stable_color,
        size   = eq_size,
        stroke = eq_stroke
      )
    }

    unstable <- eq_data[eq_data$stability == "unstable", ]
    if (nrow(unstable) > 0L) {
      layers[[length(layers) + 1L]] <- ggplot2::geom_point(
        data    = data.frame(x = line_x_position, y = unstable$y),
        mapping = ggplot2::aes(x = .data$x, y = .data$y),
        shape  = 21,
        fill   = unstable_fill,
        color  = stable_color,
        size   = eq_size,
        stroke = eq_stroke
      )
    }

    semi <- eq_data[eq_data$stability == "semi-stable", ]
    if (nrow(semi) > 0L) {
      layers[[length(layers) + 1L]] <- ggplot2::geom_point(
        data    = data.frame(x = line_x_position, y = semi$y),
        mapping = ggplot2::aes(x = .data$x, y = .data$y),
        shape  = 23,
        fill   = "grey50",
        color  = stable_color,
        size   = eq_size,
        stroke = eq_stroke
      )
    }
  }

  layers
}
