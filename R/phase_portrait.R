# phase_portrait.R
#
# gg_phase_portrait() — add a 1D phase line to an existing phase portrait.
#
# The 1D phase line consists of:
#   (1) A horizontal reference line at a user-specified y-position
#       (the "phase line" itself, drawn across the x-axis range)
#   (2) Directional arrows along the line showing the sign of dy/dt
#   (3) Filled circles at stable equilibria (dy/dt changes - to +, i.e.
#       the derivative changes sign from negative to positive going left
#       to right — wait, stable means dy/dt > 0 below and < 0 above,
#       so attracting equilibria) and open circles at unstable equilibria
#
# For a 1D system this is designed to be added to a gg_flow_field() call
# with system = "one.dim", appearing as an overlay on the phase portrait.
# It can also be used standalone on a blank ggplot for the classic
# horizontal number-line presentation.
#
# Equilibrium classification for 1D:
#   Stable (attracting):   dy/dt > 0 just below, dy/dt < 0 just above
#   Unstable (repelling):  dy/dt < 0 just below, dy/dt > 0 just above
#   Semi-stable:           dy/dt has the same sign on both sides


# ---------------------------------------------------------------------------
# .find_equilibria_1d()
# ---------------------------------------------------------------------------
#
# Finds approximate equilibrium locations for a 1D ODE by scanning the
# derivative on a fine grid and detecting sign changes.
#
# @param deriv       Normalized (Convention A) 1D ODE function.
# @param ylim        Numeric vector of length 2: range to search.
# @param parameters  Parameter vector/list.
# @param n_points    Integer: grid resolution for sign-change detection.
# @return A data frame with columns: y, stability.
#   stability is one of "stable", "unstable", or "semi-stable".
#
.find_equilibria_1d <- function(deriv, ylim, parameters, n_points = 500L) {

  ys <- seq(ylim[[1L]], ylim[[2L]], length.out = n_points)
  fs <- vapply(ys, function(yi) {
    .eval_ode(deriv, t = 0, y = c(yi), parameters = parameters)[[1L]]
  }, numeric(1L))

  # Detect sign changes between adjacent points
  eq_rows <- list()

  for (i in seq_len(length(fs) - 1L)) {
    f_left  <- fs[[i]]
    f_right <- fs[[i + 1L]]

    # Sign change (or exact zero)
    if (f_left == 0) {
      # Exact zero at grid point i
      stability <- .classify_equilibrium_1d(fs, i, n_points)
      eq_rows[[length(eq_rows) + 1L]] <- data.frame(
        y         = ys[[i]],
        stability = stability
      )
    } else if (sign(f_left) != sign(f_right)) {
      # Linear interpolation to find the zero crossing
      y_eq <- ys[[i]] - f_left * (ys[[i + 1L]] - ys[[i]]) /
        (f_right - f_left)
      stability <- .classify_equilibrium_1d(fs, i, n_points)
      eq_rows[[length(eq_rows) + 1L]] <- data.frame(
        y         = y_eq,
        stability = stability
      )
    }
  }

  if (length(eq_rows) == 0L) {
    return(data.frame(y = numeric(0L), stability = character(0L)))
  }

  do.call(rbind, eq_rows)
}


# ---------------------------------------------------------------------------
# .classify_equilibrium_1d()
# ---------------------------------------------------------------------------
#
# Classifies a 1D equilibrium at index i in the derivative vector fs.
#
# @param fs        Numeric vector of derivative values.
# @param i         Integer index of the sign change.
# @param n_points  Total grid length (for boundary checking).
# @return Character: "stable", "unstable", or "semi-stable".
#
.classify_equilibrium_1d <- function(fs, i, n_points) {
  left_sign  <- if (i > 1L)          sign(fs[[i - 1L]]) else sign(fs[[i]])
  right_sign <- if (i < n_points - 1L) sign(fs[[i + 1L]]) else sign(fs[[i + 1L]])

  if (left_sign > 0 && right_sign < 0) {
    "stable"
  } else if (left_sign < 0 && right_sign > 0) {
    "unstable"
  } else {
    "semi-stable"
  }
}


# ---------------------------------------------------------------------------
# gg_phase_portrait()
# ---------------------------------------------------------------------------

#' Add a 1D phase line to a phase portrait
#'
#' Computes and adds a one-dimensional phase line to an existing
#' [ggplot2::ggplot()] object. The phase line consists of:
#' \itemize{
#'   \item A horizontal reference line at `line_position` on the y-axis
#'   \item Directional arrows along the line showing the sign of \eqn{dy/dt}
#'   \item Filled circles at stable (attracting) equilibria
#'   \item Open circles at unstable (repelling) equilibria
#'   \item Half-filled circles at semi-stable equilibria
#' }
#'
#' Designed to be added to a [ggphasr::gg_flow_field()] call with
#' `system = "one.dim"`:
#'
#' ```r
#' gg_flow_field(my_ode, xlim = c(0, 5), ylim = c(-1, 12),
#'               system = "one.dim") +
#'   gg_phase_portrait(my_ode, ylim = c(-1, 12))
#' ```
#'
#' @param deriv A function describing the 1D ODE system, in Convention A
#'   (`f(t, y, parameters)`) or Convention B (`f(y, parameters)`).
#' @param ylim Numeric vector of length 2. Range of the y-axis (state
#'   variable). Should match the `ylim` of the parent plot.
#' @param xlim Numeric vector of length 2. Range of the x-axis. Used to
#'   determine arrow placement. Should match the `xlim` of the parent
#'   plot. Default: `c(0, 1)`.
#' @param parameters Parameter vector or list passed to `deriv`.
#' @param line_position Numeric. y-coordinate at which to draw the phase
#'   line. Default: `0` (on the x-axis).
#' @param n_arrows Integer. Number of directional arrows along the phase
#'   line. Default: `15`.
#' @param n_search Integer. Grid resolution for equilibrium detection.
#'   Default: `500`. Increase if equilibria are missed.
#' @param arrow_color Character. Color of the directional arrows.
#'   Default: `"grey40"`.
#' @param arrow_size Numeric. Size of arrow heads in lines. Default: `0.4`.
#' @param stable_color Character. Fill color for stable equilibrium points.
#'   Default: `"black"`.
#' @param unstable_color Character. Fill color for unstable equilibrium
#'   points. Default: `"white"`.
#' @param eq_size Numeric. Size of equilibrium points. Default: `4`.
#' @param eq_stroke Numeric. Border width of equilibrium points.
#'   Default: `1`.
#' @param line_color Character. Color of the phase line itself.
#'   Default: `"grey40"`.
#' @param line_linewidth Numeric. Width of the phase line. Default: `0.6`.
#'
#' @return A list of [ggplot2] layer objects. Add to a ggplot with `+`.
#'
#' @examples
#' # Logistic growth: stable equilibrium at K=10, unstable at 0
#' gg_flow_field(
#'   ode_logistic,
#'   xlim       = c(0, 6),
#'   ylim       = c(-1, 12),
#'   system     = "one.dim",
#'   parameters = c(r = 1, K = 10)
#' ) +
#'   gg_phase_portrait(
#'     ode_logistic,
#'     ylim       = c(-1, 12),
#'     xlim       = c(0, 6),
#'     parameters = c(r = 1, K = 10)
#'   )
#'
#' # Example 01: dy/dt = 4 - y^2, equilibria at y = +-2
#' gg_flow_field(
#'   ode_example_01,
#'   xlim   = c(0, 4),
#'   ylim   = c(-4, 4),
#'   system = "one.dim"
#' ) +
#'   gg_phase_portrait(
#'     ode_example_01,
#'     ylim = c(-4, 4),
#'     xlim = c(0, 4)
#'   )
#'
#' @seealso [ggphasr::gg_flow_field()], [ggphasr::gg_nullclines()],
#'   [ggphasr::gg_time_series()]
#' @export
gg_phase_portrait <- function(deriv,
                               ylim,
                               xlim             = c(0, 1),
                               parameters       = NULL,
                               line_position    = 0,
                               n_arrows         = 15L,
                               n_search         = 500L,
                               arrow_color      = "grey40",
                               arrow_size       = 0.4,
                               stable_color     = "black",
                               unstable_color   = "white",
                               eq_size          = 4,
                               eq_stroke        = 1,
                               line_color       = "grey40",
                               line_linewidth    = 0.6) {

  # ── Input validation ─────────────────────────────────────────────────────
  if (!is.numeric(ylim) || length(ylim) != 2L || ylim[[1L]] >= ylim[[2L]]) {
    rlang::abort("`ylim` must be a numeric vector of length 2 with ylim[1] < ylim[2].")
  }
  if (!is.numeric(xlim) || length(xlim) != 2L || xlim[[1L]] >= xlim[[2L]]) {
    rlang::abort("`xlim` must be a numeric vector of length 2 with xlim[1] < xlim[2].")
  }
  n_arrows <- as.integer(n_arrows)
  n_search <- as.integer(n_search)

  # ── Normalize and validate ODE ───────────────────────────────────────────
  norm <- .normalize_ode(deriv, system = "one.dim")
  .validate_ode(norm, system = "one.dim", parameters = parameters)

  # ── Find equilibria ──────────────────────────────────────────────────────
  eq_data <- .find_equilibria_1d(norm, ylim, parameters, n_search)

  # ── Build directional arrow data ─────────────────────────────────────────
  # Place n_arrows evenly along the y range; skip points too close to equil.
  arrow_ys  <- seq(ylim[[1L]], ylim[[2L]], length.out = n_arrows + 2L)
  arrow_ys  <- arrow_ys[-c(1L, length(arrow_ys))]  # trim endpoints

  # Remove arrow positions within 2% of ylim range of any equilibrium
  tol <- diff(ylim) * 0.02
  if (nrow(eq_data) > 0L) {
    too_close <- vapply(arrow_ys, function(ay) {
      any(abs(ay - eq_data$y) < tol)
    }, logical(1L))
    arrow_ys <- arrow_ys[!too_close]
  }

  # Evaluate dy/dt at each arrow position
  arrow_fs <- vapply(arrow_ys, function(yi) {
    .eval_ode(norm, t = 0, y = c(yi), parameters = parameters)[[1L]]
  }, numeric(1L))

  # Arrow segment: x fixed at midpoint of xlim, y to y + scaled_f
  x_mid     <- mean(xlim)
  x_range   <- diff(xlim)
  # Scale arrow length to 8% of x range
  max_f     <- max(abs(arrow_fs), na.rm = TRUE)
  if (max_f < .Machine$double.eps) max_f <- 1
  arrow_len <- x_range * 0.08

  arrow_data <- data.frame(
    x    = x_mid,
    y    = arrow_ys,
    xend = x_mid + arrow_len * sign(arrow_fs),
    yend = arrow_ys
  )
  # Remove zero-derivative arrows
  arrow_data <- arrow_data[abs(arrow_fs) > .Machine$double.eps, ]

  # ── Build layers ─────────────────────────────────────────────────────────
  layers <- list()

  # Phase line (horizontal reference)
  layers[[1L]] <- ggplot2::geom_hline(
    yintercept = line_position,
    color      = line_color,
    linewidth  = line_linewidth,
    linetype   = "solid"
  )

  # Directional arrows
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

  # Equilibrium points
  if (nrow(eq_data) > 0L) {
    # Stable: filled circle
    stable <- eq_data[eq_data$stability == "stable", ]
    if (nrow(stable) > 0L) {
      layers[[length(layers) + 1L]] <- ggplot2::geom_point(
        data    = data.frame(x = x_mid, y = stable$y),
        mapping = ggplot2::aes(x = .data$x, y = .data$y),
        shape  = 21,
        fill   = stable_color,
        color  = stable_color,
        size   = eq_size,
        stroke = eq_stroke
      )
    }

    # Unstable: open circle
    unstable <- eq_data[eq_data$stability == "unstable", ]
    if (nrow(unstable) > 0L) {
      layers[[length(layers) + 1L]] <- ggplot2::geom_point(
        data    = data.frame(x = x_mid, y = unstable$y),
        mapping = ggplot2::aes(x = .data$x, y = .data$y),
        shape  = 21,
        fill   = unstable_color,
        color  = stable_color,
        size   = eq_size,
        stroke = eq_stroke
      )
    }

    # Semi-stable: half-filled (use a triangle point as approximation)
    semi <- eq_data[eq_data$stability == "semi-stable", ]
    if (nrow(semi) > 0L) {
      layers[[length(layers) + 1L]] <- ggplot2::geom_point(
        data    = data.frame(x = x_mid, y = semi$y),
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
