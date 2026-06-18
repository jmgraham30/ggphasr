# trajectory.R
#
# gg_trajectory() — numerically integrate and plot solution trajectories
# for a 1D or 2D autonomous ODE system.
#
# Uses deSolve::ode() for numerical integration (default method: "lsoda",
# which automatically selects between stiff and non-stiff solvers).
#
# Supports:
#   - Multiple initial conditions in a single call
#   - Forward and/or backward integration
#   - Arrow heads showing direction of flow along the trajectory
#   - Color mapping by initial condition index or a user-supplied variable
#
# Internal helpers:
#   .integrate_trajectory()   — integrates one IC forward or backward
#   .build_trajectory_data()  — integrates all ICs, returns tidy data frame


# ---------------------------------------------------------------------------
# .integrate_trajectory()
# ---------------------------------------------------------------------------
#
# Integrates the ODE from a single initial condition in one direction.
#
# @param deriv       Normalized (Convention A) ODE function.
# @param y0          Numeric vector: initial condition.
# @param t_start     Numeric: start time (always 0 for autonomous systems).
# @param t_end       Numeric: end time (positive = forward, negative = backward).
# @param t_steps     Integer: number of output time steps.
# @param parameters  Parameter vector/list passed to `deriv`.
# @param method      deSolve integration method. Default: "lsoda".
# @return A data frame with columns: time, x (or y for 1D), and y (2D only).
#   Returns NULL if integration fails.
#
.integrate_trajectory <- function(deriv,
                                   y0,
                                   t_start,
                                   t_end,
                                   t_steps,
                                   parameters,
                                   method = "lsoda") {

  times <- seq(t_start, t_end, length.out = t_steps)

  result <- tryCatch(
    deSolve::ode(
      y     = y0,
      times = times,
      func  = deriv,
      parms = parameters,
      method = method
    ),
    error   = function(e) NULL,
    warning = function(w) NULL
  )

  if (is.null(result)) return(NULL)

  df <- as.data.frame(result)

  if (ncol(df) == 2L) {
    # 1D system: columns are time, y1
    names(df) <- c("time", "y")
    df$x <- df$time   # for plotting, x-axis is time
  } else {
    # 2D system: columns are time, y1, y2
    names(df) <- c("time", "x", "y")
  }

  df
}


# ---------------------------------------------------------------------------
# .build_trajectory_data()
# ---------------------------------------------------------------------------
#
# Integrates all initial conditions (forward and/or backward) and assembles
# a single tidy data frame with an `ic_id` column identifying each trajectory.
#
# @param deriv       Normalized (Convention A) ODE function.
# @param system      "two.dim" or "one.dim".
# @param y0_list     List of numeric vectors, one per initial condition.
# @param t_end       Numeric: forward integration end time.
# @param t_start_back Numeric: backward integration end time (negative).
#   NULL means no backward integration.
# @param t_steps     Integer: number of output time steps per direction.
# @param parameters  Parameter vector/list.
# @param method      deSolve integration method.
# @return A data frame with columns: time, x, y, ic_id, direction.
#   ic_id is an integer identifying the initial condition.
#   direction is "forward" or "backward".
#
.build_trajectory_data <- function(deriv,
                                    system,
                                    y0_list,
                                    t_end,
                                    t_start_back,
                                    t_steps,
                                    parameters,
                                    method) {

  all_segs <- list()

  for (i in seq_along(y0_list)) {
    y0 <- y0_list[[i]]

    # Forward integration
    fwd <- .integrate_trajectory(deriv, y0, 0, t_end, t_steps,
                                  parameters, method)
    if (!is.null(fwd)) {
      fwd$ic_id     <- i
      fwd$direction <- "forward"
      all_segs[[length(all_segs) + 1L]] <- fwd
    }

    # Backward integration (if requested)
    if (!is.null(t_start_back)) {
      bwd <- .integrate_trajectory(deriv, y0, 0, t_start_back, t_steps,
                                    parameters, method)
      if (!is.null(bwd)) {
        bwd$ic_id     <- i
        bwd$direction <- "backward"
        all_segs[[length(all_segs) + 1L]] <- bwd
      }
    }
  }

  if (length(all_segs) == 0L) {
    rlang::abort(
      "All trajectory integrations failed. ",
      "Check that `deriv`, `parameters`, and the initial conditions are valid."
    )
  }

  do.call(rbind, all_segs)
}


# ---------------------------------------------------------------------------
# gg_trajectory()
# ---------------------------------------------------------------------------

#' Add solution trajectories to a phase plane plot
#'
#' Numerically integrates and plots one or more solution trajectories of a
#' one- or two-dimensional autonomous ODE system. Returns a list of
#' [ggplot2] layer objects that can be added to a [ggplot2::ggplot()] object
#' with `+`.
#'
#' @param deriv A function describing the ODE system, in Convention A or B.
#'   See [ggphasr] for details.
#' @param y0 Initial condition(s). One of:
#'   \itemize{
#'     \item A numeric vector of length 1 (1D) or 2 (2D) for a single
#'       initial condition.
#'     \item A numeric matrix with one row per initial condition (columns
#'       are state variables).
#'     \item A list of numeric vectors, one per initial condition.
#'   }
#' @param xlim Numeric vector of length 2. x-axis range. Used to clip
#'   trajectories that leave the plot area.
#' @param ylim Numeric vector of length 2. y-axis range.
#' @param system Character: `"two.dim"` (default) or `"one.dim"`.
#' @param parameters Parameter vector or list passed to `deriv`.
#' @param t_end Numeric. End time for forward integration. Default: `10`.
#' @param t_start_back Numeric or `NULL`. End time for backward integration
#'   (should be negative). If `NULL` (default), only forward integration
#'   is performed. Set to a negative value (e.g., `-10`) to also integrate
#'   backward from each initial condition.
#' @param t_steps Integer. Number of time steps per integration direction.
#'   Default: `500`. Increase for smoother curves on long integrations.
#' @param method Character. deSolve integration method. Default: `"lsoda"`
#'   (automatically switches between stiff and non-stiff solvers).
#' @param color Character or `NULL`. Trajectory color. If `NULL` and there
#'   are multiple initial conditions, each trajectory gets a different color
#'   from a discrete palette. If a single color string (e.g., `"black"`),
#'   all trajectories share that color. Default: `NULL`.
#' @param linewidth Numeric. Trajectory line width. Default: `0.7`.
#' @param add_arrows Logical. If `TRUE` (default), adds an arrow head at the
#'   midpoint of each trajectory segment showing the direction of flow.
#' @param arrow_size Numeric. Size of the direction arrow heads in lines.
#'   Default: `0.3`.
#' @param add_start_point Logical. If `TRUE` (default), marks each initial
#'   condition with a filled circle.
#' @param start_point_size Numeric. Size of the initial condition point.
#'   Default: `2`.
#'
#' @return A list of [ggplot2] layer objects. Add to a ggplot with `+`.
#'
#' @examples
#' # Single trajectory on a Lotka-Volterra phase plane
#' lv_params <- c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1)
#'
#' gg_flow_field(ode_lotka_volterra,
#'               xlim = c(0, 5), ylim = c(0, 5),
#'               parameters = lv_params) +
#'   gg_trajectory(ode_lotka_volterra,
#'                 y0         = c(1, 1),
#'                 xlim       = c(0, 5),
#'                 ylim       = c(0, 5),
#'                 parameters = lv_params)
#'
#' # Multiple initial conditions from a matrix
#' ics <- matrix(c(0.5, 0.5,
#'                 1.0, 2.0,
#'                 3.0, 1.0), ncol = 2, byrow = TRUE)
#'
#' gg_flow_field(ode_lotka_volterra,
#'               xlim = c(0, 5), ylim = c(0, 5),
#'               parameters = lv_params) +
#'   gg_nullclines(ode_lotka_volterra,
#'                 xlim = c(0, 5), ylim = c(0, 5),
#'                 parameters = lv_params) +
#'   gg_trajectory(ode_lotka_volterra,
#'                 y0         = ics,
#'                 xlim       = c(0, 5),
#'                 ylim       = c(0, 5),
#'                 parameters = lv_params,
#'                 t_end      = 20)
#'
#' # Forward and backward integration near a saddle point
#' gg_flow_field(ode_example_08,
#'               xlim = c(-3, 3), ylim = c(-3, 3)) +
#'   gg_trajectory(ode_example_08,
#'                 y0           = c(0.1, 2),
#'                 xlim         = c(-3, 3),
#'                 ylim         = c(-3, 3),
#'                 t_end        = 3,
#'                 t_start_back = -3)
#'
#' @seealso [ggphasr::gg_flow_field()], [ggphasr::gg_nullclines()]
#' @export
gg_trajectory <- function(deriv,
                           y0,
                           xlim,
                           ylim,
                           system          = c("two.dim", "one.dim"),
                           parameters      = NULL,
                           t_end           = 10,
                           t_start_back    = NULL,
                           t_steps         = 500L,
                           method          = "lsoda",
                           color           = NULL,
                           linewidth       = 0.7,
                           add_arrows      = TRUE,
                           arrow_size      = 0.3,
                           add_start_point = TRUE,
                           start_point_size = 2) {

  # ── Input validation ─────────────────────────────────────────────────────
  system  <- match.arg(system)
  t_steps <- as.integer(t_steps)

  if (!is.numeric(xlim) || length(xlim) != 2L || xlim[[1L]] >= xlim[[2L]]) {
    rlang::abort("`xlim` must be a numeric vector of length 2 with xlim[1] < xlim[2].")
  }
  if (!is.numeric(ylim) || length(ylim) != 2L || ylim[[1L]] >= ylim[[2L]]) {
    rlang::abort("`ylim` must be a numeric vector of length 2 with ylim[1] < ylim[2].")
  }
  if (t_end <= 0) {
    rlang::abort("`t_end` must be positive.")
  }
  if (!is.null(t_start_back) && t_start_back >= 0) {
    rlang::abort("`t_start_back` must be negative (backward integration goes to negative time).")
  }

  # ── Normalize initial conditions to a list ───────────────────────────────
  y0_list <- .parse_initial_conditions(y0, system)

  # ── Normalize and validate ODE ───────────────────────────────────────────
  norm <- .normalize_ode(deriv, system = system)
  .validate_ode(norm, system = system, parameters = parameters)

  # ── Integrate all trajectories ───────────────────────────────────────────
  traj_data <- .build_trajectory_data(
    deriv        = norm,
    system       = system,
    y0_list      = y0_list,
    t_end        = t_end,
    t_start_back = t_start_back,
    t_steps      = t_steps,
    parameters   = parameters,
    method       = method
  )

  traj_data$ic_id <- as.factor(traj_data$ic_id)

  # ── Determine color mapping ──────────────────────────────────────────────
  use_color_scale <- is.null(color) && length(y0_list) > 1L
  fixed_color     <- if (!use_color_scale) {
    if (is.null(color)) "black" else color
  } else NULL

  # ── Build layers ─────────────────────────────────────────────────────────
  layers <- list()

  arrow_spec <- if (add_arrows) {
    grid::arrow(length = grid::unit(arrow_size, "lines"),
                type   = "open",
                ends   = "last")
  } else {
    NULL
  }

  if (use_color_scale) {
    layers[[1L]] <- ggplot2::geom_path(
      data    = traj_data,
      mapping = ggplot2::aes(x     = .data$x,
                             y     = .data$y,
                             group = interaction(.data$ic_id,
                                                 .data$direction),
                             color = .data$ic_id),
      linewidth = linewidth,
      arrow     = arrow_spec,
      na.rm     = TRUE,
      lineend   = "round"
    )
    layers[[2L]] <- ggplot2::scale_color_brewer(
      palette = "Dark2",
      name    = "IC"
    )
  } else {
    layers[[1L]] <- ggplot2::geom_path(
      data    = traj_data,
      mapping = ggplot2::aes(x     = .data$x,
                             y     = .data$y,
                             group = interaction(.data$ic_id,
                                                 .data$direction)),
      color     = fixed_color,
      linewidth = linewidth,
      arrow     = arrow_spec,
      na.rm     = TRUE,
      lineend   = "round"
    )
  }

  # ── Initial condition markers ────────────────────────────────────────────
  if (add_start_point) {
    # Extract the first row of each ic_id (the initial condition point)
    ic_points <- do.call(rbind, lapply(y0_list, function(ic) {
      if (system == "one.dim") {
        data.frame(x = 0, y = ic[[1L]])
      } else {
        data.frame(x = ic[[1L]], y = ic[[2L]])
      }
    }))
    ic_points$ic_id <- as.factor(seq_len(nrow(ic_points)))

    if (use_color_scale) {
      layers[[length(layers) + 1L]] <- ggplot2::geom_point(
        data    = ic_points,
        mapping = ggplot2::aes(x     = .data$x,
                               y     = .data$y,
                               color = .data$ic_id),
        size  = start_point_size,
        na.rm = TRUE
      )
    } else {
      layers[[length(layers) + 1L]] <- ggplot2::geom_point(
        data    = ic_points,
        mapping = ggplot2::aes(x = .data$x, y = .data$y),
        color = fixed_color,
        size  = start_point_size,
        na.rm = TRUE
      )
    }
  }

  layers
}


# ---------------------------------------------------------------------------
# .parse_initial_conditions()
# ---------------------------------------------------------------------------
#
# Normalizes the `y0` argument of gg_trajectory() to a list of numeric
# vectors, regardless of whether the user passed a vector, matrix, or list.
#
# @param y0     Numeric vector, matrix, or list.
# @param system "two.dim" or "one.dim".
# @return A list of numeric vectors.
#
.parse_initial_conditions <- function(y0, system) {

  expected_len <- if (system == "one.dim") 1L else 2L

  if (is.matrix(y0)) {
    # Each row is one initial condition
    if (ncol(y0) != expected_len) {
      rlang::abort(
        paste0(
          "When `y0` is a matrix, it must have ", expected_len,
          " column(s) (one per state variable). Got ", ncol(y0), "."
        )
      )
    }
    return(lapply(seq_len(nrow(y0)), function(i) as.numeric(y0[i, ])))
  }

  if (is.list(y0)) {
    bad <- which(vapply(y0, length, integer(1L)) != expected_len)
    if (length(bad) > 0L) {
      rlang::abort(
        paste0(
          "Each initial condition in `y0` must have length ", expected_len,
          ". Element(s) ", paste(bad, collapse = ", "), " do not."
        )
      )
    }
    return(lapply(y0, as.numeric))
  }

  # Plain numeric vector: single initial condition
  if (!is.numeric(y0) || length(y0) != expected_len) {
    rlang::abort(
      paste0(
        "For `system = \"", system, "\"`, a single initial condition must be ",
        "a numeric vector of length ", expected_len,
        ". Got length ", length(y0), "."
      )
    )
  }

  list(as.numeric(y0))
}
