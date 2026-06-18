# time_series.R
#
# gg_time_series() — plot solution curves as state variable(s) vs. time.
#
# Produces a standalone ggplot object with:
#   - One facet per state variable for 2D systems
#   - A single panel for 1D systems
#   - Multiple initial conditions shown as different colored lines
#   - Optional time series for both forward and backward integration
#
# Uses the same .integrate_trajectory() and .parse_initial_conditions()
# infrastructure as gg_trajectory().


#' Plot time series of ODE solutions
#'
#' Numerically integrates a one- or two-dimensional autonomous ODE system
#' from one or more initial conditions and plots the state variable(s) as
#' functions of time. Returns a complete [ggplot2::ggplot()] object.
#'
#' For 2D systems, the plot is faceted: one panel shows \eqn{x(t)} and the
#' other shows \eqn{y(t)}, sharing a common time axis. This avoids y-axis
#' scaling problems when the two state variables have different magnitudes.
#'
#' @param deriv A function describing the ODE system, in Convention A or B.
#'   See [ggphasr] for details.
#' @param y0 Initial condition(s). A numeric vector (single IC), a matrix
#'   with one row per IC, or a list of numeric vectors. Same format as
#'   [ggphasr::gg_trajectory()].
#' @param t_end Numeric. End time for integration. Default: `10`.
#' @param system Character: `"two.dim"` (default) or `"one.dim"`.
#' @param parameters Parameter vector or list passed to `deriv`.
#' @param t_start_back Numeric or `NULL`. If supplied (negative), also
#'   integrates backward to this time. Default: `NULL`.
#' @param t_steps Integer. Number of time steps. Default: `500`.
#' @param method Character. deSolve integration method. Default: `"lsoda"`.
#' @param color Character or `NULL`. If `NULL` (default) and multiple ICs
#'   are supplied, each IC gets a distinct color. If a single color string,
#'   all lines share that color.
#' @param linewidth Numeric. Line width. Default: `0.7`.
#' @param xlab Character. x-axis (time) label. Default: `"Time"`.
#' @param ylab Character. y-axis label for 1D systems. Default: `"y"`.
#'   For 2D systems the facet strip labels are used instead.
#' @param var_labels Character vector of length 2. Labels for the two state
#'   variables in the facet strips of a 2D plot. Default: `c("x", "y")`.
#' @param title Character or `NULL`. Plot title. Default: `NULL`.
#' @param add_legend Logical. Whether to show a legend when multiple ICs
#'   are used. Default: `TRUE`.
#'
#' @return A [ggplot2::ggplot()] object.
#'
#' @examples
#' # 1D time series: logistic growth from several initial conditions
#' gg_time_series(
#'   ode_logistic,
#'   y0         = list(c(0.5), c(3), c(7), c(12)),
#'   t_end      = 8,
#'   system     = "one.dim",
#'   parameters = c(r = 1, K = 10),
#'   title      = "Logistic growth"
#' )
#'
#' # 2D time series: Lotka-Volterra — prey and predator in separate panels
#' gg_time_series(
#'   ode_lotka_volterra,
#'   y0         = list(c(1, 1), c(3, 2)),
#'   t_end      = 20,
#'   parameters = c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1),
#'   var_labels = c("Prey", "Predator"),
#'   title      = "Lotka-Volterra time series"
#' )
#'
#' # Combine with phase portrait to see the same trajectory two ways
#' lv_params <- c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1)
#'
#' # Phase plane view
#' p_phase <- gg_flow_field(ode_lotka_volterra,
#'                           xlim = c(0,5), ylim = c(0,5),
#'                           parameters = lv_params) +
#'   gg_trajectory(ode_lotka_volterra, y0 = c(1, 1),
#'                 xlim = c(0,5), ylim = c(0,5),
#'                 parameters = lv_params, color = "black")
#'
#' # Time domain view
#' p_time <- gg_time_series(ode_lotka_volterra, y0 = c(1, 1),
#'                           t_end = 20, parameters = lv_params)
#'
#' @seealso [ggphasr::gg_trajectory()], [ggphasr::gg_flow_field()]
#' @export
gg_time_series <- function(deriv,
                            y0,
                            t_end        = 10,
                            system       = c("two.dim", "one.dim"),
                            parameters   = NULL,
                            t_start_back = NULL,
                            t_steps      = 500L,
                            method       = "lsoda",
                            color        = NULL,
                            linewidth    = 0.7,
                            xlab         = "Time",
                            ylab         = "y",
                            var_labels   = c("x", "y"),
                            title        = NULL,
                            add_legend   = TRUE) {

  # ── Input validation ─────────────────────────────────────────────────────
  system  <- match.arg(system)
  t_steps <- as.integer(t_steps)

  if (t_end <= 0) {
    rlang::abort("`t_end` must be positive.")
  }
  if (!is.null(t_start_back) && t_start_back >= 0) {
    rlang::abort("`t_start_back` must be negative.")
  }
  if (system == "two.dim" && length(var_labels) != 2L) {
    rlang::abort("`var_labels` must be a character vector of length 2 for 2D systems.")
  }

  # ── Normalize ICs and ODE ────────────────────────────────────────────────
  y0_list <- .parse_initial_conditions(y0, system)
  norm    <- .normalize_ode(deriv, system = system)
  .validate_ode(norm, system = system, parameters = parameters)

  # ── Dummy xlim/ylim for .build_trajectory_data() ─────────────────────────
  # (not used for clipping in time series, but the helper requires them)
  dummy_xlim <- c(0, t_end)
  dummy_ylim <- c(-1e6, 1e6)

  # ── Integrate ────────────────────────────────────────────────────────────
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

  # ── Reshape for plotting ──────────────────────────────────────────────────
  # For 2D: pivot to long format with a `variable` column for faceting
  use_color_scale <- is.null(color) && length(y0_list) > 1L
  fixed_color     <- if (!use_color_scale) {
    if (is.null(color)) "black" else color
  } else NULL

  if (system == "two.dim") {

    # Long format: one row per (time, ic_id, direction, variable)
    long_x <- traj_data[, c("time", "x", "ic_id", "direction")]
    long_x$value    <- long_x$x
    long_x$variable <- var_labels[[1L]]
    long_x$x        <- NULL

    long_y <- traj_data[, c("time", "y", "ic_id", "direction")]
    long_y$value    <- long_y$y
    long_y$variable <- var_labels[[2L]]
    long_y$y        <- NULL

    # Rename x column in long_x (it was already dropped above via NULL)
    # Rebuild with consistent columns
    plot_data <- rbind(
      data.frame(time      = long_x$time,
                 value     = long_x$value,
                 ic_id     = long_x$ic_id,
                 direction = long_x$direction,
                 variable  = long_x$variable),
      data.frame(time      = long_y$time,
                 value     = long_y$value,
                 ic_id     = long_y$ic_id,
                 direction = long_y$direction,
                 variable  = long_y$variable)
    )

    # Keep facet order consistent with var_labels
    plot_data$variable <- factor(plot_data$variable,
                                  levels = var_labels)

    # ── Build 2D faceted plot ───────────────────────────────────────────
    p <- ggplot2::ggplot()

    if (use_color_scale) {
      p <- p +
        ggplot2::geom_line(
          data    = plot_data,
          mapping = ggplot2::aes(
            x     = .data$time,
            y     = .data$value,
            color = .data$ic_id,
            group = interaction(.data$ic_id, .data$direction)
          ),
          linewidth = linewidth
        ) +
        ggplot2::scale_color_brewer(palette = "Dark2", name = "IC")
    } else {
      p <- p +
        ggplot2::geom_line(
          data    = plot_data,
          mapping = ggplot2::aes(
            x     = .data$time,
            y     = .data$value,
            group = interaction(.data$ic_id, .data$direction)
          ),
          color     = fixed_color,
          linewidth = linewidth
        )
    }

    p <- p +
      ggplot2::facet_wrap(
        ggplot2::vars(.data$variable),
        ncol   = 1,
        scales = "free_y",
        strip.position = "left"
      ) +
      ggplot2::labs(x = xlab, y = NULL, title = title) +
      ggplot2::theme(
        strip.placement  = "outside",
        strip.background = ggplot2::element_blank(),
        strip.text       = ggplot2::element_text(size = ggplot2::rel(0.95))
      )

  } else {
    # ── 1D: single panel, time on x, y on y ─────────────────────────────
    p <- ggplot2::ggplot()

    if (use_color_scale) {
      p <- p +
        ggplot2::geom_line(
          data    = traj_data,
          mapping = ggplot2::aes(
            x     = .data$time,
            y     = .data$y,
            color = .data$ic_id,
            group = interaction(.data$ic_id, .data$direction)
          ),
          linewidth = linewidth
        ) +
        ggplot2::scale_color_brewer(palette = "Dark2", name = "IC")
    } else {
      p <- p +
        ggplot2::geom_line(
          data    = traj_data,
          mapping = ggplot2::aes(
            x     = .data$time,
            y     = .data$y,
            group = interaction(.data$ic_id, .data$direction)
          ),
          color     = fixed_color,
          linewidth = linewidth
        )
    }

    p <- p + ggplot2::labs(x = xlab, y = ylab, title = title)
  }

  # ── Legend ───────────────────────────────────────────────────────────────
  if (!add_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  # ── Theme ────────────────────────────────────────────────────────────────
  p <- p + theme_phase_plane()

  p
}
