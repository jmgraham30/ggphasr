# phase_plane.R
#
# gg_phase_plane() — high-level wrapper that produces a complete phase plane
# analysis in a single call.
#
# By default produces:
#   - Flow field (gg_flow_field)
#   - Nullclines (gg_nullclines)
#   - Trajectories from a default grid of initial conditions (gg_trajectory)
#   - Equilibrium points annotated with stability type (if find_equilibria=TRUE)
#
# Returns a named list:
#   $plot          — the ggplot object
#   $equilibria    — data frame from classify_equilibrium() (or NULL)
#
# The plot component is a standard ggplot object and can be further
# customized with + as usual.


# ---------------------------------------------------------------------------
# .default_ics()
# ---------------------------------------------------------------------------
#
# Generates a default grid of initial conditions for gg_phase_plane().
# Places ICs at the intersections of a coarse grid, excluding points
# very close to the axes (which are often equilibria or boundaries).
#
# @param xlim     Numeric vector of length 2.
# @param ylim     Numeric vector of length 2.
# @param n_ic     Integer: number of ICs along each axis. Default: 4.
# @param system   "two.dim" or "one.dim".
# @return A list of numeric vectors.
#
.default_ics <- function(xlim, ylim, n_ic, system) {

  if (system == "one.dim") {
    # For 1D: evenly spaced initial y values, x (time) always starts at 0
    ys <- seq(ylim[[1L]] + diff(ylim) * 0.1,
              ylim[[2L]] - diff(ylim) * 0.1,
              length.out = n_ic)
    return(lapply(ys, function(yi) c(yi)))
  }

  # 2D: grid of (x, y) points with a small inset from the boundaries
  xs <- seq(xlim[[1L]] + diff(xlim) * 0.1,
            xlim[[2L]] - diff(xlim) * 0.1,
            length.out = n_ic)
  ys <- seq(ylim[[1L]] + diff(ylim) * 0.1,
            ylim[[2L]] - diff(ylim) * 0.1,
            length.out = n_ic)

  grid <- expand.grid(x = xs, y = ys)
  lapply(seq_len(nrow(grid)), function(i) c(grid$x[[i]], grid$y[[i]]))
}


# ---------------------------------------------------------------------------
# .equilibrium_shapes()
# ---------------------------------------------------------------------------
#
# Returns a named vector mapping classification strings to ggplot2 point
# shapes, for use with scale_shape_manual() in equilibrium annotation.
#
.equilibrium_shapes <- function() {
  c(
    "Stable node"              = 21L,
    "Unstable node"            = 21L,
    "Stable spiral"            = 21L,
    "Unstable spiral"          = 21L,
    "Center"                   = 22L,
    "Saddle"                   = 23L,
    "Non-isolated equilibrium" = 24L,
    "Stable"                   = 21L,
    "Unstable"                 = 21L,
    "Inconclusive (df/dy = 0)" = 22L
  )
}

.equilibrium_fills <- function() {
  c(
    "Stable node"              = "black",
    "Unstable node"            = "white",
    "Stable spiral"            = "black",
    "Unstable spiral"          = "white",
    "Center"                   = "white",
    "Saddle"                   = "grey60",
    "Non-isolated equilibrium" = "grey80",
    "Stable"                   = "black",
    "Unstable"                 = "white",
    "Inconclusive (df/dy = 0)" = "grey60"
  )
}


# ---------------------------------------------------------------------------
# gg_phase_plane()
# ---------------------------------------------------------------------------

#' Complete phase plane analysis in a single call
#'
#' A high-level wrapper that produces a complete phase plane portrait for a
#' one- or two-dimensional autonomous ODE system. By default generates a
#' flow field, nullclines, and trajectories from an evenly-spaced grid of
#' initial conditions. Optionally finds, classifies, and annotates all
#' equilibria automatically.
#'
#' Returns a named list so that both the plot and the equilibrium table are
#' immediately accessible:
#'
#' ```r
#' result <- gg_phase_plane(ode_lotka_volterra, ...)
#' result$plot        # the ggplot object
#' result$equilibria  # data frame of classified equilibria
#' ```
#'
#' The `$plot` component is a standard [ggplot2::ggplot()] object and can
#' be further customized with `+`.
#'
#' @param deriv A function describing the ODE system, in Convention A or B.
#'   See [ggphasr-package] for details.
#' @param xlim Numeric vector of length 2. x-axis range. Required.
#' @param ylim Numeric vector of length 2. y-axis range. Required.
#' @param system Character: `"two.dim"` (default) or `"one.dim"`.
#' @param parameters Parameter vector or list passed to `deriv`.
#'
#'
#' @param n_points Integer. Flow field grid resolution. Default: `21`.
#' @param arrow_type Character: `"equal"` (default) or `"proportional"`.
#' @param arrow_color Character. Flow field arrow color. Default: `"grey70"`.
#'
#'
#' @param show_nullclines Logical. Whether to draw nullclines.
#'   Default: `TRUE`.
#' @param nullcline_n_points Integer. Nullcline grid resolution.
#'   Default: `250`.
#'
#'
#' @param y0 Initial condition(s) for trajectories. A numeric vector,
#'   matrix, or list as accepted by [ggphasr::gg_trajectory()]. If `NULL`
#'   (default), a regular grid of `n_ic x n_ic` initial conditions is
#'   used automatically.
#' @param show_trajectories Logical. Whether to draw trajectories.
#'   Default: `TRUE`.
#' @param n_ic Integer. Number of auto-generated initial conditions per
#'   axis (ignored when `y0` is supplied). Default: `4` (giving 16 ICs
#'   for 2D systems, 4 for 1D).
#' @param t_end Numeric. Forward integration time. Default: `10`.
#' @param t_start_back Numeric or `NULL`. Backward integration time.
#'   Default: `NULL`.
#' @param trajectory_color Character or `NULL`. If `NULL` (default) and
#'   multiple ICs are used, each trajectory gets a distinct color. If a
#'   color string, all trajectories share that color.
#'
#'
#' @param find_equilibria Logical. Whether to automatically find, classify,
#'   and annotate equilibria. Default: `TRUE`. Uses a grid search over
#'   `xlim` x `ylim`.
#' @param eq_n_grid Integer. Grid resolution for equilibrium search.
#'   Default: `10`.
#' @param eq_grid_y0 List or `NULL`. Custom starting points for the
#'   equilibrium search. If `NULL`, a regular grid is used.
#' @param legend_position Character string or numeric vector of length 2.
#'   Controls the position of all legends (equilibrium types, nullclines)
#'   in the plot. One of `"right"` (default), `"left"`, `"top"`,
#'   `"bottom"`, `"none"`, `"inside"` (compact legend in the top-right
#'   corner of the panel), or a numeric vector `c(x, y)` with values
#'   in `[0, 1]` for a custom inside position. Passing `"inside"` is
#'   the most effective way to reclaim plot space when the external
#'   legend is too large.
#'
#'
#' @param xlab Character. x-axis label. Default: `"x"` (2D) or `"t"` (1D).
#' @param ylab Character. y-axis label. Default: `"y"`.
#' @param title Character or `NULL`. Plot title. Default: `NULL`.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{`plot`}{A [ggplot2::ggplot()] object.}
#'     \item{`equilibria`}{A data frame of classified equilibria (from
#'       [ggphasr::classify_equilibrium()]), or `NULL` if
#'       `find_equilibria = FALSE` or no equilibria were found.}
#'   }
#'
#' @examples
#' # Minimal call: produces everything automatically
#' result <- gg_phase_plane(
#'   ode_lotka_volterra,
#'   xlim       = c(0, 5),
#'   ylim       = c(0, 5),
#'   parameters = c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1)
#' )
#' result$plot
#' result$equilibria[, c("x", "y", "classification")]
#'
#' # 1D system
#' result <- gg_phase_plane(
#'   ode_logistic,
#'   xlim       = c(0, 8),
#'   ylim       = c(-1, 12),
#'   system     = "one.dim",
#'   parameters = c(r = 1, K = 10)
#' )
#' result$plot
#'
#' # Suppress equilibrium search for speed
#' result <- gg_phase_plane(
#'   ode_van_der_pol,
#'   xlim            = c(-3, 3),
#'   ylim            = c(-4, 4),
#'   parameters      = c(mu = 1),
#'   find_equilibria = FALSE,
#'   t_end           = 20
#' )
#' result$plot
#'
#' # Supply custom initial conditions
#' ics <- matrix(c(0.5,0.5, 1,2, 3,1, 2,3), ncol=2, byrow=TRUE)
#' result <- gg_phase_plane(
#'   ode_lotka_volterra,
#'   xlim       = c(0, 5),
#'   ylim       = c(0, 5),
#'   parameters = c(alpha=1, beta=0.5, delta=0.5, gamma=1),
#'   y0         = ics
#' )
#' result$plot
#'
#' # Further customize the returned plot
#' result <- gg_phase_plane(
#'   ode_competition,
#'   xlim       = c(0, 12),
#'   ylim       = c(0, 12),
#'   parameters = c(r1=1, r2=1, K1=10, K2=10, a12=0.5, a21=0.5)
#' )
#' result$plot +
#'   ggplot2::labs(x = "Species 1", y = "Species 2",
#'                 title = "Competition model")
#'
#' @seealso [ggphasr::gg_flow_field()], [ggphasr::gg_nullclines()],
#'   [ggphasr::gg_trajectory()], [ggphasr::find_equilibrium()],
#'   [ggphasr::classify_equilibrium()]
#' @export
gg_phase_plane <- function(deriv,
                            xlim,
                            ylim,
                            system            = c("two.dim", "one.dim"),
                            parameters        = NULL,
                            # flow field
                            n_points          = 21L,
                            arrow_type        = c("equal", "proportional"),
                            arrow_color       = "grey70",
                            # nullclines
                            show_nullclines   = TRUE,
                            nullcline_n_points = 250L,
                            # trajectories
                            y0                = NULL,
                            show_trajectories = TRUE,
                            n_ic              = 4L,
                            t_end             = 10,
                            t_start_back      = NULL,
                            trajectory_color  = "grey30",
                            # equilibria
                            find_equilibria   = TRUE,
                            eq_n_grid         = 10L,
                            eq_grid_y0        = NULL,
                            legend_position   = "right",
                            # labels
                            xlab              = NULL,
                            ylab              = "y",
                            title             = NULL) {

  # ── Input validation ─────────────────────────────────────────────────────
  system     <- match.arg(system)
  arrow_type <- match.arg(arrow_type)

  if (!is.numeric(xlim) || length(xlim) != 2L || xlim[[1L]] >= xlim[[2L]]) {
    rlang::abort("`xlim` must be a numeric vector of length 2 with xlim[1] < xlim[2].")
  }
  if (!is.numeric(ylim) || length(ylim) != 2L || ylim[[1L]] >= ylim[[2L]]) {
    rlang::abort("`ylim` must be a numeric vector of length 2 with ylim[1] < ylim[2].")
  }

  if (is.null(xlab)) xlab <- if (system == "one.dim") "t" else "x"
  n_ic <- as.integer(n_ic)

  # ── Normalize ODE once (reused throughout) ───────────────────────────────
  norm <- .normalize_ode(deriv, system = system)
  .validate_ode(norm, system = system, parameters = parameters)

  # ── 1. Flow field (base plot) ─────────────────────────────────────────────
  p <- gg_flow_field(
    deriv       = deriv,
    xlim        = xlim,
    ylim        = ylim,
    system      = system,
    parameters  = parameters,
    n_points    = n_points,
    arrow_type  = arrow_type,
    arrow_color = arrow_color,
    xlab        = xlab,
    ylab        = ylab,
    title       = title
  )

  # ── 2. Nullclines ─────────────────────────────────────────────────────────
  if (show_nullclines) {
    p <- p + gg_nullclines(
      deriv      = deriv,
      xlim       = xlim,
      ylim       = ylim,
      system     = system,
      parameters = parameters,
      n_points   = nullcline_n_points
    )
  }

  # ── 3. Trajectories ───────────────────────────────────────────────────────
  if (show_trajectories) {

    ic_list <- if (is.null(y0)) {
      .default_ics(xlim, ylim, n_ic, system)
    } else {
      .parse_initial_conditions(y0, system)
    }

    if (length(ic_list) > 0L) {
      p <- p + gg_trajectory(
        deriv           = deriv,
        y0              = ic_list,
        xlim            = xlim,
        ylim            = ylim,
        system          = system,
        parameters      = parameters,
        t_end           = t_end,
        t_start_back    = t_start_back,
        color           = trajectory_color,
        add_arrows      = TRUE,
        add_start_point = FALSE
      )
    }
  }

  # ── 4. Equilibria ─────────────────────────────────────────────────────────
  eq_table <- NULL

  if (find_equilibria) {

    # Build starting points for root-finder
    search_y0 <- if (!is.null(eq_grid_y0)) {
      eq_grid_y0
    } else {
      NULL   # triggers grid search inside find_equilibrium()
    }

    eq_list <- tryCatch(
      suppressWarnings(
        find_equilibrium(
          deriv      = deriv,
          y0         = search_y0,
          system     = system,
          parameters = parameters,
          xlim       = if (system == "two.dim") xlim else NULL,
          ylim       = ylim,
          n_grid     = eq_n_grid
        )
      ),
      error = function(e) list()
    )

    if (length(eq_list) > 0L) {

      # Classify each equilibrium
      eq_table <- do.call(rbind, lapply(eq_list, function(eq) {
        tryCatch(
          classify_equilibrium(deriv, equilibrium = eq,
                                system = system, parameters = parameters),
          error = function(e) NULL
        )
      }))
      eq_table <- eq_table[!vapply(eq_list,
                                    is.null, logical(1L)), , drop = FALSE]

      if (!is.null(eq_table) && nrow(eq_table) > 0L) {

        # Coordinates for plotting
        if (system == "one.dim") {
          eq_plot <- data.frame(
            x              = mean(xlim),   # place on phase line x position
            y              = eq_table$y,
            classification = eq_table$classification
          )
        } else {
          eq_plot <- data.frame(
            x              = eq_table$x,
            y              = eq_table$y,
            classification = eq_table$classification
          )
        }

        # Filter to equilibria within the plot bounds
        in_bounds <- eq_plot$x >= xlim[[1L]] & eq_plot$x <= xlim[[2L]] &
          eq_plot$y >= ylim[[1L]] & eq_plot$y <= ylim[[2L]]
        eq_plot <- eq_plot[in_bounds, , drop = FALSE]

        if (nrow(eq_plot) > 0L) {

          all_shapes <- .equilibrium_shapes()
          all_fills  <- .equilibrium_fills()
          classes    <- unique(eq_plot$classification)

          used_shapes <- all_shapes[classes]
          used_fills  <- all_fills[classes]

          p <- p +
            ggplot2::geom_point(
              data    = eq_plot,
              mapping = ggplot2::aes(x     = .data$x,
                                     y     = .data$y,
                                     shape = .data$classification,
                                     fill  = .data$classification),
              color  = "black",
              size   = 4,
              stroke = 1
            ) +
            ggplot2::scale_shape_manual(
              name   = "Equilibrium type",
              values = used_shapes
            ) +
            ggplot2::scale_fill_manual(
              name   = "Equilibrium type",
              values = used_fills
            )


        }
      }
    }
  }

  # ── 5. Legend position ────────────────────────────────────────────────────
  # Applied after all layers are added so it overrides any per-layer legend
  # settings from gg_nullclines(), gg_trajectory(), etc.
  p <- p + .legend_theme(legend_position)

  # ── 6. 1D phase portrait overlay ─────────────────────────────────────────
  if (system == "one.dim") {
    p <- p + gg_phase_portrait(
      deriv      = deriv,
      ylim       = ylim,
      xlim       = xlim,
      parameters = parameters
    )
  }

  # ── Return ────────────────────────────────────────────────────────────────
  structure(
    list(
      plot       = p,
      equilibria = eq_table
    ),
    class = "ggphasr_result"
  )
}


# ---------------------------------------------------------------------------
# S3 methods for ggphasr_result
# ---------------------------------------------------------------------------

#' Print method for gg_phase_plane() output
#'
#' Prints the plot and, if present, a summary of the classified equilibria.
#'
#' @param x A `ggphasr_result` object returned by [ggphasr::gg_phase_plane()].
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.ggphasr_result <- function(x, ...) {
  print(x$plot)
  if (!is.null(x$equilibria) && nrow(x$equilibria) > 0L) {
    cat("\nEquilibria:\n")
    cols <- intersect(c("x", "y", "classification", "tr", "det"),
                      names(x$equilibria))
    print(x$equilibria[, cols, drop = FALSE], digits = 4, row.names = FALSE)
  }
  invisible(x)
}

#' Add ggplot2 layers to a ggphasr_result object
#'
#' Adds a ggplot2 layer, scale, or theme element to the `$plot` component
#' of a [ggphasr::gg_phase_plane()] result, returning an updated
#' `ggphasr_result` object. This is the recommended way to further
#' customize a phase plane plot while preserving the `$equilibria` table:
#'
#' ```r
#' result <- gg_phase_plane(ode_lotka_volterra, ...)
#'
#' # Use add_layer() to customize and keep the ggphasr_result structure:
#' result2 <- add_layer(result, ggplot2::labs(title = "My plot"))
#' result2$plot        # updated plot
#' result2$equilibria  # equilibria preserved
#'
#' # Or extract $plot and add layers with + directly:
#' result$plot + ggplot2::labs(title = "My plot")
#' ```
#'
#' @param result A `ggphasr_result` object from [ggphasr::gg_phase_plane()].
#' @param layer A ggplot2 layer, scale, theme, or list thereof.
#' @return A `ggphasr_result` object with the updated `$plot`.
#' @export
add_layer <- function(result, layer) {
  if (!inherits(result, "ggphasr_result")) {
    rlang::abort("`result` must be a `ggphasr_result` object from gg_phase_plane().")
  }
  result$plot <- result$plot + layer
  result
}
