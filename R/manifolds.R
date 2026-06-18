# manifolds.R
#
# gg_manifolds() — draw stable and/or unstable manifolds of a saddle point
# for a 2D autonomous ODE system.
#
# Algorithm:
#   1. Obtain the Jacobian at the saddle (from a pre-computed
#      classify_equilibrium() result or by computing it internally).
#   2. Compute eigenvalues and eigenvectors of the Jacobian.
#   3. Identify the stable eigenvector (negative eigenvalue) and
#      unstable eigenvector (positive eigenvalue).
#   4. Place four seed points: +epsilon and -epsilon along each eigenvector.
#   5. Integrate forward from the two unstable seeds -> unstable manifold.
#   6. Integrate backward from the two stable seeds  -> stable manifold.
#   7. Return the resulting trajectory paths as ggplot2 layers.
#
# Internal helpers:
#   .get_manifold_eigenvectors()  — extracts stable/unstable eigenvectors
#   .integrate_manifold_branch()  — integrates one manifold branch


# ---------------------------------------------------------------------------
# .get_manifold_eigenvectors()
# ---------------------------------------------------------------------------

#' @keywords internal
.get_manifold_eigenvectors <- function(deriv,
                                        equilibrium,
                                        parameters,
                                        eq_classified = NULL,
                                        h             = 1e-6) {

  if (!is.null(eq_classified)) {
    J <- matrix(
      c(eq_classified$jacobian_11, eq_classified$jacobian_21,
        eq_classified$jacobian_12, eq_classified$jacobian_22),
      nrow = 2L, ncol = 2L
    )
  } else {
    J <- .numerical_jacobian(deriv, equilibrium, parameters, h = h)
  }

  eig  <- eigen(J)
  vals <- eig$values
  vecs <- eig$vectors

  re_vals <- Re(vals)
  idx_s   <- which.min(re_vals)
  idx_u   <- which.max(re_vals)

  lambda_s <- vals[[idx_s]]
  lambda_u <- vals[[idx_u]]

  vec_s <- Re(vecs[, idx_s])
  vec_u <- Re(vecs[, idx_u])
  vec_s <- vec_s / sqrt(sum(vec_s^2))
  vec_u <- vec_u / sqrt(sum(vec_u^2))

  list(
    stable_vec   = vec_s,
    unstable_vec = vec_u,
    lambda_s     = lambda_s,
    lambda_u     = lambda_u,
    J            = J
  )
}


# ---------------------------------------------------------------------------
# .integrate_manifold_branch()
# ---------------------------------------------------------------------------

#' @keywords internal
.integrate_manifold_branch <- function(deriv,
                                        seed,
                                        t_end,
                                        t_steps,
                                        parameters,
                                        method = "lsoda") {

  result <- .integrate_trajectory(
    deriv      = deriv,
    y0         = seed,
    t_start    = 0,
    t_end      = t_end,
    t_steps    = t_steps,
    parameters = parameters,
    method     = method
  )

  if (is.null(result)) return(NULL)
  result[, c("x", "y")]
}


# ---------------------------------------------------------------------------
# .build_manifold_data()
# ---------------------------------------------------------------------------
#
# Integrates all branches of one manifold (stable or unstable) and returns
# a single tidy data frame with atomic character columns branch and
# manifold_type. Using a for-loop with explicit index avoids the
# sapply(seeds, identical, seed) pattern which can produce list-type columns.
#
# @param norm          Normalized ODE function.
# @param seeds         List of length-2 numeric seed vectors.
# @param t_end         Numeric: integration end time (negative = backward).
# @param t_steps       Integer: number of integration steps.
# @param parameters    Parameter vector/list.
# @param method        deSolve method string.
# @param label         Character: "Stable manifold" or "Unstable manifold".
# @param prefix        Character: single-letter branch prefix ("s" or "u").
# @return A data frame with columns x, y, branch, manifold_type, or NULL.
#
.build_manifold_data <- function(norm, seeds, t_end, t_steps,
                                  parameters, method, label, prefix) {

  segs <- vector("list", length(seeds))

  for (i in seq_along(seeds)) {
    df <- .integrate_manifold_branch(norm, seeds[[i]], t_end,
                                      t_steps, parameters, method)
    if (!is.null(df)) {
      df$branch        <- paste0(prefix, i)    # plain character: "s1", "s2"
      df$manifold_type <- label                # plain character
      segs[[i]]        <- df
    }
  }

  segs <- Filter(Negate(is.null), segs)
  if (length(segs) == 0L) return(NULL)

  result <- do.call(rbind, segs)
  # Explicitly coerce to atomic types to prevent list-column issues in ggplot2
  result$branch        <- as.character(result$branch)
  result$manifold_type <- as.character(result$manifold_type)
  result$x             <- as.numeric(result$x)
  result$y             <- as.numeric(result$y)
  result
}


# ---------------------------------------------------------------------------
# gg_manifolds()
# ---------------------------------------------------------------------------

#' Draw stable and unstable manifolds of a saddle point
#'
#' Computes and draws the stable and/or unstable manifolds of a saddle
#' point for a two-dimensional autonomous ODE system. Returns a list of
#' [ggplot2] layer objects that can be added to a phase plane plot with `+`.
#'
#' The manifolds are computed by:
#' \enumerate{
#'   \item Finding the eigenvectors of the Jacobian at the saddle point.
#'   \item Placing seed points at distance `epsilon` from the saddle along
#'     each eigenvector (both + and - directions, giving four seeds total).
#'   \item Integrating **forward** from the two unstable seeds to trace the
#'     unstable manifold.
#'   \item Integrating **backward** from the two stable seeds to trace the
#'     stable manifold.
#' }
#'
#' If a [ggphasr::classify_equilibrium()] result is supplied via
#' `eq_classified`, its Jacobian entries are reused to avoid redundant
#' computation. Otherwise the Jacobian is computed internally via finite
#' differences.
#'
#' @param deriv A function describing the 2D ODE system, in Convention A
#'   or B. See [ggphasr] for details.
#' @param equilibrium Numeric vector of length 2 giving the saddle point
#'   coordinates. Typically obtained from [ggphasr::find_equilibrium()].
#' @param parameters Parameter vector or list passed to `deriv`.
#' @param eq_classified A data frame returned by
#'   [ggphasr::classify_equilibrium()], or `NULL` (default). When supplied,
#'   the Jacobian is read from this object rather than recomputed.
#' @param draw_stable Logical. Whether to draw the stable manifold.
#'   Default: `TRUE`.
#' @param draw_unstable Logical. Whether to draw the unstable manifold.
#'   Default: `TRUE`.
#' @param t_manifold Numeric. Integration time for each manifold branch.
#'   Default: `10`.
#' @param t_steps Integer. Number of integration steps per branch.
#'   Default: `500`.
#' @param epsilon Numeric. Distance from the saddle to the seed points.
#'   Default: `1e-4`.
#' @param method Character. deSolve integration method. Default: `"lsoda"`.
#' @param stable_color Character. Color of the stable manifold.
#'   Default: `"#4575b4"` (blue).
#' @param unstable_color Character. Color of the unstable manifold.
#'   Default: `"#d73027"` (red).
#' @param linewidth Numeric. Line width for both manifolds. Default: `1`.
#' @param stable_linetype Character or integer. Line type for the stable
#'   manifold. Default: `"dashed"`.
#' @param unstable_linetype Character or integer. Line type for the unstable
#'   manifold. Default: `"solid"`.
#' @param add_arrows Logical. Whether to add direction arrows. Default: `TRUE`.
#' @param arrow_size Numeric. Arrow head size in lines. Default: `0.35`.
#' @param add_legend Logical. Whether to add a legend. Default: `FALSE`.
#'   Note: setting `TRUE` adds `scale_color_manual()` and
#'   `scale_linetype_manual()` layers, which will conflict if the parent
#'   plot already has color or linetype scales.
#' @param h Numeric. Finite-difference step for Jacobian. Default: `1e-6`.
#'
#' @return A list of [ggplot2] layer objects. Add to a ggplot with `+`.
#'   Returns an empty list (invisibly) with a warning if the equilibrium
#'   is not a saddle point.
#'
#' @examples
#' # Saddle at the origin of example_08: dx/dt = x, dy/dt = -y
#' gg_flow_field(ode_example_08,
#'               xlim = c(-3, 3), ylim = c(-3, 3)) +
#'   gg_manifolds(ode_example_08, equilibrium = c(0, 0))
#'
#' # Full workflow: find, classify, then draw manifolds
#' eq    <- find_equilibrium(ode_example_11, y0 = c(0.8, 0.8))
#' eq_cl <- classify_equilibrium(ode_example_11, equilibrium = eq[[1L]])
#'
#' gg_flow_field(ode_example_11,
#'               xlim = c(0, 4), ylim = c(0, 3)) +
#'   gg_nullclines(ode_example_11, xlim = c(0, 4), ylim = c(0, 3)) +
#'   gg_manifolds(ode_example_11,
#'                equilibrium   = eq[[1L]],
#'                eq_classified = eq_cl,
#'                t_manifold    = 5)
#'
#' @seealso [ggphasr::find_equilibrium()], [ggphasr::classify_equilibrium()],
#'   [ggphasr::gg_flow_field()], [ggphasr::gg_trajectory()]
#' @export
gg_manifolds <- function(deriv,
                          equilibrium,
                          parameters        = NULL,
                          eq_classified     = NULL,
                          draw_stable       = TRUE,
                          draw_unstable     = TRUE,
                          t_manifold        = 10,
                          t_steps           = 500L,
                          epsilon           = 1e-4,
                          method            = "lsoda",
                          stable_color      = "#4575b4",
                          unstable_color    = "#d73027",
                          linewidth         = 1,
                          stable_linetype   = "dashed",
                          unstable_linetype = "solid",
                          add_arrows        = TRUE,
                          arrow_size        = 0.35,
                          add_legend        = FALSE,
                          h                 = 1e-6) {

  # ── Input validation ─────────────────────────────────────────────────────
  equilibrium <- as.numeric(equilibrium)
  if (length(equilibrium) != 2L) {
    rlang::abort("`equilibrium` must be a numeric vector of length 2.")
  }
  if (!draw_stable && !draw_unstable) {
    rlang::abort("At least one of `draw_stable` or `draw_unstable` must be TRUE.")
  }
  if (t_manifold <= 0) {
    rlang::abort("`t_manifold` must be positive.")
  }
  t_steps <- as.integer(t_steps)

  # ── Normalize ODE ────────────────────────────────────────────────────────
  norm <- .normalize_ode(deriv, system = "two.dim")
  .validate_ode(norm, system = "two.dim", parameters = parameters)

  # ── Check classification if pre-computed result supplied ─────────────────
  if (!is.null(eq_classified)) {
    if (!inherits(eq_classified, "data.frame") ||
        !"classification" %in% names(eq_classified)) {
      rlang::abort(
        "`eq_classified` must be a data frame returned by classify_equilibrium()."
      )
    }
    if (eq_classified$classification[[1L]] != "Saddle") {
      rlang::warn(
        paste0(
          "The supplied equilibrium is classified as \"",
          eq_classified$classification[[1L]],
          "\", not \"Saddle\". ",
          "Manifolds are only meaningful for saddle points. ",
          "Returning empty layer list."
        )
      )
      return(invisible(list()))
    }
  }

  # ── Get eigenvectors ──────────────────────────────────────────────────────
  eig_info <- tryCatch(
    .get_manifold_eigenvectors(norm, equilibrium, parameters,
                                eq_classified, h),
    error = function(e) {
      rlang::abort(
        paste0("Failed to compute eigenvectors at the equilibrium: ",
               conditionMessage(e))
      )
    }
  )

  # ── Check that the equilibrium is actually a saddle (no eq_classified) ───
  if (is.null(eq_classified)) {
    re_s <- Re(eig_info$lambda_s)
    re_u <- Re(eig_info$lambda_u)
    if (!(re_s < 0 && re_u > 0)) {
      rlang::warn(
        paste0(
          "The equilibrium does not appear to be a saddle point ",
          "(eigenvalues: ", round(eig_info$lambda_s, 4),
          ", ", round(eig_info$lambda_u, 4), "). ",
          "Manifolds are only meaningful for saddle points. ",
          "Returning empty layer list."
        )
      )
      return(invisible(list()))
    }
  }

  # ── Build seed points ─────────────────────────────────────────────────────
  vec_s <- eig_info$stable_vec
  vec_u <- eig_info$unstable_vec

  seeds_stable   <- list(equilibrium + epsilon * vec_s,
                          equilibrium - epsilon * vec_s)
  seeds_unstable <- list(equilibrium + epsilon * vec_u,
                          equilibrium - epsilon * vec_u)

  # ── Arrow spec ────────────────────────────────────────────────────────────
  arrow_spec <- if (add_arrows) {
    grid::arrow(length = grid::unit(arrow_size, "lines"),
                type   = "open", ends = "last")
  } else NULL

  layers <- list()

  # ── Stable manifold ───────────────────────────────────────────────────────
  if (draw_stable) {
    stable_data <- .build_manifold_data(
      norm, seeds_stable, -t_manifold, t_steps,
      parameters, method, "Stable manifold", "s"
    )

    if (!is.null(stable_data)) {
      if (add_legend) {
        layers[[length(layers) + 1L]] <- ggplot2::geom_path(
          data    = stable_data,
          mapping = ggplot2::aes(x        = .data$x,
                                 y        = .data$y,
                                 group    = .data$branch,
                                 color    = .data$manifold_type,
                                 linetype = .data$manifold_type),
          linewidth = linewidth,
          arrow     = arrow_spec,
          na.rm     = TRUE,
          lineend   = "round"
        )
      } else {
        layers[[length(layers) + 1L]] <- ggplot2::geom_path(
          data    = stable_data,
          mapping = ggplot2::aes(x     = .data$x,
                                 y     = .data$y,
                                 group = .data$branch),
          color     = stable_color,
          linetype  = stable_linetype,
          linewidth = linewidth,
          arrow     = arrow_spec,
          na.rm     = TRUE,
          lineend   = "round"
        )
      }
    }
  }

  # ── Unstable manifold ─────────────────────────────────────────────────────
  if (draw_unstable) {
    unstable_data <- .build_manifold_data(
      norm, seeds_unstable, t_manifold, t_steps,
      parameters, method, "Unstable manifold", "u"
    )

    if (!is.null(unstable_data)) {
      if (add_legend) {
        layers[[length(layers) + 1L]] <- ggplot2::geom_path(
          data    = unstable_data,
          mapping = ggplot2::aes(x        = .data$x,
                                 y        = .data$y,
                                 group    = .data$branch,
                                 color    = .data$manifold_type,
                                 linetype = .data$manifold_type),
          linewidth = linewidth,
          arrow     = arrow_spec,
          na.rm     = TRUE,
          lineend   = "round"
        )
      } else {
        layers[[length(layers) + 1L]] <- ggplot2::geom_path(
          data    = unstable_data,
          mapping = ggplot2::aes(x     = .data$x,
                                 y     = .data$y,
                                 group = .data$branch),
          color     = unstable_color,
          linetype  = unstable_linetype,
          linewidth = linewidth,
          arrow     = arrow_spec,
          na.rm     = TRUE,
          lineend   = "round"
        )
      }
    }
  }

  # ── Color and linetype scales (legend path only) ──────────────────────────
  if (add_legend && length(layers) > 0L) {
    color_vals    <- c()
    linetype_vals <- c()
    if (draw_stable) {
      color_vals[["Stable manifold"]]    <- stable_color
      linetype_vals[["Stable manifold"]] <- stable_linetype
    }
    if (draw_unstable) {
      color_vals[["Unstable manifold"]]    <- unstable_color
      linetype_vals[["Unstable manifold"]] <- unstable_linetype
    }
    layers[[length(layers) + 1L]] <- ggplot2::scale_color_manual(
      name = NULL, values = color_vals
    )
    layers[[length(layers) + 1L]] <- ggplot2::scale_linetype_manual(
      name = NULL, values = linetype_vals
    )
  }

  # ── Saddle point marker ───────────────────────────────────────────────────
  layers[[length(layers) + 1L]] <- ggplot2::geom_point(
    data    = data.frame(x = equilibrium[[1L]], y = equilibrium[[2L]]),
    mapping = ggplot2::aes(x = .data$x, y = .data$y),
    shape  = 23,
    fill   = "white",
    color  = "black",
    size   = 3,
    stroke = 1
  )

  layers
}
