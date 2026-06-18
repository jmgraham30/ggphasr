# find_equilibrium.R
# classify_equilibrium.R
#
# Numerical equilibrium finding and classification for ggphasr.
#
# find_equilibrium():
#   Locates equilibria of a 1D or 2D ODE system using rootSolve::multiroot()
#   (Newton-Raphson). Accepts either a single initial guess or a grid of
#   starting points to search for multiple equilibria simultaneously.
#
# classify_equilibrium():
#   Classifies a known equilibrium point using the Jacobian matrix computed
#   via finite differences. Returns a tidy data frame row with the
#   equilibrium location, Jacobian, eigenvalues, and classification string.
#
# Internal helpers:
#   .numerical_jacobian()   — finite-difference Jacobian at a point
#   .classify_2d()          — trace-determinant classification for 2D systems
#   .classify_1d()          — derivative-sign classification for 1D systems
#   .deduplicate_equilibria() — removes duplicate roots from grid searches


# ---------------------------------------------------------------------------
# .numerical_jacobian()
# ---------------------------------------------------------------------------
#
# Computes the Jacobian matrix of a normalized ODE at a point y using
# forward finite differences.
#
# For an n-dimensional system f(y), the Jacobian J has entries:
#   J[i,j] = df_i/dy_j ≈ (f_i(y + h*e_j) - f_i(y)) / h
#
# where e_j is the j-th standard basis vector and h is a small step size.
#
# @param deriv       Normalized (Convention A) ODE function.
# @param y           Numeric vector: the point at which to evaluate J.
# @param parameters  Parameter vector/list.
# @param h           Step size for finite differences. Default: 1e-6.
# @return A numeric matrix of dimensions length(y) x length(y).
#
.numerical_jacobian <- function(deriv, y, parameters, h = 1e-6) {
  n  <- length(y)
  f0 <- .eval_ode(deriv, t = 0, y = y, parameters = parameters)
  J  <- matrix(0, nrow = n, ncol = n)

  for (j in seq_len(n)) {
    y_perturb      <- y
    y_perturb[[j]] <- y_perturb[[j]] + h
    fj             <- .eval_ode(deriv, t = 0, y = y_perturb,
                                parameters = parameters)
    J[, j]         <- (fj - f0) / h
  }

  J
}


# ---------------------------------------------------------------------------
# .classify_2d()
# ---------------------------------------------------------------------------
#
# Classifies a 2D equilibrium using the trace-determinant plane method.
#
# Given eigenvalues lambda1, lambda2 of the Jacobian at the equilibrium:
#   tr  = lambda1 + lambda2  (trace of J)
#   det = lambda1 * lambda2  (determinant of J)
#   disc = tr^2 - 4*det      (discriminant)
#
# Classification:
#   det < 0                          -> Saddle
#   det > 0, disc > 0, tr < 0       -> Stable node
#   det > 0, disc > 0, tr > 0       -> Unstable node
#   det > 0, disc < 0, tr < 0       -> Stable spiral
#   det > 0, disc < 0, tr > 0       -> Unstable spiral
#   det > 0, disc < 0, tr = 0       -> Center (neutrally stable)
#   det = 0                         -> Non-isolated equilibrium
#
# @param eigenvalues Complex vector of length 2.
# @param tol         Tolerance for treating values as zero. Default: 1e-8.
# @return Character string classification.
#
.classify_2d <- function(eigenvalues, tol = 1e-8) {
  tr   <- sum(Re(eigenvalues))
  det  <- prod(Re(eigenvalues)) + prod(Im(eigenvalues))
  disc <- tr^2 - 4 * det
  # Recompute det directly from eigenvalue product
  det  <- Re(eigenvalues[[1L]]) * Re(eigenvalues[[2L]]) -
    Im(eigenvalues[[1L]]) * Im(eigenvalues[[2L]])

  if (abs(det) < tol) return("Non-isolated equilibrium")

  if (det < -tol) return("Saddle")

  # det > 0
  both_real <- all(abs(Im(eigenvalues)) < tol)

  if (both_real) {
    if (tr < -tol) return("Stable node")
    if (tr >  tol) return("Unstable node")
    return("Non-isolated equilibrium")
  } else {
    # Complex eigenvalues (spiral or center)
    if (tr < -tol) return("Stable spiral")
    if (tr >  tol) return("Unstable spiral")
    return("Center")
  }
}


# ---------------------------------------------------------------------------
# .classify_1d()
# ---------------------------------------------------------------------------
#
# Classifies a 1D equilibrium from the sign of df/dy at the equilibrium.
# df/dy < 0 -> stable; df/dy > 0 -> unstable; df/dy = 0 -> inconclusive.
#
# @param jacobian_val Scalar: df/dy at the equilibrium (J[1,1]).
# @param tol          Tolerance. Default: 1e-8.
# @return Character string: "Stable", "Unstable", or "Inconclusive (df/dy = 0)".
#
.classify_1d <- function(jacobian_val, tol = 1e-8) {
  if (jacobian_val < -tol) return("Stable")
  if (jacobian_val >  tol) return("Unstable")
  return("Inconclusive (df/dy = 0)")
}


# ---------------------------------------------------------------------------
# .deduplicate_equilibria()
# ---------------------------------------------------------------------------
#
# Removes duplicate equilibrium points from a list of roots found by
# grid-search multi-start. Two equilibria are considered duplicates if
# they are within `tol` of each other in Euclidean distance.
#
# @param roots List of numeric vectors (equilibrium locations).
# @param tol   Distance tolerance. Default: 1e-4.
# @return A deduplicated list of numeric vectors.
#
.deduplicate_equilibria <- function(roots, tol = 1e-4) {
  if (length(roots) <= 1L) return(roots)

  keep <- rep(TRUE, length(roots))
  for (i in seq_len(length(roots) - 1L)) {
    if (!keep[[i]]) next
    for (j in seq(i + 1L, length(roots))) {
      if (!keep[[j]]) next
      dist <- sqrt(sum((roots[[i]] - roots[[j]])^2))
      if (dist < tol) keep[[j]] <- FALSE
    }
  }

  roots[keep]
}


# ---------------------------------------------------------------------------
# find_equilibrium()
# ---------------------------------------------------------------------------

#' Find equilibria of an ODE system numerically
#'
#' Locates one or more equilibrium points of a one- or two-dimensional
#' autonomous ODE system using Newton-Raphson root-finding via
#' [rootSolve::multiroot()]. Accepts either a single initial guess or a
#' grid of starting points to search for multiple equilibria.
#'
#' @param deriv A function describing the ODE system, in Convention A or B.
#'   See [ggphasr] for details.
#' @param y0 Initial guess(es) for the root-finder. One of:
#'   \itemize{
#'     \item A numeric vector of length 1 (1D) or 2 (2D): a single starting
#'       point.
#'     \item A numeric matrix with one row per starting point.
#'     \item A list of numeric vectors.
#'     \item `NULL`: triggers automatic grid search over `xlim` × `ylim`
#'       (requires `xlim` and `ylim` to be supplied).
#'   }
#' @param system Character: `"two.dim"` (default) or `"one.dim"`.
#' @param parameters Parameter vector or list passed to `deriv`.
#' @param xlim Numeric vector of length 2. x-axis search range. Required
#'   when `y0 = NULL`. Default: `NULL`.
#' @param ylim Numeric vector of length 2. y-axis search range. Required
#'   when `y0 = NULL`. Default: `NULL`.
#' @param n_grid Integer. Grid resolution per axis for automatic search
#'   when `y0 = NULL`. Default: `10` (giving up to 100 starting points
#'   for 2D systems).
#' @param tol Numeric. Convergence tolerance passed to
#'   [rootSolve::multiroot()]. Default: `1e-8`.
#' @param dedup_tol Numeric. Distance tolerance for removing duplicate
#'   roots. Two equilibria closer than this are treated as the same point.
#'   Default: `1e-4`.
#'
#' @return A list of numeric vectors, each giving the coordinates of one
#'   equilibrium point. If only one equilibrium is found, a list of length
#'   1 is returned (not an unwrapped vector) for consistency.
#'
#' @details
#' Newton-Raphson is sensitive to the initial guess — it may fail to
#' converge or converge to a non-equilibrium point for badly chosen
#' starting values. When `y0 = NULL`, the function tries `n_grid^2` starting
#' points on a regular grid, collects all converged roots, and deduplicates
#' them. This is more reliable than a single guess but still not exhaustive.
#'
#' Failed convergences are silently discarded; only roots where the
#' residual is below `tol * 100` are retained.
#'
#' @examples
#' # Single guess: find the interior equilibrium of Lotka-Volterra
#' find_equilibrium(
#'   ode_lotka_volterra,
#'   y0         = c(1.5, 1.5),
#'   parameters = c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1)
#' )
#'
#' # Automatic grid search: find all equilibria of example 11
#' find_equilibrium(
#'   ode_example_11,
#'   y0    = NULL,
#'   xlim  = c(0, 4),
#'   ylim  = c(0, 4)
#' )
#'
#' # 1D system: find equilibria of the logistic equation
#' find_equilibrium(
#'   ode_logistic,
#'   y0         = NULL,
#'   system     = "one.dim",
#'   ylim       = c(-1, 12),
#'   parameters = c(r = 1, K = 10)
#' )
#'
#' @seealso [ggphasr::classify_equilibrium()]
#' @export
find_equilibrium <- function(deriv,
                              y0         = NULL,
                              system     = c("two.dim", "one.dim"),
                              parameters = NULL,
                              xlim       = NULL,
                              ylim       = NULL,
                              n_grid     = 10L,
                              tol        = 1e-8,
                              dedup_tol  = 1e-4) {

  # ── Input validation ─────────────────────────────────────────────────────
  system <- match.arg(system)
  n_grid <- as.integer(n_grid)

  if (is.null(y0) && is.null(ylim)) {
    rlang::abort(
      "When `y0 = NULL`, `ylim` must be supplied to define the search range."
    )
  }
  if (is.null(y0) && system == "two.dim" && is.null(xlim)) {
    rlang::abort(
      "When `y0 = NULL` and `system = \"two.dim\"`, `xlim` must be supplied."
    )
  }

  # ── Normalize ODE ────────────────────────────────────────────────────────
  norm <- .normalize_ode(deriv, system = system)
  .validate_ode(norm, system = system, parameters = parameters)

  # Wrapper for rootSolve: takes a vector y, returns derivatives as vector
  f_root <- function(y) {
    .eval_ode(norm, t = 0, y = y, parameters = parameters)
  }

  # ── Build list of starting points ────────────────────────────────────────
  if (is.null(y0)) {
    # Automatic grid search
    if (system == "one.dim") {
      grid_y  <- seq(ylim[[1L]], ylim[[2L]], length.out = n_grid)
      starts  <- lapply(grid_y, function(yi) c(yi))
    } else {
      grid_x  <- seq(xlim[[1L]], xlim[[2L]], length.out = n_grid)
      grid_y  <- seq(ylim[[1L]], ylim[[2L]], length.out = n_grid)
      grid    <- expand.grid(x = grid_x, y = grid_y)
      starts  <- lapply(seq_len(nrow(grid)),
                         function(i) c(grid$x[[i]], grid$y[[i]]))
    }
  } else {
    starts <- .parse_initial_conditions(y0, system)
  }

  # ── Run root-finder from each starting point ──────────────────────────────
  roots <- list()

  for (start in starts) {
    result <- tryCatch(
      rootSolve::multiroot(f = f_root, start = start, ctol = tol,
                            rtol = tol, atol = tol),
      error   = function(e) NULL,
      warning = function(w) NULL
    )

    if (is.null(result)) next

    # Check that the residual is small (genuine root, not failed convergence)
    residual <- tryCatch(
      max(abs(f_root(result$root))),
      error = function(e) Inf
    )

    if (is.finite(residual) && residual < tol * 100) {
      roots[[length(roots) + 1L]] <- result$root
    }
  }

  if (length(roots) == 0L) {
    rlang::warn(
      "No equilibria found. Try different starting points or a wider search range."
    )
    return(list())
  }

  # ── Deduplicate ───────────────────────────────────────────────────────────
  roots <- .deduplicate_equilibria(roots, tol = dedup_tol)

  # Sort for consistent output: by x first, then y
  if (length(roots) > 1L) {
    order_idx <- order(vapply(roots, function(r) r[[1L]], numeric(1L)))
    roots     <- roots[order_idx]
  }

  roots
}


# ---------------------------------------------------------------------------
# classify_equilibrium()
# ---------------------------------------------------------------------------

#' Classify an equilibrium point of an ODE system
#'
#' Computes the Jacobian matrix of a one- or two-dimensional ODE system at
#' a known equilibrium point using forward finite differences, then
#' classifies the equilibrium using eigenvalue analysis. Returns a tidy
#' data frame row that can be combined with results from multiple equilibria
#' using [rbind()].
#'
#' @param deriv A function describing the ODE system, in Convention A or B.
#' @param equilibrium Numeric vector of length 1 (1D) or 2 (2D) giving the
#'   equilibrium coordinates. Typically obtained from [find_equilibrium()].
#' @param system Character: `"two.dim"` (default) or `"one.dim"`.
#' @param parameters Parameter vector or list passed to `deriv`.
#' @param h Numeric. Step size for finite-difference Jacobian computation.
#'   Default: `1e-6`.
#' @param tol Numeric. Tolerance for treating eigenvalue real/imaginary
#'   parts as zero. Default: `1e-8`.
#'
#' @return A [data.frame()] with one row and the following columns:
#'   \describe{
#'     \item{`x`}{x-coordinate of the equilibrium (2D) or `NA` (1D).}
#'     \item{`y`}{y-coordinate of the equilibrium (both 1D and 2D).}
#'     \item{`classification`}{Character. One of: `"Stable node"`,
#'       `"Unstable node"`, `"Stable spiral"`, `"Unstable spiral"`,
#'       `"Center"`, `"Saddle"`, `"Non-isolated equilibrium"` (2D);
#'       or `"Stable"`, `"Unstable"`, `"Inconclusive (df/dy = 0)"` (1D).}
#'     \item{`tr`}{Trace of the Jacobian (2D only; `NA` for 1D).}
#'     \item{`det`}{Determinant of the Jacobian (2D only; `NA` for 1D).}
#'     \item{`jacobian_11`}{J[1,1] — always present.}
#'     \item{`jacobian_12`}{J[1,2] (2D only; `NA` for 1D).}
#'     \item{`jacobian_21`}{J[2,1] (2D only; `NA` for 1D).}
#'     \item{`jacobian_22`}{J[2,2] (2D only; `NA` for 1D).}
#'     \item{`lambda_1_re`}{Real part of eigenvalue 1 (2D only; `NA` for 1D).}
#'     \item{`lambda_1_im`}{Imaginary part of eigenvalue 1 (2D only).}
#'     \item{`lambda_2_re`}{Real part of eigenvalue 2 (2D only).}
#'     \item{`lambda_2_im`}{Imaginary part of eigenvalue 2 (2D only).}
#'   }
#'
#' @details
#' Results from multiple equilibria can be combined into a summary table:
#'
#' ```r
#' eq_list <- find_equilibrium(ode_example_11, y0 = NULL,
#'                              xlim = c(0,4), ylim = c(0,4))
#' results <- do.call(rbind, lapply(eq_list, function(eq) {
#'   classify_equilibrium(ode_example_11, equilibrium = eq)
#' }))
#' ```
#'
#' @examples
#' # Classify the interior equilibrium of Lotka-Volterra (expect: Center)
#' classify_equilibrium(
#'   ode_lotka_volterra,
#'   equilibrium = c(2, 2),
#'   parameters  = c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1)
#' )
#'
#' # Classify all equilibria of example 11
#' eq_list <- find_equilibrium(ode_example_11, y0 = NULL,
#'                              xlim = c(0, 4), ylim = c(0, 4))
#' do.call(rbind, lapply(eq_list, function(eq) {
#'   classify_equilibrium(ode_example_11, equilibrium = eq)
#' }))
#'
#' @seealso [ggphasr::find_equilibrium()]
#' @export
classify_equilibrium <- function(deriv,
                                  equilibrium,
                                  system     = c("two.dim", "one.dim"),
                                  parameters = NULL,
                                  h          = 1e-6,
                                  tol        = 1e-8) {

  # ── Input validation ─────────────────────────────────────────────────────
  system      <- match.arg(system)
  expected    <- if (system == "one.dim") 1L else 2L
  equilibrium <- as.numeric(equilibrium)

  if (length(equilibrium) != expected) {
    rlang::abort(
      paste0(
        "For `system = \"", system, "\"`, `equilibrium` must have length ",
        expected, ". Got length ", length(equilibrium), "."
      )
    )
  }

  # ── Normalize ODE ────────────────────────────────────────────────────────
  norm <- .normalize_ode(deriv, system = system)

  # ── Compute Jacobian ─────────────────────────────────────────────────────
  J <- .numerical_jacobian(norm, equilibrium, parameters, h = h)

  # ── Classify ─────────────────────────────────────────────────────────────
  if (system == "one.dim") {

    classification <- .classify_1d(J[[1L, 1L]], tol = tol)

    result <- data.frame(
      x              = NA_real_,
      y              = equilibrium[[1L]],
      classification = classification,
      tr             = NA_real_,
      det            = NA_real_,
      jacobian_11    = J[[1L, 1L]],
      jacobian_12    = NA_real_,
      jacobian_21    = NA_real_,
      jacobian_22    = NA_real_,
      lambda_1_re    = NA_real_,
      lambda_1_im    = NA_real_,
      lambda_2_re    = NA_real_,
      lambda_2_im    = NA_real_,
      stringsAsFactors = FALSE
    )

  } else {

    eig            <- eigen(J, only.values = TRUE)$values
    classification <- .classify_2d(eig, tol = tol)
    tr_val         <- sum(diag(J))
    det_val        <- det(J)

    result <- data.frame(
      x              = equilibrium[[1L]],
      y              = equilibrium[[2L]],
      classification = classification,
      tr             = tr_val,
      det            = det_val,
      jacobian_11    = J[[1L, 1L]],
      jacobian_12    = J[[1L, 2L]],
      jacobian_21    = J[[2L, 1L]],
      jacobian_22    = J[[2L, 2L]],
      lambda_1_re    = Re(eig[[1L]]),
      lambda_1_im    = Im(eig[[1L]]),
      lambda_2_re    = Re(eig[[2L]]),
      lambda_2_im    = Im(eig[[2L]]),
      stringsAsFactors = FALSE
    )
  }

  result
}
