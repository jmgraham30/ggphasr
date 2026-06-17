# ode_utils.R
#
# Internal utilities for ODE function handling.
#
# Every user-facing function in ggphasr accepts a `deriv` argument — a
# function describing the ODE system. Two calling conventions are supported:
#
#   Convention A (deSolve-compatible, primary):
#     f(t, y, parameters) -> list(c(dy1, dy2, ...))
#     - First argument is time `t` (even for autonomous systems)
#     - Second argument is the state vector `y`
#     - Third argument is a parameters vector or list
#     - Return value is a list whose first element is the derivatives vector
#     - This is the convention used by phaseR and deSolve::ode()
#
#   Convention B (simplified, secondary):
#     2D: f(x, y, parameters = NULL) -> c(dx, dy)
#     1D: f(y, parameters = NULL)    -> dy (scalar)
#     - No time argument (autonomous systems only)
#     - Returns a plain numeric vector, not a list
#     - Simpler to write for users new to deSolve
#
# Convention B functions are wrapped into Convention A by .normalize_ode()
# before being passed to deSolve or any internal computation. This means
# all downstream code in ggphasr works exclusively with Convention A, and
# the dual-convention support is fully contained in this file.
#
# Key exported (internal) functions:
#   .detect_ode_convention()  — returns "A" or "B"
#   .normalize_ode()          — wraps Convention B into Convention A
#   .validate_ode()           — checks that a function is a valid ODE for ggphasr
#   .eval_ode()               — evaluates a (normalized) ODE at a point


# ---------------------------------------------------------------------------
# .detect_ode_convention()
# ---------------------------------------------------------------------------
#
# Inspects the formal argument names of `deriv` to determine which calling
# convention is being used.
#
# Detection rules (applied in order):
#   1. If the first argument is named "t"  -> Convention A
#   2. If the first argument is named "y" and there is no "x" argument
#      (i.e., 1D system in Convention B)   -> Convention B (1D)
#   3. If the first two arguments are "x" and "y"
#      (i.e., 2D system in Convention B)   -> Convention B (2D)
#   4. Otherwise -> error: unrecognized convention
#
# @param deriv A function representing the ODE system.
# @return A character string: "A", "B1" (Convention B, 1D), or "B2"
#   (Convention B, 2D).
#
.detect_ode_convention <- function(deriv) {

  if (!is.function(deriv)) {
    rlang::abort(
      "`deriv` must be a function.",
      call = rlang::caller_env()
    )
  }

  arg_names <- names(formals(deriv))

  if (length(arg_names) == 0L) {
    rlang::abort(
      paste0(
        "`deriv` must have at least one argument. ",
        "See `?ggphasr` for supported ODE conventions."
      ),
      call = rlang::caller_env()
    )
  }

  first_arg <- arg_names[1L]

  # Convention A: deSolve-compatible — first argument is time "t"
  if (first_arg == "t") {
    return("A")
  }

  # Convention B (1D): first argument is "y", no "x" argument present
  if (first_arg == "y" && !("x" %in% arg_names)) {
    return("B1")
  }

  # Convention B (2D): first two arguments are "x" and "y"
  if (length(arg_names) >= 2L && first_arg == "x" && arg_names[2L] == "y") {
    return("B2")
  }

  # Unrecognized
  rlang::abort(
    paste0(
      "Could not detect the ODE convention used by `deriv`. ",
      "Expected first argument to be `t` (Convention A), ",
      "`y` (Convention B, 1D), or `x` with second argument `y` ",
      "(Convention B, 2D). Got: `",
      paste(arg_names, collapse = "`, `"),
      "`."
    ),
    call = rlang::caller_env()
  )
}


# ---------------------------------------------------------------------------
# .normalize_ode()
# ---------------------------------------------------------------------------
#
# Wraps a Convention B function into Convention A so that all downstream
# code can work with a single interface. Convention A functions are returned
# unchanged.
#
# @param deriv   A function in Convention A, B1, or B2.
# @param system  Character string: "one.dim" or "two.dim". Used to validate
#                that the detected convention is consistent with the system
#                dimension being plotted.
# @return A function with signature f(t, y, parameters) -> list(c(...)),
#   suitable for passing to deSolve::ode() or .eval_ode().
#
.normalize_ode <- function(deriv, system = c("two.dim", "one.dim")) {

  system     <- match.arg(system)
  convention <- .detect_ode_convention(deriv)

  # Validate convention against system dimension
  if (system == "one.dim" && convention == "B2") {
    rlang::abort(
      paste0(
        "`deriv` appears to be a 2D Convention B function (arguments `x`, `y`) ",
        "but `system = \"one.dim\"` was specified. ",
        "For 1D systems, use f(y, parameters) returning a scalar."
      ),
      call = rlang::caller_env()
    )
  }

  if (system == "two.dim" && convention == "B1") {
    rlang::abort(
      paste0(
        "`deriv` appears to be a 1D Convention B function (argument `y` only) ",
        "but `system = \"two.dim\"` was specified. ",
        "For 2D systems, use f(x, y, parameters) returning c(dx, dy)."
      ),
      call = rlang::caller_env()
    )
  }

  # Convention A: pass through unchanged
  if (convention == "A") {
    return(deriv)
  }

  # Convention B1 (1D): wrap f(y, parameters) -> list(dy)
  if (convention == "B1") {
    normalized <- function(t, y, parameters) {
      list(c(deriv(y[1L], parameters)))
    }
    return(normalized)
  }

  # Convention B2 (2D): wrap f(x, y, parameters) -> list(c(dx, dy))
  if (convention == "B2") {
    normalized <- function(t, y, parameters) {
      list(c(deriv(y[1L], y[2L], parameters)))
    }
    return(normalized)
  }
}


# ---------------------------------------------------------------------------
# .validate_ode()
# ---------------------------------------------------------------------------
#
# Performs a test evaluation of `deriv` at a known point to catch common
# errors early (wrong return type, wrong output length, errors thrown by the
# function itself) before any grid computation or deSolve call.
#
# Called internally by gg_flow_field(), gg_trajectory(), etc., at the top
# of each function so that error messages reference the user-facing call.
#
# @param deriv       A *normalized* (Convention A) ODE function.
# @param system      "one.dim" or "two.dim".
# @param parameters  Parameter vector/list to use in the test call.
# @return Invisibly returns TRUE if validation passes; otherwise aborts with
#   an informative error message.
#
.validate_ode <- function(deriv,
                          system     = c("two.dim", "one.dim"),
                          parameters = NULL) {

  system   <- match.arg(system)
  test_y   <- if (system == "one.dim") c(1.0) else c(1.0, 1.0)
  expected <- if (system == "one.dim") 1L     else 2L

  result <- tryCatch(
    deriv(t = 0, y = test_y, parameters = parameters),
    error = function(e) {
      rlang::abort(
        paste0(
          "`deriv` threw an error during a test evaluation:\n",
          conditionMessage(e)
        ),
        call = rlang::caller_env()
      )
    }
  )

  # Must return a list
  if (!is.list(result)) {
    rlang::abort(
      paste0(
        "`deriv` must return a list (e.g., `list(c(dy1, dy2))`), ",
        "but returned an object of class `", class(result)[1L], "`."
      ),
      call = rlang::caller_env()
    )
  }

  # First element must be a numeric vector of the right length
  derivs <- result[[1L]]

  if (!is.numeric(derivs)) {
    rlang::abort(
      paste0(
        "The first element of the list returned by `deriv` must be numeric, ",
        "but got class `", class(derivs)[1L], "`."
      ),
      call = rlang::caller_env()
    )
  }

  if (length(derivs) != expected) {
    rlang::abort(
      paste0(
        "For a `system = \"", system, "\"` ODE, `deriv` must return a list ",
        "whose first element has length ", expected, ", ",
        "but got length ", length(derivs), "."
      ),
      call = rlang::caller_env()
    )
  }

  invisible(TRUE)
}


# ---------------------------------------------------------------------------
# .eval_ode()
# ---------------------------------------------------------------------------
#
# Evaluates a normalized (Convention A) ODE at a single point (t, y) and
# returns the derivatives as a plain numeric vector (not a list). Used by
# the grid-computation functions inside gg_flow_field() and gg_nullclines().
#
# @param deriv       A normalized (Convention A) ODE function.
# @param t           Scalar time value (irrelevant for autonomous systems,
#                    but required by the deSolve interface).
# @param y           Numeric vector of state values, length 1 (1D) or 2 (2D).
# @param parameters  Parameter vector/list passed through to `deriv`.
# @return A plain numeric vector of derivatives.
#
.eval_ode <- function(deriv, t = 0, y, parameters = NULL) {
  deriv(t = t, y = y, parameters = parameters)[[1L]]
}
