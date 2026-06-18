# ode_systems_1d.R
#
# Built-in one-dimensional ODE systems for ggphasr.
#
# All functions use Convention A (deSolve-compatible):
#   f(t, y, parameters) -> list(c(dy/dt))
#
# The `parameters` argument is always a named numeric vector so that users
# can pass parameters by name rather than relying on positional indexing.
# Default parameter values are provided so each function works out-of-the-box
# for quick exploration without specifying parameters.
#
# Functions in this file:
#   ode_exponential()      — exponential growth / decay
#   ode_logistic()         — logistic growth (Verhulst)
#   ode_monomolecular()    — monomolecular (saturating) growth
#   ode_von_bertalanffy()  — von Bertalanffy growth


# ---------------------------------------------------------------------------
# ode_exponential()
# ---------------------------------------------------------------------------

#' Exponential growth model
#'
#' The ODE for exponential (Malthusian) population growth:
#'
#' \deqn{\frac{dy}{dt} = r \, y}
#'
#' where \eqn{y(t)} is the population size (or any quantity growing
#' proportionally to itself) and \eqn{r} is the intrinsic growth rate.
#' When \eqn{r > 0} the population grows without bound; when \eqn{r < 0}
#' it decays to zero; when \eqn{r = 0} it is constant.
#'
#' The analytic solution is \eqn{y(t) = y_0 \, e^{rt}}.
#'
#' @param t Numeric scalar. Time (not used directly, included for
#'   deSolve compatibility).
#' @param y Numeric vector of length 1. Current state: \eqn{y[1]} is the
#'   population size.
#' @param parameters Named numeric vector with element:
#'   \describe{
#'     \item{`r`}{Intrinsic growth rate. Default: `0.5`.}
#'   }
#'
#' @return A list with one element: a numeric vector of length 1 containing
#'   \eqn{dy/dt}, as required by [deSolve::ode()].
#'
#' @examples
#' # Evaluate at y = 2 with default parameters (r = 0.5)
#' # dy/dt = 0.5 * 2 = 1
#' ode_exponential(t = 0, y = c(2), parameters = c(r = 0.5))
#'
#' # Use with gg_flow_field() for a 1D phase portrait
#' \dontrun{
#' gg_phase_portrait(ode_exponential, ylim = c(-3, 3))
#' }
#'
#' @seealso [ode_logistic()], [ode_monomolecular()]
#' @export
ode_exponential <- function(t, y, parameters = c(r = 0.5)) {
  r <- .get_param(parameters, "r", 1L)
  list(c(r * y[[1L]]))
}


# ---------------------------------------------------------------------------
# ode_logistic()
# ---------------------------------------------------------------------------

#' Logistic growth model
#'
#' The ODE for logistic (Verhulst) population growth:
#'
#' \deqn{\frac{dy}{dt} = r \, y \left(1 - \frac{y}{K}\right)}
#'
#' where \eqn{y(t)} is the population size, \eqn{r > 0} is the intrinsic
#' growth rate, and \eqn{K > 0} is the carrying capacity. The model has
#' two equilibria: an unstable equilibrium at \eqn{y = 0} and a stable
#' equilibrium at \eqn{y = K}.
#'
#' The analytic solution is
#' \eqn{y(t) = K / (1 + ((K - y_0)/y_0) \, e^{-rt})}.
#'
#' @param t Numeric scalar. Time (not used directly, included for
#'   deSolve compatibility).
#' @param y Numeric vector of length 1. Current state: \eqn{y[1]} is the
#'   population size.
#' @param parameters Named numeric vector with elements:
#'   \describe{
#'     \item{`r`}{Intrinsic growth rate. Default: `1`.}
#'     \item{`K`}{Carrying capacity. Default: `10`.}
#'   }
#'
#' @return A list with one element: a numeric vector of length 1 containing
#'   \eqn{dy/dt}, as required by [deSolve::ode()].
#'
#' @examples
#' # Evaluate at y = 5 with default parameters (r = 1, K = 10)
#' # dy/dt = 1 * 5 * (1 - 5/10) = 2.5
#' ode_logistic(t = 0, y = c(5), parameters = c(r = 1, K = 10))
#'
#' \dontrun{
#' gg_phase_portrait(ode_logistic, ylim = c(-2, 12))
#' }
#'
#' @seealso [ode_exponential()], [ode_monomolecular()]
#' @export
ode_logistic <- function(t, y, parameters = c(r = 1, K = 10)) {
  r <- .get_param(parameters, "r", 1L)
  K <- .get_param(parameters, "K", 2L)
  list(c(r * y[[1L]] * (1 - y[[1L]] / K)))
}


# ---------------------------------------------------------------------------
# ode_monomolecular()
# ---------------------------------------------------------------------------

#' Monomolecular growth model
#'
#' The ODE for monomolecular (saturating, or Mitscherlich) growth:
#'
#' \deqn{\frac{dy}{dt} = r \left(K - y\right)}
#'
#' where \eqn{y(t)} is the quantity of interest (e.g., biomass), \eqn{r > 0}
#' is the rate constant, and \eqn{K > 0} is the asymptote (maximum attainable
#' value). Growth decelerates monotonically as \eqn{y} approaches \eqn{K}.
#' There is a single stable equilibrium at \eqn{y = K}.
#'
#' Unlike logistic growth, there is no inflection point — growth rate is
#' greatest at \eqn{y = 0} and decreases linearly to zero at \eqn{y = K}.
#'
#' The analytic solution is \eqn{y(t) = K(1 - e^{-rt}) + y_0 \, e^{-rt}}.
#'
#' @param t Numeric scalar. Time (not used directly, included for
#'   deSolve compatibility).
#' @param y Numeric vector of length 1. Current state: \eqn{y[1]} is the
#'   current value.
#' @param parameters Named numeric vector with elements:
#'   \describe{
#'     \item{`r`}{Rate constant. Default: `1`.}
#'     \item{`K`}{Asymptote (maximum value). Default: `10`.}
#'   }
#'
#' @return A list with one element: a numeric vector of length 1 containing
#'   \eqn{dy/dt}, as required by [deSolve::ode()].
#'
#' @examples
#' # Evaluate at y = 4 with default parameters (r = 1, K = 10)
#' # dy/dt = 1 * (10 - 4) = 6
#' ode_monomolecular(t = 0, y = c(4), parameters = c(r = 1, K = 10))
#'
#' \dontrun{
#' gg_phase_portrait(ode_monomolecular, ylim = c(0, 12))
#' }
#'
#' @seealso [ode_logistic()], [ode_von_bertalanffy()]
#' @export
ode_monomolecular <- function(t, y, parameters = c(r = 1, K = 10)) {
  r <- .get_param(parameters, "r", 1L)
  K <- .get_param(parameters, "K", 2L)
  list(c(r * (K - y[[1L]])))
}


# ---------------------------------------------------------------------------
# ode_von_bertalanffy()
# ---------------------------------------------------------------------------

#' Von Bertalanffy growth model
#'
#' The ODE for von Bertalanffy growth, widely used in fisheries biology and
#' ecology to model growth of organisms:
#'
#' \deqn{\frac{dy}{dt} = \alpha \, y^{2/3} - \beta \, y}
#'
#' where \eqn{y(t)} is body mass (or a proportional measure of body size),
#' \eqn{\alpha > 0} is the anabolism coefficient (growth rate proportional to
#' surface area, hence the \eqn{2/3} power), and \eqn{\beta > 0} is the
#' catabolism coefficient (decay proportional to mass).
#'
#' The model has a single non-trivial stable equilibrium at
#' \eqn{y^* = (\alpha / \beta)^3}, which corresponds to the asymptotic
#' body size \eqn{W_\infty} in fisheries notation.
#'
#' @param t Numeric scalar. Time (not used directly, included for
#'   deSolve compatibility).
#' @param y Numeric vector of length 1. Current state: \eqn{y[1]} is body
#'   mass (must be non-negative).
#' @param parameters Named numeric vector with elements:
#'   \describe{
#'     \item{`alpha`}{Anabolism coefficient. Default: `1`.}
#'     \item{`beta`}{Catabolism coefficient. Default: `0.5`.}
#'   }
#'
#' @return A list with one element: a numeric vector of length 1 containing
#'   \eqn{dy/dt}, as required by [deSolve::ode()].
#'
#' @details
#' Note that \eqn{y^{2/3}} is not defined for \eqn{y < 0}. Although the
#' biological interpretation restricts \eqn{y \geq 0}, no guard is applied
#' internally — users should restrict phase portrait axes accordingly.
#'
#' @examples
#' # Evaluate at y = 8 with default parameters (alpha = 1, beta = 0.5)
#' # dy/dt = 1 * 8^(2/3) - 0.5 * 8 = 4 - 4 = 0  (equilibrium)
#' ode_von_bertalanffy(t = 0, y = c(8), parameters = c(alpha = 1, beta = 0.5))
#'
#' \dontrun{
#' gg_phase_portrait(ode_von_bertalanffy, ylim = c(0, 12))
#' }
#'
#' @seealso [ode_monomolecular()], [ode_logistic()]
#' @export
ode_von_bertalanffy <- function(t, y, parameters = c(alpha = 1, beta = 0.5)) {
  alpha <- .get_param(parameters, "alpha", 1L)
  beta  <- .get_param(parameters, "beta",  2L)
  list(c(alpha * y[[1L]]^(2/3) - beta * y[[1L]]))
}
