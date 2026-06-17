# ode_systems_2d.R
#
# Built-in two-dimensional autonomous ODE systems for ggphasr.
#
# All functions use Convention A (deSolve-compatible):
#   f(t, y, parameters) -> list(c(dy1/dt, dy2/dt))
#
# Parameter extraction uses the internal helper .get_param() from
# ode_utils.R, which handles both named and positional parameter vectors.
# Default parameter values are provided so each function works immediately
# for classroom exploration without specifying parameters explicitly.
#
# Functions in this file:
#   ode_lotka_volterra()   — predator-prey (Lotka-Volterra)
#   ode_sir()              — SIR epidemic (proportions or counts)
#   ode_van_der_pol()      — Van der Pol oscillator
#   ode_simple_pendulum()  — simple pendulum (with/without damping)
#   ode_competition()      — two-species Lotka-Volterra competition
#   ode_toggle()           — genetic toggle switch
#   ode_morris_lecar()     — Morris-Lecar neuron model
#   ode_lindemann()        — Lindemann mechanism (chemical kinetics)


# ---------------------------------------------------------------------------
# ode_lotka_volterra()
# ---------------------------------------------------------------------------

#' Lotka-Volterra predator-prey model
#'
#' The classic two-species predator-prey system:
#'
#' \deqn{\frac{dx}{dt} = \alpha x - \beta x y}
#' \deqn{\frac{dy}{dt} = \delta x y - \gamma y}
#'
#' where \eqn{x(t)} is prey abundance, \eqn{y(t)} is predator abundance,
#' \eqn{\alpha > 0} is the prey growth rate, \eqn{\beta > 0} is the
#' predation rate, \eqn{\delta > 0} is the predator growth rate per prey
#' consumed, and \eqn{\gamma > 0} is the predator death rate.
#'
#' The system has two equilibria: a trivial unstable equilibrium at
#' \eqn{(0, 0)} and a neutrally stable center at
#' \eqn{(x^*, y^*) = (\gamma/\delta,\; \alpha/\beta)}.
#' All trajectories starting in the positive quadrant are closed orbits
#' around the interior equilibrium.
#'
#' @param t Numeric scalar. Time (autonomous system; included for deSolve
#'   compatibility).
#' @param y Numeric vector of length 2. State vector: \eqn{y[1]} = prey
#'   (\eqn{x}), \eqn{y[2]} = predator.
#' @param parameters Named numeric vector with elements:
#'   \describe{
#'     \item{`alpha`}{Prey growth rate. Default: `1`.}
#'     \item{`beta`}{Predation rate. Default: `0.5`.}
#'     \item{`delta`}{Predator growth rate per prey consumed. Default: `0.5`.}
#'     \item{`gamma`}{Predator death rate. Default: `1`.}
#'   }
#'
#' @return A list with one element: a numeric vector of length 2 containing
#'   \eqn{(dx/dt, dy/dt)}, as required by [deSolve::ode()].
#'
#' @examples
#' # Interior equilibrium at (gamma/delta, alpha/beta) = (2, 2)
#' ode_lotka_volterra(t = 0, y = c(2, 2),
#'                   parameters = c(alpha = 1, beta = 0.5,
#'                                  delta = 0.5, gamma = 1))
#'
#' \dontrun{
#' gg_flow_field(ode_lotka_volterra, xlim = c(0, 5), ylim = c(0, 5))
#' }
#'
#' @seealso [ggphasr::ode_competition()]
#' @export
ode_lotka_volterra <- function(t, y,
                                parameters = c(alpha = 1, beta  = 0.5,
                                               delta = 0.5, gamma = 1)) {
  alpha <- .get_param(parameters, "alpha", 1L)
  beta  <- .get_param(parameters, "beta",  2L)
  delta <- .get_param(parameters, "delta", 3L)
  gamma <- .get_param(parameters, "gamma", 4L)
  x <- y[[1L]]; v <- y[[2L]]
  list(c(
    alpha * x - beta  * x * v,
    delta * x * v - gamma * v
  ))
}


# ---------------------------------------------------------------------------
# ode_sir()
# ---------------------------------------------------------------------------

#' SIR epidemic model
#'
#' The classic Kermack-McKendrick SIR (Susceptible-Infected-Recovered)
#' epidemic model. Supports two formulations via the `scale` parameter:
#'
#' **Proportions** (`scale = "proportions"`, default):
#' \deqn{\frac{dS}{dt} = -\beta S I}
#' \deqn{\frac{dI}{dt} = \beta S I - \gamma I}
#'
#' where \eqn{S + I + R = 1} at all times. Here \eqn{\beta} is the
#' transmission rate and \eqn{\gamma} is the recovery rate.
#' The basic reproduction number is \eqn{R_0 = \beta / \gamma}.
#'
#' **Counts** (`scale = "counts"`):
#' \deqn{\frac{dS}{dt} = -\beta S I / N}
#' \deqn{\frac{dI}{dt} = \beta S I / N - \gamma I}
#'
#' where \eqn{N} is the total (constant) population size. The
#' proportions formulation is recovered by setting \eqn{N = 1}.
#'
#' Note: only \eqn{S} and \eqn{I} are tracked as state variables.
#' \eqn{R} can be recovered as \eqn{R = N - S - I}.
#'
#' @param t Numeric scalar. Time (autonomous system; included for deSolve
#'   compatibility).
#' @param y Numeric vector of length 2. State vector: \eqn{y[1]} = \eqn{S}
#'   (susceptible), \eqn{y[2]} = \eqn{I} (infected).
#' @param parameters Named numeric vector with elements:
#'   \describe{
#'     \item{`beta`}{Transmission rate. Default: `0.5`.}
#'     \item{`gamma`}{Recovery rate. Default: `0.1`.}
#'     \item{`N`}{Total population size (used only when
#'       `scale = "counts"`). Default: `1000`.}
#'     \item{`scale`}{Character: `"proportions"` (default) or `"counts"`.
#'       Can be passed as a named element of the `parameters` list, or as
#'       a separate argument.}
#'   }
#'
#' @return A list with one element: a numeric vector of length 2 containing
#'   \eqn{(dS/dt, dI/dt)}, as required by [deSolve::ode()].
#'
#' @details
#' The `scale` argument is passed as part of `parameters` to maintain
#' deSolve compatibility. Because mixing character and numeric values in
#' `c()` coerces everything to character, pass `parameters` as a `list()`
#' when including `scale`:
#' `parameters = list(beta = 0.5, gamma = 0.1, N = 1000, scale = "counts")`.
#' Passing a plain numeric `c()` vector (without `scale`) also works and is
#' the recommended form for the proportions formulation.
#' When `scale` is absent, `"proportions"` is assumed.
#'
#' @examples
#' # Proportions formulation — plain numeric vector is fine
#' ode_sir(t = 0, y = c(0.99, 0.01),
#'         parameters = c(beta = 0.5, gamma = 0.1))
#'
#' # Counts formulation — use list() to avoid character coercion
#' ode_sir(t = 0, y = c(990, 10),
#'         parameters = list(beta = 0.5, gamma = 0.1, N = 1000,
#'                           scale = "counts"))
#'
#' \dontrun{
#' gg_flow_field(ode_sir, xlim = c(0, 1), ylim = c(0, 1))
#' }
#'
#' @export
ode_sir <- function(t, y,
                    parameters = c(beta = 0.5, gamma = 0.1)) {
  # Use as.numeric() on all numeric parameters: when `scale` is included in
  # a c() vector, R coerces the entire vector to character, so explicit
  # conversion is required here.
  beta  <- as.numeric(.get_param(parameters, "beta",  1L))
  gamma <- as.numeric(.get_param(parameters, "gamma", 2L))

  # Scale: read from named parameters list; default to "proportions"
  scale <- if (!is.null(names(parameters)) &&
               "scale" %in% names(parameters)) {
    parameters[["scale"]]
  } else {
    "proportions"
  }

  N <- if (!is.null(names(parameters)) &&
           "N" %in% names(parameters)) {
    as.numeric(parameters[["N"]])
  } else {
    1
  }

  # For proportions, N = 1 (recovers the standard formulation)
  if (scale == "proportions") N <- 1

  S <- y[[1L]]; I <- y[[2L]]
  list(c(
    -beta  * S * I / N,
     beta  * S * I / N - gamma * I
  ))
}


# ---------------------------------------------------------------------------
# ode_van_der_pol()
# ---------------------------------------------------------------------------

#' Van der Pol oscillator
#'
#' The Van der Pol oscillator, a nonlinear oscillator with self-sustaining
#' oscillations (a limit cycle), written as a 2D first-order system:
#'
#' \deqn{\frac{dx}{dt} = y}
#' \deqn{\frac{dy}{dt} = \mu (1 - x^2) y - x}
#'
#' where \eqn{x(t)} is the position (displacement), \eqn{y(t)} is the
#' velocity, and \eqn{\mu \geq 0} is the nonlinearity/damping parameter.
#'
#' When \eqn{\mu = 0} the system reduces to a harmonic oscillator with
#' circular orbits. For \eqn{\mu > 0} a stable limit cycle exists; the
#' shape becomes increasingly relaxation-oscillator-like as \eqn{\mu}
#' increases. The origin \eqn{(0, 0)} is an unstable spiral for
#' \eqn{\mu > 0}.
#'
#' @param t Numeric scalar. Time (autonomous; included for deSolve
#'   compatibility).
#' @param y Numeric vector of length 2. State vector: \eqn{y[1]} = position
#'   \eqn{x}, \eqn{y[2]} = velocity.
#' @param parameters Named numeric vector with element:
#'   \describe{
#'     \item{`mu`}{Nonlinearity parameter. Default: `1`. Must be \eqn{\geq 0}.}
#'   }
#'
#' @return A list with one element: a numeric vector of length 2 containing
#'   \eqn{(dx/dt, dy/dt)}, as required by [deSolve::ode()].
#'
#' @examples
#' # At (x=2, y=0) with mu=1: dx/dt=0, dy/dt = 1*(1-4)*0 - 2 = -2
#' ode_van_der_pol(t = 0, y = c(2, 0), parameters = c(mu = 1))
#'
#' \dontrun{
#' gg_flow_field(ode_van_der_pol, xlim = c(-4, 4), ylim = c(-4, 4))
#' }
#'
#' @export
ode_van_der_pol <- function(t, y, parameters = c(mu = 1)) {
  mu <- .get_param(parameters, "mu", 1L)
  x  <- y[[1L]]; v <- y[[2L]]
  list(c(
    v,
    mu * (1 - x^2) * v - x
  ))
}


# ---------------------------------------------------------------------------
# ode_simple_pendulum()
# ---------------------------------------------------------------------------

#' Simple pendulum model
#'
#' The simple pendulum written as a 2D first-order system:
#'
#' \deqn{\frac{d\theta}{dt} = \omega}
#' \deqn{\frac{d\omega}{dt} = -\frac{g}{L} \sin(\theta) - b\,\omega}
#'
#' where \eqn{\theta} is the angle from vertical (radians), \eqn{\omega}
#' is the angular velocity, \eqn{g/L} is the ratio of gravitational
#' acceleration to pendulum length (combined into a single parameter
#' \code{gL} for convenience), and \eqn{b \geq 0} is a damping coefficient.
#'
#' Setting \eqn{b = 0} gives the undamped (conservative) pendulum with
#' heteroclinic orbits connecting the unstable equilibria at
#' \eqn{\theta = \pm\pi}. Adding damping (\eqn{b > 0}) makes the stable
#' equilibrium at \eqn{\theta = 0} a stable spiral (for small \eqn{b}) or
#' stable node (for large \eqn{b}).
#'
#' @param t Numeric scalar. Time (autonomous; included for deSolve
#'   compatibility).
#' @param y Numeric vector of length 2. State vector: \eqn{y[1]} = angle
#'   \eqn{\theta} (radians), \eqn{y[2]} = angular velocity \eqn{\omega}.
#' @param parameters Named numeric vector with elements:
#'   \describe{
#'     \item{`gL`}{Ratio \eqn{g/L} (gravitational acceleration divided by
#'       pendulum length). Default: `1`.}
#'     \item{`b`}{Damping coefficient. Default: `0` (undamped).}
#'   }
#'
#' @return A list with one element: a numeric vector of length 2 containing
#'   \eqn{(d\theta/dt, d\omega/dt)}, as required by [deSolve::ode()].
#'
#' @examples
#' # Undamped pendulum at theta = pi/4, omega = 0
#' ode_simple_pendulum(t = 0, y = c(pi/4, 0), parameters = c(gL = 1, b = 0))
#'
#' # Damped pendulum
#' ode_simple_pendulum(t = 0, y = c(pi/4, 0), parameters = c(gL = 1, b = 0.5))
#'
#' \dontrun{
#' gg_flow_field(ode_simple_pendulum,
#'               xlim = c(-pi, pi), ylim = c(-3, 3))
#' }
#'
#' @export
ode_simple_pendulum <- function(t, y, parameters = c(gL = 1, b = 0)) {
  gL <- .get_param(parameters, "gL", 1L)
  b  <- .get_param(parameters, "b",  2L)
  theta <- y[[1L]]; omega <- y[[2L]]
  list(c(
    omega,
    -gL * sin(theta) - b * omega
  ))
}


# ---------------------------------------------------------------------------
# ode_competition()
# ---------------------------------------------------------------------------

#' Two-species Lotka-Volterra competition model
#'
#' The classic two-species interspecific competition model:
#'
#' \deqn{\frac{dN_1}{dt} = r_1 N_1 \left(1 - \frac{N_1 + \alpha_{12} N_2}{K_1}\right)}
#' \deqn{\frac{dN_2}{dt} = r_2 N_2 \left(1 - \frac{N_2 + \alpha_{21} N_1}{K_2}\right)}
#'
#' where \eqn{N_1, N_2} are species abundances, \eqn{r_1, r_2} are
#' intrinsic growth rates, \eqn{K_1, K_2} are carrying capacities, and
#' \eqn{\alpha_{12}, \alpha_{21}} are interspecific competition coefficients
#' (the effect of species 2 on species 1, and vice versa).
#'
#' The four possible outcomes (competitive exclusion of species 1 or 2,
#' stable coexistence, or an unstable equilibrium with priority effects)
#' depend on the relative magnitudes of \eqn{K_1}, \eqn{K_2},
#' \eqn{\alpha_{12}}, and \eqn{\alpha_{21}}, making this an excellent
#' teaching example for phase plane analysis.
#'
#' @param t Numeric scalar. Time (autonomous; included for deSolve
#'   compatibility).
#' @param y Numeric vector of length 2. State vector: \eqn{y[1]} = \eqn{N_1},
#'   \eqn{y[2]} = \eqn{N_2}.
#' @param parameters Named numeric vector with elements:
#'   \describe{
#'     \item{`r1`}{Growth rate of species 1. Default: `1`.}
#'     \item{`r2`}{Growth rate of species 2. Default: `1`.}
#'     \item{`K1`}{Carrying capacity of species 1. Default: `10`.}
#'     \item{`K2`}{Carrying capacity of species 2. Default: `10`.}
#'     \item{`a12`}{Effect of species 2 on species 1 (\eqn{\alpha_{12}}).
#'       Default: `0.5`.}
#'     \item{`a21`}{Effect of species 1 on species 2 (\eqn{\alpha_{21}}).
#'       Default: `0.5`.}
#'   }
#'
#' @return A list with one element: a numeric vector of length 2 containing
#'   \eqn{(dN_1/dt, dN_2/dt)}, as required by [deSolve::ode()].
#'
#' @examples
#' ode_competition(t = 0, y = c(5, 5),
#'                parameters = c(r1 = 1, r2 = 1,
#'                               K1 = 10, K2 = 10,
#'                               a12 = 0.5, a21 = 0.5))
#'
#' \dontrun{
#' gg_flow_field(ode_competition, xlim = c(0, 15), ylim = c(0, 15))
#' }
#'
#' @seealso [ggphasr::ode_lotka_volterra()]
#' @export
ode_competition <- function(t, y,
                             parameters = c(r1  = 1,   r2  = 1,
                                            K1  = 10,  K2  = 10,
                                            a12 = 0.5, a21 = 0.5)) {
  r1  <- .get_param(parameters, "r1",  1L)
  r2  <- .get_param(parameters, "r2",  2L)
  K1  <- .get_param(parameters, "K1",  3L)
  K2  <- .get_param(parameters, "K2",  4L)
  a12 <- .get_param(parameters, "a12", 5L)
  a21 <- .get_param(parameters, "a21", 6L)
  N1 <- y[[1L]]; N2 <- y[[2L]]
  list(c(
    r1 * N1 * (1 - (N1 + a12 * N2) / K1),
    r2 * N2 * (1 - (N2 + a21 * N1) / K2)
  ))
}


# ---------------------------------------------------------------------------
# ode_toggle()
# ---------------------------------------------------------------------------

#' Genetic toggle switch model
#'
#' The Gardner et al. (2000) mutual repressor model of a synthetic genetic
#' toggle switch:
#'
#' \deqn{\frac{du}{dt} = \frac{\alpha_1}{1 + v^\beta} - u}
#' \deqn{\frac{dv}{dt} = \frac{\alpha_2}{1 + u^\gamma} - v}
#'
#' where \eqn{u} and \eqn{v} are the concentrations of two mutually
#' repressing proteins. \eqn{\alpha_1, \alpha_2 > 0} are the effective
#' synthesis rates (incorporating promoter strength and degradation), and
#' \eqn{\beta, \gamma > 0} are the cooperativity (Hill) coefficients.
#'
#' For sufficiently large \eqn{\alpha_1} and \eqn{\alpha_2} with
#' cooperativity \eqn{> 1}, the system is bistable: two stable equilibria
#' (each corresponding to one protein being dominant) separated by an
#' unstable saddle point. This bistability is the basis of a synthetic
#' biological memory device.
#'
#' @param t Numeric scalar. Time (autonomous; included for deSolve
#'   compatibility).
#' @param y Numeric vector of length 2. State vector: \eqn{y[1]} = \eqn{u},
#'   \eqn{y[2]} = \eqn{v}.
#' @param parameters Named numeric vector with elements:
#'   \describe{
#'     \item{`alpha1`}{Effective synthesis rate of \eqn{u}. Default: `3`.}
#'     \item{`alpha2`}{Effective synthesis rate of \eqn{v}. Default: `3`.}
#'     \item{`beta`}{Cooperativity coefficient in \eqn{u} equation.
#'       Default: `2`.}
#'     \item{`gamma`}{Cooperativity coefficient in \eqn{v} equation.
#'       Default: `2`.}
#'   }
#'
#' @return A list with one element: a numeric vector of length 2 containing
#'   \eqn{(du/dt, dv/dt)}, as required by [deSolve::ode()].
#'
#' @references
#' Gardner TS, Cantor CR, Collins JJ (2000). Construction of a genetic
#' toggle switch in Escherichia coli. *Nature* 403: 339–342.
#' \doi{10.1038/35002131}
#'
#' @examples
#' ode_toggle(t = 0, y = c(2, 0.5),
#'            parameters = c(alpha1 = 3, alpha2 = 3, beta = 2, gamma = 2))
#'
#' \dontrun{
#' gg_flow_field(ode_toggle, xlim = c(0, 4), ylim = c(0, 4))
#' }
#'
#' @export
ode_toggle <- function(t, y,
                        parameters = c(alpha1 = 3, alpha2 = 3,
                                       beta   = 2, gamma  = 2)) {
  alpha1 <- .get_param(parameters, "alpha1", 1L)
  alpha2 <- .get_param(parameters, "alpha2", 2L)
  beta   <- .get_param(parameters, "beta",   3L)
  gamma  <- .get_param(parameters, "gamma",  4L)
  u <- y[[1L]]; v <- y[[2L]]
  list(c(
    alpha1 / (1 + v^beta)  - u,
    alpha2 / (1 + u^gamma) - v
  ))
}


# ---------------------------------------------------------------------------
# ode_morris_lecar()
# ---------------------------------------------------------------------------

#' Morris-Lecar neuron model
#'
#' The Morris-Lecar (1981) conductance-based neuron model, a 2D reduction of
#' Hodgkin-Huxley dynamics for a barnacle muscle fiber:
#'
#' \deqn{C \frac{dV}{dt} = I - g_{Ca} M_\infty(V)(V - V_{Ca}) - g_K N (V - V_K) - g_L (V - V_L)}
#' \deqn{\frac{dN}{dt} = \phi \frac{N_\infty(V) - N}{\tau_N(V)}}
#'
#' where the voltage-dependent steady-state functions are:
#' \deqn{M_\infty(V) = \frac{1}{2}\left(1 + \tanh\!\left(\frac{V - V_1}{V_2}\right)\right)}
#' \deqn{N_\infty(V) = \frac{1}{2}\left(1 + \tanh\!\left(\frac{V - V_3}{V_4}\right)\right)}
#' \deqn{\tau_N(V) = \left(\cosh\!\left(\frac{V - V_3}{2 V_4}\right)\right)^{-1}}
#'
#' The state variables are membrane potential \eqn{V} (mV) and \eqn{N},
#' the probability that a K\eqn{^+} channel is open. Depending on the
#' applied current \eqn{I} and other parameters, the system can exhibit
#' a stable resting state, a stable limit cycle (repetitive firing), or
#' bistability, making it a rich teaching example for bifurcation analysis.
#'
#' @param t Numeric scalar. Time (autonomous; included for deSolve
#'   compatibility).
#' @param y Numeric vector of length 2. State vector: \eqn{y[1]} = membrane
#'   potential \eqn{V} (mV), \eqn{y[2]} = K\eqn{^+} channel open probability
#'   \eqn{N} (dimensionless, \eqn{0 \leq N \leq 1}).
#' @param parameters Named numeric vector with elements:
#'   \describe{
#'     \item{`I`}{Applied current (\eqn{\mu}A/cm\eqn{^2}). Default: `0.`}
#'     \item{`C`}{Membrane capacitance (\eqn{\mu}F/cm\eqn{^2}). Default: `20`.}
#'     \item{`gCa`}{Maximum Ca\eqn{^{2+}} conductance (mS/cm\eqn{^2}).
#'       Default: `4.4`.}
#'     \item{`gK`}{Maximum K\eqn{^+} conductance (mS/cm\eqn{^2}).
#'       Default: `8`.}
#'     \item{`gL`}{Leak conductance (mS/cm\eqn{^2}). Default: `2`.}
#'     \item{`VCa`}{Ca\eqn{^{2+}} reversal potential (mV). Default: `120`.}
#'     \item{`VK`}{K\eqn{^+} reversal potential (mV). Default: `-84`.}
#'     \item{`VL`}{Leak reversal potential (mV). Default: `-60`.}
#'     \item{`V1`}{Voltage at half-activation of Ca\eqn{^{2+}} (mV).
#'       Default: `-1.2`.}
#'     \item{`V2`}{Slope of Ca\eqn{^{2+}} activation (mV). Default: `18`.}
#'     \item{`V3`}{Voltage at half-activation of K\eqn{^+} (mV).
#'       Default: `2`.}
#'     \item{`V4`}{Slope of K\eqn{^+} activation (mV). Default: `30`.}
#'     \item{`phi`}{Reference frequency (dimensionless). Default: `0.04`.}
#'   }
#'
#' @return A list with one element: a numeric vector of length 2 containing
#'   \eqn{(dV/dt, dN/dt)}, as required by [deSolve::ode()].
#'
#' @references
#' Morris C, Lecar H (1981). Voltage oscillations in the barnacle giant
#' muscle fiber. *Biophysical Journal* 35(1): 193-213.
#' \doi{10.1016/S0006-3495(81)84782-0}
#'
#' @examples
#' # At resting potential with zero applied current
#' ode_morris_lecar(t = 0, y = c(-60, 0),
#'                 parameters = c(I = 0, C = 20,
#'                                gCa = 4.4, gK = 8, gL = 2,
#'                                VCa = 120, VK = -84, VL = -60,
#'                                V1 = -1.2, V2 = 18, V3 = 2, V4 = 30,
#'                                phi = 0.04))
#'
#' \dontrun{
#' gg_flow_field(ode_morris_lecar,
#'               xlim = c(-80, 60), ylim = c(0, 0.6))
#' }
#'
#' @export
ode_morris_lecar <- function(t, y,
                              parameters = c(I   =  0,    C   = 20,
                                             gCa =  4.4,  gK  =  8,   gL  = 2,
                                             VCa =  120,  VK  = -84,  VL  = -60,
                                             V1  = -1.2,  V2  = 18,
                                             V3  =  2,    V4  = 30,
                                             phi =  0.04)) {
  I   <- .get_param(parameters, "I",   1L)
  C   <- .get_param(parameters, "C",   2L)
  gCa <- .get_param(parameters, "gCa", 3L)
  gK  <- .get_param(parameters, "gK",  4L)
  gL  <- .get_param(parameters, "gL",  5L)
  VCa <- .get_param(parameters, "VCa", 6L)
  VK  <- .get_param(parameters, "VK",  7L)
  VL  <- .get_param(parameters, "VL",  8L)
  V1  <- .get_param(parameters, "V1",  9L)
  V2  <- .get_param(parameters, "V2",  10L)
  V3  <- .get_param(parameters, "V3",  11L)
  V4  <- .get_param(parameters, "V4",  12L)
  phi <- .get_param(parameters, "phi", 13L)

  V <- y[[1L]]; N <- y[[2L]]

  # Voltage-gated steady-state functions
  M_inf   <- 0.5 * (1 + tanh((V - V1) / V2))
  N_inf   <- 0.5 * (1 + tanh((V - V3) / V4))
  tau_N   <- 1 / cosh((V - V3) / (2 * V4))

  # Ionic currents
  I_Ca <- gCa * M_inf * (V - VCa)
  I_K  <- gK  * N     * (V - VK)
  I_L  <- gL           * (V - VL)

  list(c(
    (I - I_Ca - I_K - I_L) / C,
    phi * (N_inf - N) / tau_N
  ))
}


# ---------------------------------------------------------------------------
# ode_lindemann()
# ---------------------------------------------------------------------------

#' Lindemann mechanism (chemical kinetics)
#'
#' The Lindemann-Christiansen mechanism for unimolecular gas-phase reactions,
#' tracking the concentrations of reactant \eqn{A} and activated
#' intermediate \eqn{A^*}:
#'
#' Reaction scheme:
#' \deqn{A + A \xrightarrow{k_1} A^* + A \quad \text{(activation)}}
#' \deqn{A^* + A \xrightarrow{k_{-1}} A + A \quad \text{(deactivation)}}
#' \deqn{A^* \xrightarrow{k_2} P \quad \text{(product formation)}}
#'
#' The resulting ODE system for \eqn{[A]} and \eqn{[A^*]} is:
#'
#' \deqn{\frac{d[A]}{dt}   = -k_1 [A]^2 + k_{-1} [A^*][A]}
#' \deqn{\frac{d[A^*]}{dt} =  k_1 [A]^2 - k_{-1} [A^*][A] - k_2 [A^*]}
#'
#' @param t Numeric scalar. Time (autonomous; included for deSolve
#'   compatibility).
#' @param y Numeric vector of length 2. State vector: \eqn{y[1]} = \eqn{[A]},
#'   \eqn{y[2]} = \eqn{[A^*]}.
#' @param parameters Named numeric vector with elements:
#'   \describe{
#'     \item{`k1`}{Activation rate constant. Default: `1`.}
#'     \item{`k_1`}{Deactivation rate constant. Default: `1`.}
#'     \item{`k2`}{Product formation rate constant. Default: `0.5`.}
#'   }
#'
#' @return A list with one element: a numeric vector of length 2 containing
#'   \eqn{(d[A]/dt, d[A^*]/dt)}, as required by [deSolve::ode()].
#'
#' @examples
#' ode_lindemann(t = 0, y = c(2, 0.1),
#'              parameters = c(k1 = 1, k_1 = 1, k2 = 0.5))
#'
#' \dontrun{
#' gg_flow_field(ode_lindemann, xlim = c(0, 3), ylim = c(0, 1))
#' }
#'
#' @export
ode_lindemann <- function(t, y,
                           parameters = c(k1 = 1, k_1 = 1, k2 = 0.5)) {
  k1  <- .get_param(parameters, "k1",  1L)
  k_1 <- .get_param(parameters, "k_1", 2L)
  k2  <- .get_param(parameters, "k2",  3L)
  A  <- y[[1L]]; As <- y[[2L]]   # As = [A*]
  list(c(
    -k1 * A^2 + k_1 * As * A,
     k1 * A^2 - k_1 * As * A - k2 * As
  ))
}
