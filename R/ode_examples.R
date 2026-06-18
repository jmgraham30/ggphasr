# ode_examples.R
#
# Fifteen generic textbook ODE example systems, ported from phaseR.
#
# These are the same systems as phaseR's example1() through example15(),
# renamed to ode_example_01() through ode_example_15() following the
# ggphasr naming convention (zero-padded for consistent sorting).
#
# Dimension key:
#   Examples 01-05  one-dimensional (1D)
#   Examples 06-15  two-dimensional (2D)
#
# All functions use Convention A (deSolve-compatible):
#   1D: f(t, y, parameters) -> list(c(dy/dt))
#   2D: f(t, y, parameters) -> list(c(dx/dt, dy/dt))
#
# None of the example systems require parameters; the `parameters` argument
# is accepted for deSolve compatibility but ignored.
#
# Source: Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis
# of Autonomous ODE Systems. The R Journal 6(2): 43-51.
# <https://doi.org/10.32614/RJ-2014-023>


# ---------------------------------------------------------------------------
# Shared roxygen parameter block (avoids repetition across 15 functions)
# ---------------------------------------------------------------------------
# @param t Numeric scalar. Time (autonomous system; included for deSolve
#   compatibility).
# @param parameters Not used. Accepted for deSolve compatibility.


# ===========================================================================
# 1D examples (01-05)
# ===========================================================================

#' Example ODE system 1 (1D)
#'
#' The derivative function of example one-dimensional ODE system 1 from
#' Grayling (2014):
#'
#' \deqn{\frac{dy}{dt} = 4 - y^2}
#'
#' Equilibria at \eqn{y^* = \pm 2}: \eqn{y^* = 2} is stable,
#' \eqn{y^* = -2} is unstable.
#'
#' @param t Numeric scalar. Time (autonomous; included for deSolve
#'   compatibility).
#' @param y Numeric vector of length 1. Current state.
#' @param parameters Not used. Accepted for deSolve compatibility.
#'
#' @return A list with one element: a numeric vector of length 1 containing
#'   \eqn{dy/dt}, as required by [deSolve::ode()].
#'
#' @references
#' Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
#' Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
#' \doi{10.32614/RJ-2014-023}
#'
#' @examples
#' ode_example_01(t = 0, y = c(1), parameters = NULL)
#'
#' \dontrun{
#' gg_phase_portrait(ode_example_01, ylim = c(-4, 4))
#' }
#'
#' @export
ode_example_01 <- function(t, y, parameters = NULL) {
  list(c(4 - y[[1L]]^2))
}


#' Example ODE system 2 (1D)
#'
#' The derivative function of example one-dimensional ODE system 2 from
#' Grayling (2014):
#'
#' \deqn{\frac{dy}{dt} = y(1 - y)(2 - y)}
#'
#' Three equilibria: \eqn{y^* = 0} (unstable), \eqn{y^* = 1} (stable),
#' \eqn{y^* = 2} (unstable).
#'
#' @inheritParams ode_example_01
#'
#' @return A list with one element: a numeric vector of length 1 containing
#'   \eqn{dy/dt}, as required by [deSolve::ode()].
#'
#' @references
#' Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
#' Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
#' \doi{10.32614/RJ-2014-023}
#'
#' @examples
#' ode_example_02(t = 0, y = c(0.5), parameters = NULL)
#'
#' \dontrun{
#' gg_phase_portrait(ode_example_02, ylim = c(-0.5, 2.5))
#' }
#'
#' @export
ode_example_02 <- function(t, y, parameters = NULL) {
  list(c(y[[1L]] * (1 - y[[1L]]) * (2 - y[[1L]])))
}


#' Example ODE system 3 (1D)
#'
#' The derivative function of example one-dimensional ODE system 3 from
#' Grayling (2014):
#'
#' \deqn{\frac{dy}{dt} = y^2 - 1}
#'
#' Equilibria at \eqn{y^* = \pm 1}: \eqn{y^* = 1} is unstable,
#' \eqn{y^* = -1} is stable.
#'
#' @inheritParams ode_example_01
#'
#' @return A list with one element: a numeric vector of length 1 containing
#'   \eqn{dy/dt}, as required by [deSolve::ode()].
#'
#' @references
#' Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
#' Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
#' \doi{10.32614/RJ-2014-023}
#'
#' @examples
#' ode_example_03(t = 0, y = c(0), parameters = NULL)
#'
#' \dontrun{
#' gg_phase_portrait(ode_example_03, ylim = c(-3, 3))
#' }
#'
#' @export
ode_example_03 <- function(t, y, parameters = NULL) {
  list(c(y[[1L]]^2 - 1))
}


#' Example ODE system 4 (1D)
#'
#' The derivative function of example one-dimensional ODE system 4 from
#' Grayling (2014):
#'
#' \deqn{\frac{dy}{dt} = y(y - 1)(y + 1)}
#'
#' Three equilibria: \eqn{y^* = -1} (unstable), \eqn{y^* = 0} (stable),
#' \eqn{y^* = 1} (unstable).
#'
#' @inheritParams ode_example_01
#'
#' @return A list with one element: a numeric vector of length 1 containing
#'   \eqn{dy/dt}, as required by [deSolve::ode()].
#'
#' @references
#' Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
#' Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
#' \doi{10.32614/RJ-2014-023}
#'
#' @examples
#' ode_example_04(t = 0, y = c(0.5), parameters = NULL)
#'
#' \dontrun{
#' gg_phase_portrait(ode_example_04, ylim = c(-2, 2))
#' }
#'
#' @export
ode_example_04 <- function(t, y, parameters = NULL) {
  list(c(y[[1L]] * (y[[1L]] - 1) * (y[[1L]] + 1)))
}


#' Example ODE system 5 (1D)
#'
#' The derivative function of example one-dimensional ODE system 5 from
#' Grayling (2014):
#'
#' \deqn{\frac{dy}{dt} = \sin(y)}
#'
#' Infinitely many equilibria at \eqn{y^* = k\pi} for integer \eqn{k}.
#' Even multiples of \eqn{\pi} are unstable; odd multiples are stable.
#' A useful example for illustrating periodic equilibrium structure.
#'
#' @inheritParams ode_example_01
#'
#' @return A list with one element: a numeric vector of length 1 containing
#'   \eqn{dy/dt}, as required by [deSolve::ode()].
#'
#' @references
#' Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
#' Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
#' \doi{10.32614/RJ-2014-023}
#'
#' @examples
#' ode_example_05(t = 0, y = c(pi/2), parameters = NULL)
#'
#' \dontrun{
#' gg_phase_portrait(ode_example_05, ylim = c(-2*pi, 2*pi))
#' }
#'
#' @export
ode_example_05 <- function(t, y, parameters = NULL) {
  list(c(sin(y[[1L]])))
}


# ===========================================================================
# 2D examples (06-15)
# ===========================================================================

#' Example ODE system 6 (2D)
#'
#' The derivative function of example two-dimensional ODE system 6 from
#' Grayling (2014):
#'
#' \deqn{\frac{dx}{dt} = y, \qquad \frac{dy}{dt} = -x}
#'
#' The harmonic oscillator. All trajectories are closed ellipses (circles
#' when plotted on equal-aspect axes) centered on the origin, which is a
#' neutrally stable center.
#'
#' @param t Numeric scalar. Time (autonomous; included for deSolve
#'   compatibility).
#' @param y Numeric vector of length 2. State vector: \eqn{y[1] = x},
#'   \eqn{y[2] = y}.
#' @param parameters Not used. Accepted for deSolve compatibility.
#'
#' @return A list with one element: a numeric vector of length 2 containing
#'   \eqn{(dx/dt, dy/dt)}, as required by [deSolve::ode()].
#'
#' @references
#' Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
#' Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
#' \doi{10.32614/RJ-2014-023}
#'
#' @examples
#' ode_example_06(t = 0, y = c(1, 0), parameters = NULL)
#'
#' \dontrun{
#' gg_flow_field(ode_example_06, xlim = c(-3, 3), ylim = c(-3, 3))
#' }
#'
#' @export
ode_example_06 <- function(t, y, parameters = NULL) {
  list(c(y[[2L]], -y[[1L]]))
}


#' Example ODE system 7 (2D)
#'
#' The derivative function of example two-dimensional ODE system 7 from
#' Grayling (2014):
#'
#' \deqn{\frac{dx}{dt} = -x, \qquad \frac{dy}{dt} = -y}
#'
#' A stable node at the origin. All trajectories decay exponentially to
#' \eqn{(0, 0)} along straight-line paths through the origin.
#'
#' @inheritParams ode_example_06
#'
#' @return A list with one element: a numeric vector of length 2 containing
#'   \eqn{(dx/dt, dy/dt)}, as required by [deSolve::ode()].
#'
#' @references
#' Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
#' Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
#' \doi{10.32614/RJ-2014-023}
#'
#' @examples
#' ode_example_07(t = 0, y = c(2, 1), parameters = NULL)
#'
#' \dontrun{
#' gg_flow_field(ode_example_07, xlim = c(-3, 3), ylim = c(-3, 3))
#' }
#'
#' @export
ode_example_07 <- function(t, y, parameters = NULL) {
  list(c(-y[[1L]], -y[[2L]]))
}


#' Example ODE system 8 (2D)
#'
#' The derivative function of example two-dimensional ODE system 8 from
#' Grayling (2014):
#'
#' \deqn{\frac{dx}{dt} = x, \qquad \frac{dy}{dt} = -y}
#'
#' A saddle point at the origin. The x-axis is the unstable manifold
#' (trajectories move away along it) and the y-axis is the stable manifold
#' (trajectories approach along it).
#'
#' @inheritParams ode_example_06
#'
#' @return A list with one element: a numeric vector of length 2 containing
#'   \eqn{(dx/dt, dy/dt)}, as required by [deSolve::ode()].
#'
#' @references
#' Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
#' Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
#' \doi{10.32614/RJ-2014-023}
#'
#' @examples
#' ode_example_08(t = 0, y = c(1, 1), parameters = NULL)
#'
#' \dontrun{
#' gg_flow_field(ode_example_08, xlim = c(-3, 3), ylim = c(-3, 3))
#' }
#'
#' @export
ode_example_08 <- function(t, y, parameters = NULL) {
  list(c(y[[1L]], -y[[2L]]))
}


#' Example ODE system 9 (2D)
#'
#' The derivative function of example two-dimensional ODE system 9 from
#' Grayling (2014):
#'
#' \deqn{\frac{dx}{dt} = -x + y, \qquad \frac{dy}{dt} = -x - y}
#'
#' A stable spiral at the origin (eigenvalues \eqn{-1 \pm i}).
#' Trajectories spiral inward to \eqn{(0, 0)}.
#'
#' @inheritParams ode_example_06
#'
#' @return A list with one element: a numeric vector of length 2 containing
#'   \eqn{(dx/dt, dy/dt)}, as required by [deSolve::ode()].
#'
#' @references
#' Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
#' Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
#' \doi{10.32614/RJ-2014-023}
#'
#' @examples
#' ode_example_09(t = 0, y = c(2, 0), parameters = NULL)
#'
#' \dontrun{
#' gg_flow_field(ode_example_09, xlim = c(-3, 3), ylim = c(-3, 3))
#' }
#'
#' @export
ode_example_09 <- function(t, y, parameters = NULL) {
  list(c(-y[[1L]] + y[[2L]], -y[[1L]] - y[[2L]]))
}


#' Example ODE system 10 (2D)
#'
#' The derivative function of example two-dimensional ODE system 10 from
#' Grayling (2014):
#'
#' \deqn{\frac{dx}{dt} = x + y, \qquad \frac{dy}{dt} = x - y}
#'
#' A saddle point at the origin (eigenvalues \eqn{\pm\sqrt{2}}).
#'
#' @inheritParams ode_example_06
#'
#' @return A list with one element: a numeric vector of length 2 containing
#'   \eqn{(dx/dt, dy/dt)}, as required by [deSolve::ode()].
#'
#' @references
#' Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
#' Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
#' \doi{10.32614/RJ-2014-023}
#'
#' @examples
#' ode_example_10(t = 0, y = c(1, 0), parameters = NULL)
#'
#' \dontrun{
#' gg_flow_field(ode_example_10, xlim = c(-3, 3), ylim = c(-3, 3))
#' }
#'
#' @export
ode_example_10 <- function(t, y, parameters = NULL) {
  list(c(y[[1L]] + y[[2L]], y[[1L]] - y[[2L]]))
}


#' Example ODE system 11 (2D)
#'
#' The derivative function of example two-dimensional ODE system 11 from
#' Grayling (2014):
#'
#' \deqn{\frac{dx}{dt} = x(3 - x - 2y), \qquad \frac{dy}{dt} = y(2 - x - y)}
#'
#' A nonlinear two-species competition system. Equilibria at \eqn{(0,0)},
#' \eqn{(3,0)}, \eqn{(0,2)}, and the coexistence point \eqn{(1,1)}.
#'
#' @inheritParams ode_example_06
#'
#' @return A list with one element: a numeric vector of length 2 containing
#'   \eqn{(dx/dt, dy/dt)}, as required by [deSolve::ode()].
#'
#' @references
#' Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
#' Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
#' \doi{10.32614/RJ-2014-023}
#'
#' @examples
#' ode_example_11(t = 0, y = c(1, 1), parameters = NULL)
#'
#' \dontrun{
#' gg_flow_field(ode_example_11, xlim = c(0, 4), ylim = c(0, 3))
#' }
#'
#' @export
ode_example_11 <- function(t, y, parameters = NULL) {
  x <- y[[1L]]; v <- y[[2L]]
  list(c(
    x * (3 - x - 2*v),
    v * (2 - x - v)
  ))
}


#' Example ODE system 12 (2D)
#'
#' The derivative function of example two-dimensional ODE system 12 from
#' Grayling (2014):
#'
#' \deqn{\frac{dx}{dt} = y + x(x^2 + y^2 - 1), \qquad
#'       \frac{dy}{dt} = -x + y(x^2 + y^2 - 1)}
#'
#' An unstable spiral at \eqn{(0,0)} and a second equilibrium at
#' \eqn{(1,1)}, a saddle point. Useful for illustrating unstable
#' equilibria and limit cycle-like behavior near the unit circle.
#'
#' @inheritParams ode_example_06
#'
#' @return A list with one element: a numeric vector of length 2 containing
#'   \eqn{(dx/dt, dy/dt)}, as required by [deSolve::ode()].
#'
#' @references
#' Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
#' Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
#' \doi{10.32614/RJ-2014-023}
#'
#' @examples
#' ode_example_12(t = 0, y = c(2, 2), parameters = NULL)
#'
#' \dontrun{
#' gg_flow_field(ode_example_12, xlim = c(-4, 4), ylim = c(-4, 4))
#' }
#'
#' @export
ode_example_12 <- function(t, y, parameters = NULL) {
  x <- y[[1L]]; v <- y[[2L]]
  r2 <- x^2 + v^2
  list(c(
    v + x*(r2 - 1),
   -x + v*(r2 - 1)
  ))
}


#' Example ODE system 13 (2D)
#'
#' The derivative function of example two-dimensional ODE system 13 from
#' Grayling (2014):
#'
#' \deqn{\frac{dx}{dt} = x(1-x) - xy, \qquad
#'       \frac{dy}{dt} = y\!\left(\frac{x}{x+0.5} - 0.5\right)}
#'
#' A predator-prey system with a saturating functional response
#' (Holling type II). Equilibria include the trivial \eqn{(0,0)},
#' the prey-only state \eqn{(1,0)}, and a coexistence equilibrium
#' that may be stable or unstable depending on parameters, giving
#' rise to a limit cycle via Hopf bifurcation.
#'
#' @inheritParams ode_example_06
#'
#' @return A list with one element: a numeric vector of length 2 containing
#'   \eqn{(dx/dt, dy/dt)}, as required by [deSolve::ode()].
#'
#' @references
#' Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
#' Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
#' \doi{10.32614/RJ-2014-023}
#'
#' @examples
#' ode_example_13(t = 0, y = c(0.5, 0.5), parameters = NULL)
#'
#' \dontrun{
#' gg_flow_field(ode_example_13, xlim = c(0, 1.5), ylim = c(0, 1))
#' }
#'
#' @export
ode_example_13 <- function(t, y, parameters = NULL) {
  x <- y[[1L]]; v <- y[[2L]]
  list(c(
    x*(1 - x) - x*v,
    v*(x/(x + 0.5) - 0.5)
  ))
}


#' Example ODE system 14 (2D)
#'
#' The derivative function of example two-dimensional ODE system 14 from
#' Grayling (2014):
#'
#' \deqn{\frac{dx}{dt} = x(1-x^2-y^2), \qquad
#'       \frac{dy}{dt} = y(1-x^2-y^2)}
#'
#' A system with a circle of equilibria at \eqn{x^2 + y^2 = 1} (the unit
#' circle) and an unstable equilibrium at the origin. All trajectories
#' starting outside the origin approach the unit circle, illustrating
#' non-isolated equilibria.
#'
#' @inheritParams ode_example_06
#'
#' @return A list with one element: a numeric vector of length 2 containing
#'   \eqn{(dx/dt, dy/dt)}, as required by [deSolve::ode()].
#'
#' @references
#' Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
#' Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
#' \doi{10.32614/RJ-2014-023}
#'
#' @examples
#' ode_example_14(t = 0, y = c(0.5, 0.5), parameters = NULL)
#'
#' \dontrun{
#' gg_flow_field(ode_example_14, xlim = c(-2, 2), ylim = c(-2, 2))
#' }
#'
#' @export
ode_example_14 <- function(t, y, parameters = NULL) {
  x <- y[[1L]]; v <- y[[2L]]
  r2 <- x^2 + v^2
  list(c(
    x*(1 - r2),
    v*(1 - r2)
  ))
}


#' Example ODE system 15 (2D)
#'
#' The derivative function of example two-dimensional ODE system 15 from
#' Grayling (2014):
#'
#' \deqn{\frac{dx}{dt} = x - y - x(x^2+y^2), \qquad
#'       \frac{dy}{dt} = x + y - y(x^2+y^2)}
#'
#' In polar coordinates this becomes \eqn{dr/dt = r(1-r^2)},
#' \eqn{d\theta/dt = 1}, revealing a stable limit cycle at \eqn{r = 1}
#' (the unit circle) and an unstable spiral at the origin. A classic
#' example for illustrating limit cycles.
#'
#' @inheritParams ode_example_06
#'
#' @return A list with one element: a numeric vector of length 2 containing
#'   \eqn{(dx/dt, dy/dt)}, as required by [deSolve::ode()].
#'
#' @references
#' Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
#' Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
#' \doi{10.32614/RJ-2014-023}
#'
#' @examples
#' ode_example_15(t = 0, y = c(0.5, 0), parameters = NULL)
#'
#' \dontrun{
#' gg_flow_field(ode_example_15, xlim = c(-2, 2), ylim = c(-2, 2))
#' }
#'
#' @export
ode_example_15 <- function(t, y, parameters = NULL) {
  x <- y[[1L]]; v <- y[[2L]]
  r2 <- x^2 + v^2
  list(c(
    x - v - x*r2,
    x + v - v*r2
  ))
}
