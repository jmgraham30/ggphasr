# ggphasr-package.R
#
# Package-level documentation and global variable declarations.

#' ggphasr: Phase Plane Analysis of ODE Systems Using ggplot2
#'
#' @description
#' `ggphasr` provides tools for qualitative analysis of one- and
#' two-dimensional autonomous ordinary differential equation (ODE) systems
#' using phase plane methods. All visualizations are produced with
#' \pkg{ggplot2} and its extensions, making plots fully customizable with
#' standard \pkg{ggplot2} syntax.
#'
#' @section ODE conventions:
#' All functions accept ODE systems in one of two calling conventions:
#'
#' **Convention A** (deSolve-compatible, primary):
#' ```r
#' f <- function(t, y, parameters) {
#'   list(c(dy1, dy2))
#' }
#' ```
#'
#' **Convention B** (simplified):
#' ```r
#' # 2D:
#' f <- function(x, y, parameters = NULL) c(dx, dy)
#' # 1D:
#' f <- function(y, parameters = NULL) dy
#' ```
#'
#' Convention B functions are automatically detected and wrapped into
#' Convention A internally. All built-in ODE systems use Convention A.
#'
#' @section Main functions:
#'
#' **Plotting (return ggplot objects or layer lists):**
#' \describe{
#'   \item{[gg_flow_field()]}{Direction/velocity field arrows on a grid}
#'   \item{[gg_nullclines()]}{Zero-isocline curves (x- and y-nullclines)}
#'   \item{[gg_trajectory()]}{Numerically integrated solution paths}
#'   \item{[gg_phase_portrait()]}{1D phase line with directional arrows}
#'   \item{[gg_time_series()]}{Time-series plots of state variables}
#'   \item{[gg_manifolds()]}{Stable/unstable manifolds of saddle points}
#'   \item{[gg_phase_plane()]}{All-in-one phase plane analysis wrapper}
#' }
#'
#' **Analysis (return data structures, no plots):**
#' \describe{
#'   \item{[find_equilibrium()]}{Numerical root-finding for equilibria}
#'   \item{[classify_equilibrium()]}{Trace-determinant stability classification}
#' }
#'
#' **Built-in ODE systems:**
#' \describe{
#'   \item{1D models}{[ode_exponential()], [ode_logistic()],
#'     [ode_monomolecular()], [ode_von_bertalanffy()]}
#'   \item{2D models}{[ode_lotka_volterra()], [ode_sir()],
#'     [ode_van_der_pol()], [ode_simple_pendulum()], [ode_competition()],
#'     [ode_toggle()], [ode_morris_lecar()], [ode_lindemann()]}
#'   \item{Textbook examples}{[ode_example_01()] through [ode_example_15()]}
#' }
#'
#' @section Typical workflow:
#' ```r
#' library(ggphasr)
#' library(ggplot2)
#'
#' lv_params <- c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1)
#'
#' # Quick all-in-one analysis
#' result <- gg_phase_plane(ode_lotka_volterra,
#'                           xlim = c(0, 5), ylim = c(0, 5),
#'                           parameters = lv_params)
#' result$plot
#' result$equilibria
#'
#' # Composable layer-by-layer approach
#' gg_flow_field(ode_lotka_volterra,
#'               xlim = c(0, 5), ylim = c(0, 5),
#'               parameters = lv_params) +
#'   gg_nullclines(ode_lotka_volterra,
#'                 xlim = c(0, 5), ylim = c(0, 5),
#'                 parameters = lv_params) +
#'   gg_trajectory(ode_lotka_volterra,
#'                 y0 = c(1, 1), xlim = c(0, 5), ylim = c(0, 5),
#'                 parameters = lv_params)
#' ```
#'
#' @references
#' Grayling MJ (2014). phaseR: An R Package for Phase Plane Analysis of
#' Autonomous ODE Systems. *The R Journal* 6(2): 43-51.
#' \doi{10.32614/RJ-2014-023}
#'
#' @keywords internal
"_PACKAGE"

# ---------------------------------------------------------------------------
# Global variable declarations
# ---------------------------------------------------------------------------
# Suppress R CMD CHECK notes about .data pronoun from rlang/ggplot2.
# All aes() mappings in ggphasr use .data$column to avoid ambiguity,
# but R CMD CHECK does not recognize .data as a known binding.

utils::globalVariables(".data")
