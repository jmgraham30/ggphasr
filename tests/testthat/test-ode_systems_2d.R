# tests/testthat/test-ode_systems_2d.R
#
# Unit tests for the eight built-in 2D ODE systems.
#
# Each function is tested for:
#   (1) Correct return structure (list containing a length-2 numeric vector)
#   (2) Correct derivative values at analytically known points
#   (3) Correct behavior at equilibrium points (both derivatives = 0)
#   (4) Named and positional parameter passing both work
#   (5) Default parameters produce a valid result
#   (6) Key qualitative properties (e.g., sign of derivatives in known regions)


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

expect_ode2_structure <- function(result) {
  expect_type(result, "list")
  expect_length(result, 1L)
  expect_type(result[[1L]], "double")
  expect_length(result[[1L]], 2L)
}


# ===========================================================================
# ode_lotka_volterra()
# ===========================================================================

test_that("ode_lotka_volterra() returns correct structure", {
  result <- ode_lotka_volterra(t = 0, y = c(1, 1),
                                parameters = c(alpha=1, beta=0.5,
                                               delta=0.5, gamma=1))
  expect_ode2_structure(result)
})

test_that("ode_lotka_volterra() interior equilibrium has zero derivatives", {
  # Equilibrium at (gamma/delta, alpha/beta) = (2, 2)
  result <- ode_lotka_volterra(t = 0, y = c(2, 2),
                                parameters = c(alpha=1, beta=0.5,
                                               delta=0.5, gamma=1))
  expect_equal(result[[1L]], c(0, 0), tolerance = 1e-10)
})

test_that("ode_lotka_volterra() origin is an equilibrium", {
  result <- ode_lotka_volterra(t = 0, y = c(0, 0),
                                parameters = c(alpha=1, beta=0.5,
                                               delta=0.5, gamma=1))
  expect_equal(result[[1L]], c(0, 0))
})

test_that("ode_lotka_volterra() computes derivatives correctly at (1,1)", {
  # alpha=1, beta=0.5, delta=0.5, gamma=1, x=1, y=1
  # dx/dt = 1*1 - 0.5*1*1 = 0.5
  # dy/dt = 0.5*1*1 - 1*1 = -0.5
  result <- ode_lotka_volterra(t = 0, y = c(1, 1),
                                parameters = c(alpha=1, beta=0.5,
                                               delta=0.5, gamma=1))
  expect_equal(result[[1L]], c(0.5, -0.5))
})

test_that("ode_lotka_volterra() works with positional parameters", {
  r_named <- ode_lotka_volterra(t=0, y=c(1,1),
                                 parameters=c(alpha=1, beta=0.5,
                                              delta=0.5, gamma=1))
  r_pos   <- ode_lotka_volterra(t=0, y=c(1,1),
                                 parameters=c(1, 0.5, 0.5, 1))
  expect_equal(r_named[[1L]], r_pos[[1L]])
})

test_that("ode_lotka_volterra() default parameters produce valid output", {
  result <- ode_lotka_volterra(t = 0, y = c(1, 1))
  expect_ode2_structure(result)
})


# ===========================================================================
# ode_sir()
# ===========================================================================

test_that("ode_sir() returns correct structure (proportions)", {
  result <- ode_sir(t = 0, y = c(0.99, 0.01),
                   parameters = c(beta = 0.5, gamma = 0.1))
  expect_ode2_structure(result)
})

test_that("ode_sir() proportions: dS/dt + dI/dt = -(gamma * I)", {
  # In the SI reduction, dS/dt + dI/dt = -gamma*I (R is absorbing)
  beta <- 0.5; gamma <- 0.1; S <- 0.8; I <- 0.1
  result <- ode_sir(t = 0, y = c(S, I),
                   parameters = c(beta = beta, gamma = gamma))
  dS <- result[[1L]][1L]; dI <- result[[1L]][2L]
  expect_equal(dS + dI, -gamma * I, tolerance = 1e-12)
})

test_that("ode_sir() proportions: disease-free equilibrium (I=0) has dS=dI=0", {
  result <- ode_sir(t = 0, y = c(1, 0),
                   parameters = c(beta = 0.5, gamma = 0.1))
  expect_equal(result[[1L]], c(0, 0))
})

test_that("ode_sir() proportions: epidemic grows when R0 > 1 and S near 1", {
  # R0 = beta/gamma = 5 > 1, so I should be increasing near S=1
  result <- ode_sir(t = 0, y = c(0.99, 0.01),
                   parameters = c(beta = 0.5, gamma = 0.1))
  expect_gt(result[[1L]][2L], 0)   # dI/dt > 0
})

test_that("ode_sir() counts formulation produces same proportional behavior when N=1", {
  params_prop   <- c(beta = 0.5, gamma = 0.1)
  params_counts <- list(beta = 0.5, gamma = 0.1, N = 1, scale = "counts")
  r_prop   <- ode_sir(t=0, y=c(0.8, 0.1), parameters = params_prop)
  r_counts <- ode_sir(t=0, y=c(0.8, 0.1), parameters = params_counts)
  expect_equal(r_prop[[1L]], r_counts[[1L]], tolerance = 1e-12)
})

test_that("ode_sir() counts: derivatives scale correctly with N", {
  # With N=1000, S=800, I=100: dS/dt = -beta*S*I/N = -0.5*800*100/1000 = -40
  result <- ode_sir(t = 0, y = c(800, 100),
                   parameters = list(beta = 0.5, gamma = 0.1,
                                     N = 1000, scale = "counts"))
  expect_equal(result[[1L]][1L], -40, tolerance = 1e-10)
})

test_that("ode_sir() default parameters produce valid output", {
  result <- ode_sir(t = 0, y = c(0.99, 0.01))
  expect_ode2_structure(result)
})


# ===========================================================================
# ode_van_der_pol()
# ===========================================================================

test_that("ode_van_der_pol() returns correct structure", {
  result <- ode_van_der_pol(t = 0, y = c(1, 0), parameters = c(mu = 1))
  expect_ode2_structure(result)
})

test_that("ode_van_der_pol() when mu=0 reduces to harmonic oscillator", {
  # mu=0: dx/dt = y, dy/dt = 0*(...)*y - x = -x
  result <- ode_van_der_pol(t = 0, y = c(2, 3), parameters = c(mu = 0))
  expect_equal(result[[1L]], c(3, -2))
})

test_that("ode_van_der_pol() computes derivatives correctly at (2, 0)", {
  # x=2, y=0, mu=1: dx/dt=0, dy/dt = 1*(1-4)*0 - 2 = -2
  result <- ode_van_der_pol(t = 0, y = c(2, 0), parameters = c(mu = 1))
  expect_equal(result[[1L]], c(0, -2))
})

test_that("ode_van_der_pol() origin is an equilibrium for any mu", {
  result <- ode_van_der_pol(t = 0, y = c(0, 0), parameters = c(mu = 2))
  expect_equal(result[[1L]], c(0, 0))
})

test_that("ode_van_der_pol() works with positional parameters", {
  r_named <- ode_van_der_pol(t=0, y=c(1,1), parameters=c(mu=2))
  r_pos   <- ode_van_der_pol(t=0, y=c(1,1), parameters=c(2))
  expect_equal(r_named[[1L]], r_pos[[1L]])
})

test_that("ode_van_der_pol() default parameters produce valid output", {
  result <- ode_van_der_pol(t = 0, y = c(1, 0))
  expect_ode2_structure(result)
})


# ===========================================================================
# ode_simple_pendulum()
# ===========================================================================

test_that("ode_simple_pendulum() returns correct structure", {
  result <- ode_simple_pendulum(t=0, y=c(pi/4, 0),
                                parameters=c(gL=1, b=0))
  expect_ode2_structure(result)
})

test_that("ode_simple_pendulum() equilibrium at theta=0, omega=0", {
  result <- ode_simple_pendulum(t=0, y=c(0, 0),
                                parameters=c(gL=1, b=0))
  expect_equal(result[[1L]], c(0, 0), tolerance = 1e-15)
})

test_that("ode_simple_pendulum() unstable equilibrium at theta=pi, omega=0", {
  result <- ode_simple_pendulum(t=0, y=c(pi, 0),
                                parameters=c(gL=1, b=0))
  # sin(pi) is numerically near 0
  expect_equal(result[[1L]][1L], 0)
  expect_equal(result[[1L]][2L], 0, tolerance = 1e-14)
})

test_that("ode_simple_pendulum() undamped: dtheta/dt = omega", {
  omega <- 1.5
  result <- ode_simple_pendulum(t=0, y=c(pi/6, omega),
                                parameters=c(gL=1, b=0))
  expect_equal(result[[1L]][1L], omega)
})

test_that("ode_simple_pendulum() damping reduces angular velocity", {
  # With damping, domega/dt should be less (more negative) than without
  undamped <- ode_simple_pendulum(t=0, y=c(pi/4, 1),
                                  parameters=c(gL=1, b=0))
  damped   <- ode_simple_pendulum(t=0, y=c(pi/4, 1),
                                  parameters=c(gL=1, b=0.5))
  expect_lt(damped[[1L]][2L], undamped[[1L]][2L])
})

test_that("ode_simple_pendulum() default parameters produce valid output", {
  result <- ode_simple_pendulum(t=0, y=c(pi/4, 0))
  expect_ode2_structure(result)
})


# ===========================================================================
# ode_competition()
# ===========================================================================

test_that("ode_competition() returns correct structure", {
  result <- ode_competition(t=0, y=c(5, 5),
                             parameters=c(r1=1, r2=1, K1=10, K2=10,
                                          a12=0.5, a21=0.5))
  expect_ode2_structure(result)
})

test_that("ode_competition() trivial equilibrium at (0,0)", {
  result <- ode_competition(t=0, y=c(0, 0),
                             parameters=c(r1=1, r2=1, K1=10, K2=10,
                                          a12=0.5, a21=0.5))
  expect_equal(result[[1L]], c(0, 0))
})

test_that("ode_competition() single-species equilibria at (K1,0) and (0,K2)", {
  r1 <- ode_competition(t=0, y=c(10, 0),
                         parameters=c(r1=1, r2=1, K1=10, K2=10,
                                      a12=0.5, a21=0.5))
  r2 <- ode_competition(t=0, y=c(0, 10),
                         parameters=c(r1=1, r2=1, K1=10, K2=10,
                                      a12=0.5, a21=0.5))
  expect_equal(r1[[1L]][1L], 0, tolerance = 1e-10)
  expect_equal(r2[[1L]][2L], 0, tolerance = 1e-10)
})

test_that("ode_competition() coexistence equilibrium has zero derivatives", {
  # With symmetric params (r1=r2, K1=K2=K, a12=a21=a),
  # coexistence equilibrium is at N1* = N2* = K/(1+a)
  K <- 10; a <- 0.5; N_star <- K / (1 + a)  # = 6.667
  result <- ode_competition(t=0, y=c(N_star, N_star),
                             parameters=c(r1=1, r2=1, K1=K, K2=K,
                                          a12=a, a21=a))
  expect_equal(result[[1L]], c(0, 0), tolerance = 1e-10)
})

test_that("ode_competition() default parameters produce valid output", {
  result <- ode_competition(t=0, y=c(5, 5))
  expect_ode2_structure(result)
})


# ===========================================================================
# ode_toggle()
# ===========================================================================

test_that("ode_toggle() returns correct structure", {
  result <- ode_toggle(t=0, y=c(2, 0.5),
                       parameters=c(alpha1=3, alpha2=3, beta=2, gamma=2))
  expect_ode2_structure(result)
})

test_that("ode_toggle() computes du/dt correctly", {
  # alpha1=3, alpha2=3, beta=2, gamma=2, u=1, v=1
  # du/dt = 3/(1+1^2) - 1 = 1.5 - 1 = 0.5
  # dv/dt = 3/(1+1^2) - 1 = 0.5
  result <- ode_toggle(t=0, y=c(1, 1),
                       parameters=c(alpha1=3, alpha2=3, beta=2, gamma=2))
  expect_equal(result[[1L]], c(0.5, 0.5))
})

test_that("ode_toggle() is symmetric under u<->v swap when alpha1=alpha2, beta=gamma", {
  r1 <- ode_toggle(t=0, y=c(2, 0.5),
                   parameters=c(alpha1=3, alpha2=3, beta=2, gamma=2))
  r2 <- ode_toggle(t=0, y=c(0.5, 2),
                   parameters=c(alpha1=3, alpha2=3, beta=2, gamma=2))
  # Swapping u and v should swap du/dt and dv/dt
  expect_equal(r1[[1L]][1L], r2[[1L]][2L])
  expect_equal(r1[[1L]][2L], r2[[1L]][1L])
})

test_that("ode_toggle() default parameters produce valid output", {
  result <- ode_toggle(t=0, y=c(1, 1))
  expect_ode2_structure(result)
})


# ===========================================================================
# ode_morris_lecar()
# ===========================================================================

# Standard parameter set from Morris & Lecar (1981)
ml_params <- c(I=0, C=20, gCa=4.4, gK=8, gL=2,
                VCa=120, VK=-84, VL=-60,
                V1=-1.2, V2=18, V3=2, V4=30, phi=0.04)

test_that("ode_morris_lecar() returns correct structure", {
  result <- ode_morris_lecar(t=0, y=c(-60, 0), parameters=ml_params)
  expect_ode2_structure(result)
})

test_that("ode_morris_lecar() dN/dt = 0 when N = N_inf(V)", {
  # When N equals its steady-state value N_inf(V), dN/dt should be 0
  V    <- -20
  V3   <- 2; V4 <- 30
  N_inf <- 0.5 * (1 + tanh((V - V3) / V4))
  result <- ode_morris_lecar(t=0, y=c(V, N_inf), parameters=ml_params)
  expect_equal(result[[1L]][2L], 0, tolerance = 1e-12)
})

test_that("ode_morris_lecar() N increases when below N_inf(V)", {
  V    <- -20
  V3   <- 2; V4 <- 30
  N_inf <- 0.5 * (1 + tanh((V - V3) / V4))
  # Set N below N_inf -> dN/dt should be positive
  result <- ode_morris_lecar(t=0, y=c(V, N_inf * 0.5), parameters=ml_params)
  expect_gt(result[[1L]][2L], 0)
})

test_that("ode_morris_lecar() N decreases when above N_inf(V)", {
  V     <- -20
  V3    <- 2; V4 <- 30
  N_inf <- 0.5 * (1 + tanh((V - V3) / V4))
  result <- ode_morris_lecar(t=0, y=c(V, min(N_inf * 1.5, 0.99)),
                              parameters=ml_params)
  expect_lt(result[[1L]][2L], 0)
})

test_that("ode_morris_lecar() applied current shifts dV/dt", {
  # Increasing I should increase dV/dt at the same state
  r_low  <- ode_morris_lecar(t=0, y=c(-60, 0.1),
                              parameters=c(ml_params[names(ml_params) != "I"],
                                           I = 0))
  r_high <- ode_morris_lecar(t=0, y=c(-60, 0.1),
                              parameters=c(ml_params[names(ml_params) != "I"],
                                           I = 100))
  expect_gt(r_high[[1L]][1L], r_low[[1L]][1L])
})

test_that("ode_morris_lecar() default parameters produce valid output", {
  result <- ode_morris_lecar(t=0, y=c(-60, 0))
  expect_ode2_structure(result)
})


# ===========================================================================
# ode_lindemann()
# ===========================================================================

test_that("ode_lindemann() returns correct structure", {
  result <- ode_lindemann(t=0, y=c(2, 0.1),
                          parameters=c(k1=1, k_1=1, k2=0.5))
  expect_ode2_structure(result)
})

test_that("ode_lindemann() mass is not conserved (product sinks A)", {
  # d[A]/dt + d[A*]/dt = -k2*[A*] < 0 (total A + A* decreases as P forms)
  k2 <- 0.5; As <- 0.1
  result <- ode_lindemann(t=0, y=c(2, As),
                          parameters=c(k1=1, k_1=1, k2=k2))
  total_rate <- sum(result[[1L]])
  expect_equal(total_rate, -k2 * As, tolerance = 1e-12)
})

test_that("ode_lindemann() computes derivatives correctly at known point", {
  # k1=1, k_1=1, k2=0.5, A=2, As=0
  # d[A]/dt  = -1*4 + 1*0*2 = -4
  # d[A*]/dt =  1*4 - 1*0*2 - 0.5*0 = 4
  result <- ode_lindemann(t=0, y=c(2, 0),
                          parameters=c(k1=1, k_1=1, k2=0.5))
  expect_equal(result[[1L]], c(-4, 4))
})

test_that("ode_lindemann() trivial equilibrium at (0, 0)", {
  result <- ode_lindemann(t=0, y=c(0, 0),
                          parameters=c(k1=1, k_1=1, k2=0.5))
  expect_equal(result[[1L]], c(0, 0))
})

test_that("ode_lindemann() works with positional parameters", {
  r_named <- ode_lindemann(t=0, y=c(2, 0.1),
                            parameters=c(k1=1, k_1=1, k2=0.5))
  r_pos   <- ode_lindemann(t=0, y=c(2, 0.1),
                            parameters=c(1, 1, 0.5))
  expect_equal(r_named[[1L]], r_pos[[1L]])
})

test_that("ode_lindemann() default parameters produce valid output", {
  result <- ode_lindemann(t=0, y=c(2, 0.1))
  expect_ode2_structure(result)
})


# ===========================================================================
# Cross-function tests
# ===========================================================================

test_that("all 2D ODE functions pass .validate_ode()", {
  cases <- list(
    list(fn=ode_lotka_volterra,  params=c(alpha=1, beta=0.5, delta=0.5, gamma=1)),
    list(fn=ode_sir,             params=c(beta=0.5, gamma=0.1)),
    list(fn=ode_van_der_pol,     params=c(mu=1)),
    list(fn=ode_simple_pendulum, params=c(gL=1, b=0)),
    list(fn=ode_competition,     params=c(r1=1, r2=1, K1=10, K2=10,
                                          a12=0.5, a21=0.5)),
    list(fn=ode_toggle,          params=c(alpha1=3, alpha2=3, beta=2, gamma=2)),
    list(fn=ode_morris_lecar,    params=ml_params),
    list(fn=ode_lindemann,       params=c(k1=1, k_1=1, k2=0.5))
  )

  for (item in cases) {
    normalized <- ggphasr:::.normalize_ode(item$fn, system = "two.dim")
    expect_true(
      ggphasr:::.validate_ode(normalized,
                               system     = "two.dim",
                               parameters = item$params)
    )
  }
})

test_that("all 2D ODE functions are time-invariant (autonomous)", {
  cases <- list(
    list(fn=ode_lotka_volterra,  y=c(1,1), params=c(alpha=1,beta=0.5,delta=0.5,gamma=1)),
    list(fn=ode_van_der_pol,     y=c(1,1), params=c(mu=1)),
    list(fn=ode_simple_pendulum, y=c(1,1), params=c(gL=1, b=0))
  )

  for (item in cases) {
    r0   <- item$fn(t=0,   y=item$y, parameters=item$params)
    r100 <- item$fn(t=100, y=item$y, parameters=item$params)
    expect_equal(r0[[1L]], r100[[1L]])
  }
})
