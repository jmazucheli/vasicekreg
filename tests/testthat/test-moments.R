# tests/testthat/test-moments.R
#
# Monte Carlo validation of the closed-form marginal mean and variance
# of the Vasicek-type families in vasicekreg.
#
# Rationale: the closed-form marginal moments for ZANVASIM, OANVASIM and
# ZOANVASIM have been verified algebraically. These tests therefore target
# the IMPLEMENTATION of the random generators (r0NVASIM, r1NVASIM,
# r01NVASIM) and of the continuous-component variance formula for NVASIM
# (the bivariate-normal expression Var(X) = Phi2(qmu, qmu; sigma) - mu^2).
# A wrong mixing probability, a swapped nu/tau, or mishandled boundary mass
# would surface here even though the algebra is correct.
#
# Design: instead of an arbitrary fixed tolerance, each check compares the
# standardized discrepancy |estimate - target| / SE against a z-threshold.
# With a fixed seed the tests are deterministic; a genuine bug produces a
# z-score in the tens or hundreds, while legitimate MC noise stays below ~4.

skip_on_cran()               # stochastic + large n: not for CRAN checks
skip_if_not_installed("mvtnorm")

# ---- tuning ---------------------------------------------------------------
# n drives precision (SE ~ 1/sqrt(n)); z_max is the pass/fail threshold on
# the standardized discrepancy. Lower n for a faster local run.
n     <- 1e6
z_max <- 5

# ---- helpers --------------------------------------------------------------

# Sample mean/variance and their Monte Carlo standard errors.
# SE(var_hat) uses the asymptotic form (mu4 - var^2)/n.
mc_moments <- function(y) {
  nn <- length(y)
  m  <- mean(y)
  v  <- stats::var(y)
  m4 <- mean((y - m)^4)
  list(
    mean   = m,
    var    = v,
    se_mean = stats::sd(y) / sqrt(nn),
    se_var  = sqrt(max(m4 - v^2, 0) / nn)
  )
}

# Closed-form variance of the continuous NVASIM component:
# Var(X) = P(Z1 <= qmu, Z2 <= qmu) - mu^2, with corr(Z1, Z2) = sigma.
var_nvasim_cf <- function(mu, sigma) {
  qmu <- stats::qnorm(mu)
  R   <- matrix(c(1, sigma, sigma, 1), 2L, 2L)
  phi2 <- mvtnorm::pmvnorm(upper = c(qmu, qmu), corr = R)
  as.numeric(phi2) - mu^2
}

# Assert that MC estimates of mean and variance are within z_max standard
# errors of the supplied closed-form targets.
expect_moment_match <- function(y, E_target, V_target, label) {
  mc <- mc_moments(y)
  z_mean <- abs(mc$mean - E_target) / mc$se_mean
  z_var  <- abs(mc$var  - V_target) / mc$se_var
  expect_lt(z_mean, z_max,
            label = sprintf("%s: mean z=%.2f (mc=%.6f, cf=%.6f)",
                            label, z_mean, mc$mean, E_target))
  expect_lt(z_var, z_max,
            label = sprintf("%s: var  z=%.2f (mc=%.6f, cf=%.6f)",
                            label, z_var, mc$var, V_target))
}

# ---- continuous component: NVASIM ----------------------------------------

test_that("NVASIM continuous mean and variance match the closed form", {
  grid <- expand.grid(
    mu    = c(0.30, 0.50, 0.70),
    sigma = c(0.30, 0.50),
    KEEP.OUT.ATTRS = FALSE
  )
  for (i in seq_len(nrow(grid))) {
    mu <- grid$mu[i]; sigma <- grid$sigma[i]
    set.seed(20260903 + i)
    y  <- rNVASIM(n, mu = mu, sigma = sigma)
    Vc <- var_nvasim_cf(mu, sigma)
    expect_moment_match(y, E_target = mu, V_target = Vc,
                        label = sprintf("NVASIM(mu=%.2f, sigma=%.2f)", mu, sigma))
  }
})

# ---- zero-augmented: ZANVASIM --------------------------------------------

test_that("ZANVASIM marginal mean and variance match the closed form", {
  grid <- data.frame(
    mu    = c(0.30, 0.50, 0.70, 0.50),
    sigma = c(0.30, 0.50, 0.30, 0.50),
    nu    = c(0.10, 0.25, 0.20, 0.15)
  )
  for (i in seq_len(nrow(grid))) {
    mu <- grid$mu[i]; sigma <- grid$sigma[i]; nu <- grid$nu[i]
    Vc <- var_nvasim_cf(mu, sigma)
    E  <- (1 - nu) * mu
    V  <- (1 - nu) * Vc + nu * (1 - nu) * mu^2
    set.seed(20260903 + 100 + i)
    y <- r0NVASIM(n, mu = mu, sigma = sigma, nu = nu)
    expect_moment_match(y, E, V,
                        label = sprintf("ZANVASIM(mu=%.2f, sigma=%.2f, nu=%.2f)",
                                        mu, sigma, nu))
    # structural check: empirical zero mass ~ nu
    expect_lt(abs(mean(y == 0) - nu) / sqrt(nu * (1 - nu) / n), z_max)
  }
})

# ---- one-augmented: OANVASIM ---------------------------------------------

test_that("OANVASIM marginal mean and variance match the closed form", {
  grid <- data.frame(
    mu    = c(0.30, 0.50, 0.70, 0.60),
    sigma = c(0.30, 0.50, 0.30, 0.40),
    nu    = c(0.10, 0.25, 0.20, 0.15)
  )
  for (i in seq_len(nrow(grid))) {
    mu <- grid$mu[i]; sigma <- grid$sigma[i]; nu <- grid$nu[i]
    Vc <- var_nvasim_cf(mu, sigma)
    E  <- nu + (1 - nu) * mu
    V  <- (1 - nu) * Vc + nu * (1 - nu) * (1 - mu)^2
    set.seed(20260903 + 200 + i)
    y <- r1NVASIM(n, mu = mu, sigma = sigma, nu = nu)
    expect_moment_match(y, E, V,
                        label = sprintf("OANVASIM(mu=%.2f, sigma=%.2f, nu=%.2f)",
                                        mu, sigma, nu))
    expect_lt(abs(mean(y == 1) - nu) / sqrt(nu * (1 - nu) / n), z_max)
  }
})

# ---- zero-and-one-augmented: ZOANVASIM -----------------------------------
# nu = P(Y = 0); tau = P(Y = 1 | Y > 0); mu = E(Y | 0 < Y < 1).

test_that("ZOANVASIM marginal mean and variance match the closed form", {
  grid <- data.frame(
    mu    = c(0.40, 0.50, 0.60, 0.55),
    sigma = c(0.30, 0.50, 0.40, 0.35),
    nu    = c(0.15, 0.20, 0.25, 0.10),
    tau   = c(0.20, 0.30, 0.25, 0.40)
  )
  for (i in seq_len(nrow(grid))) {
    mu <- grid$mu[i]; sigma <- grid$sigma[i]
    nu <- grid$nu[i]; tau <- grid$tau[i]
    Vc <- var_nvasim_cf(mu, sigma)
    E  <- (1 - nu) * (tau + (1 - tau) * mu)
    V  <- (1 - nu) * ((1 - tau) * (Vc + mu^2) + tau) -
          ((1 - nu) * (tau + (1 - tau) * mu))^2
    set.seed(20260903 + 300 + i)
    y <- r01NVASIM(n, mu = mu, sigma = sigma, nu = nu, tau = tau)
    expect_moment_match(y, E, V,
                        label = sprintf(
                          "ZOANVASIM(mu=%.2f, sigma=%.2f, nu=%.2f, tau=%.2f)",
                          mu, sigma, nu, tau))
    # structural checks on the two point masses
    p0 <- nu
    p1 <- (1 - nu) * tau
    expect_lt(abs(mean(y == 0) - p0) / sqrt(p0 * (1 - p0) / n), z_max)
    expect_lt(abs(mean(y == 1) - p1) / sqrt(p1 * (1 - p1) / n), z_max)
  }
})
