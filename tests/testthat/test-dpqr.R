test_that("NVASIM distribution functions are mutually consistent", {
    p <- c(1e-8, 0.10, 0.50, 0.90, 1 - 1e-8)
    mu <- 0.63
    sigma <- 0.27

    expect_equal(
        pNVASIM(qNVASIM(p, mu, sigma), mu, sigma),
        p,
        tolerance = 1e-10
    )
    expect_equal(
        qNVASIM(log(p), mu, sigma, log.p = TRUE),
        qNVASIM(p, mu, sigma),
        tolerance = 1e-12
    )
    expect_equal(
        qNVASIM(log(p), mu, sigma, lower.tail = FALSE, log.p = TRUE),
        qNVASIM(p, mu, sigma, lower.tail = FALSE),
        tolerance = 1e-12
    )
    expect_equal(pNVASIM(c(0, 1), mu, sigma), c(0, 1))
    expect_equal(qNVASIM(c(0, 1), mu, sigma), c(0, 1))
    expect_equal(qNVASIM(c(-Inf, 0), mu, sigma, log.p = TRUE), c(0, 1))
    expect_equal(
        integrate(function(x) dNVASIM(x, mu, sigma), 0, 1)$value,
        1,
        tolerance = 1e-7
    )
})

test_that("NVASIQ distribution functions preserve the selected quantile", {
    p <- c(1e-8, 0.10, 0.50, 0.90, 1 - 1e-8)
    mu <- 0.63
    sigma <- 0.27
    quantile <- 0.20

    expect_equal(
        qNVASIQ(quantile, mu, sigma, quantile = quantile),
        mu,
        tolerance = 1e-12
    )
    expect_equal(
        pNVASIQ(
            qNVASIQ(p, mu, sigma, quantile = quantile),
            mu, sigma, quantile = quantile
        ),
        p,
        tolerance = 1e-10
    )
    expect_equal(
        qNVASIQ(log(p), mu, sigma, quantile = quantile, log.p = TRUE),
        qNVASIQ(p, mu, sigma, quantile = quantile),
        tolerance = 1e-12
    )
    expect_equal(
        pNVASIQ(c(0, 1), mu, sigma, quantile = quantile), c(0, 1)
    )
    expect_equal(
        qNVASIQ(c(0, 1), mu, sigma, quantile = quantile), c(0, 1)
    )
    expect_equal(
        integrate(
            function(x) dNVASIQ(x, mu, sigma, quantile = quantile), 0, 1
        )$value,
        1,
        tolerance = 1e-7
    )
})

test_that("logistic-kernel distribution functions are mutually consistent", {
    p <- c(1e-8, 0.10, 0.50, 0.90, 1 - 1e-8)
    mu <- 0.63
    sigma <- 0.27
    quantile <- 0.20

    expect_equal(
        qLVASIQ(quantile, mu, sigma, quantile = quantile),
        mu,
        tolerance = 1e-12
    )
    expect_equal(
        pLVASIQ(
            qLVASIQ(p, mu, sigma, quantile = quantile),
            mu, sigma, quantile = quantile
        ),
        p,
        tolerance = 1e-10
    )
    expect_equal(
        qLVASIQ(log(p), mu, sigma, quantile = quantile, log.p = TRUE),
        qLVASIQ(p, mu, sigma, quantile = quantile),
        tolerance = 1e-12
    )
    expect_equal(
        pLVASIQ(c(0, 1), mu, sigma, quantile = quantile), c(0, 1)
    )
    expect_equal(
        qLVASIQ(c(0, 1), mu, sigma, quantile = quantile), c(0, 1)
    )
    expect_equal(
        integrate(
            function(x) dLVASIQ(x, mu, sigma, quantile = quantile), 0, 1
        )$value,
        1,
        tolerance = 1e-6
    )
})

test_that("random generators follow the standard n convention", {
    expect_length(rNVASIM(5, 0.5, 0.3), 5)
    expect_length(rNVASIQ(c(4, 8), 0.5, 0.3, 0.5), 2)
    expect_length(rLVASIQ(0, 0.5, 0.3, 0.5), 0)
})

test_that("incompatible parameter lengths are rejected", {
    expect_error(dNVASIM(c(0.2, 0.4, 0.6), c(0.3, 0.7), 0.4))
    expect_error(pNVASIQ(c(0.2, 0.4, 0.6), 0.5, c(0.3, 0.7), 0.5))
    expect_error(qLVASIQ(c(0.2, 0.4, 0.6), 0.5, c(0.3, 0.7), 0.5))
})


test_that("boundary-augmented GAMLSS aliases match numbered dpqr functions", {
    x <- c(0, 0.23, 0.71, 1)
    p <- c(0.05, 0.25, 0.75, 0.95)
    pars <- list(mu = 0.61, sigma = 0.28, nu = 0.17)

    expect_equal(
        dZANVASIM(x, pars$mu, pars$sigma, pars$nu),
        d0NVASIM(x, pars$mu, pars$sigma, pars$nu)
    )
    expect_equal(
        pZANVASIM(x, pars$mu, pars$sigma, pars$nu),
        p0NVASIM(x, pars$mu, pars$sigma, pars$nu)
    )
    expect_equal(
        qZANVASIM(p, pars$mu, pars$sigma, pars$nu),
        q0NVASIM(p, pars$mu, pars$sigma, pars$nu)
    )
    set.seed(101)
    za <- rZANVASIM(20, pars$mu, pars$sigma, pars$nu)
    set.seed(101)
    expect_equal(za, r0NVASIM(20, pars$mu, pars$sigma, pars$nu))

    expect_equal(
        dOANVASIM(x, pars$mu, pars$sigma, pars$nu),
        d1NVASIM(x, pars$mu, pars$sigma, pars$nu)
    )
    expect_equal(
        pOANVASIM(x, pars$mu, pars$sigma, pars$nu),
        p1NVASIM(x, pars$mu, pars$sigma, pars$nu)
    )
    expect_equal(
        qOANVASIM(p, pars$mu, pars$sigma, pars$nu),
        q1NVASIM(p, pars$mu, pars$sigma, pars$nu)
    )
    set.seed(102)
    oa <- rOANVASIM(20, pars$mu, pars$sigma, pars$nu)
    set.seed(102)
    expect_equal(oa, r1NVASIM(20, pars$mu, pars$sigma, pars$nu))

    expect_equal(
        dZOANVASIM(x, pars$mu, pars$sigma, pars$nu, tau = 0.24),
        d01NVASIM(x, pars$mu, pars$sigma, pars$nu, tau = 0.24)
    )
    expect_equal(
        pZOANVASIM(x, pars$mu, pars$sigma, pars$nu, tau = 0.24),
        p01NVASIM(x, pars$mu, pars$sigma, pars$nu, tau = 0.24)
    )
    expect_equal(
        qZOANVASIM(p, pars$mu, pars$sigma, pars$nu, tau = 0.24),
        q01NVASIM(p, pars$mu, pars$sigma, pars$nu, tau = 0.24)
    )
    set.seed(103)
    zoa <- rZOANVASIM(
        20, pars$mu, pars$sigma, pars$nu, tau = 0.24
    )
    set.seed(103)
    expect_equal(
        zoa,
        r01NVASIM(20, pars$mu, pars$sigma, pars$nu, tau = 0.24)
    )
})
