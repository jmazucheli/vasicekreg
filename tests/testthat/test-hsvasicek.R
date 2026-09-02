test_that("HSVASIQ distribution functions are mutually consistent", {
    p <- c(1e-8, 0.10, 0.50, 0.90, 1 - 1e-8)
    mu <- 0.63
    sigma <- 0.27
    tau <- 0.20

    expect_equal(qHSVASIQ(tau, mu, sigma, tau), mu, tolerance = 1e-12)
    expect_equal(
        pHSVASIQ(qHSVASIQ(p, mu, sigma, tau), mu, sigma, tau),
        p,
        tolerance = 1e-9
    )
    expect_equal(
        qHSVASIQ(log(p), mu, sigma, tau, log.p = TRUE),
        qHSVASIQ(p, mu, sigma, tau),
        tolerance = 1e-12
    )
    expect_equal(
        qHSVASIQ(log(p), mu, sigma, tau,
                 lower.tail = FALSE, log.p = TRUE),
        qHSVASIQ(p, mu, sigma, tau, lower.tail = FALSE),
        tolerance = 1e-12
    )
    expect_equal(pHSVASIQ(c(0, 1), mu, sigma, tau), c(0, 1))
    expect_equal(qHSVASIQ(c(0, 1), mu, sigma, tau), c(0, 1))
    expect_equal(
        integrate(function(x) dHSVASIQ(x, mu, sigma, tau), 0, 1)$value,
        1,
        tolerance = 1e-6
    )
})

test_that("HSVASIQ analytical derivatives agree with numerical derivatives", {
    y <- 0.37
    mu <- 0.69
    sigma <- 0.37
    tau <- 0.81
    loglik <- function(m, s) dHSVASIQ(y, m, s, tau, log = TRUE)

    with_global_tau(tau, {
        family <- HSVASIQ()
        expect_equal(
            family$dldm(y, mu, sigma),
            hs_first_derivative(function(m) loglik(m, sigma), mu),
            tolerance = 1e-5
        )
        expect_equal(
            family$d2ldm2(y, mu, sigma),
            hs_second_derivative(function(m) loglik(m, sigma), mu),
            tolerance = 1e-4
        )
        expect_equal(
            family$dldd(y, mu, sigma),
            hs_first_derivative(function(s) loglik(mu, s), sigma),
            tolerance = 1e-5
        )
        expect_equal(
            family$d2ldd2(y, mu, sigma),
            hs_second_derivative(function(s) loglik(mu, s), sigma),
            tolerance = 1e-4
        )
        expect_equal(
            family$d2ldmdd(y, mu, sigma),
            hs_cross_derivative(loglik, mu, sigma),
            tolerance = 1e-4
        )
    })
})

test_that("HSVASIQ family moments agree with quantile integration", {
    mu <- 0.63
    sigma <- 0.27
    tau <- 0.20

    with_global_tau(tau, {
        family <- HSVASIQ()
        numerical_mean <- integrate(
            function(p) qHSVASIQ(p, mu, sigma, tau),
            0, 1, subdivisions = 200L, rel.tol = 1e-8
        )$value
        numerical_variance <- integrate(
            function(p) {
                (qHSVASIQ(p, mu, sigma, tau) - numerical_mean)^2
            },
            0, 1, subdivisions = 200L, rel.tol = 1e-8
        )$value

        expect_equal(family$mean(mu, sigma), numerical_mean, tolerance = 1e-7)
        expect_equal(
            family$variance(mu, sigma),
            numerical_variance,
            tolerance = 1e-6
        )
        expect_identical(family$rqres[[1]][["tau"]], quote(tau))
    })
})

test_that("HSVASIQ follows vectorization and random-generation conventions", {
    expect_length(rHSVASIQ(5, 0.5, 0.3, 0.5), 5)
    expect_length(rHSVASIQ(c(4, 8), 0.5, 0.3, 0.5), 2)
    expect_error(
        dHSVASIQ(c(0.2, 0.4, 0.6), 0.5, c(0.3, 0.7), 0.5)
    )
})

test_that("HSVASIQ requires a valid global quantile level", {
    expect_error(HSVASIQ(tau = 0.25), "unused argument")
    with_global_tau(0, expect_error(HSVASIQ(), "strictly between"))
    with_global_tau(1, expect_error(HSVASIQ(), "strictly between"))
    with_global_tau(
        c(0.25, 0.75),
        expect_error(HSVASIQ(), "single number")
    )
})

test_that("HSVASIQ fits an intercept-only GAMLSS model", {
    set.seed(20260902)
    tau <- 0.25
    y <- rHSVASIQ(120, mu = 0.55, sigma = 0.25, tau = tau)
    control <- gamlss::gamlss.control(n.cyc = 75, trace = FALSE)

    fit <- with_global_tau(
        tau,
        gamlss::gamlss(
            y ~ 1,
            sigma.formula = ~ 1,
            family = HSVASIQ(),
            control = control
        )
    )

    expect_s3_class(fit, "gamlss")
    expect_equal(
        as.numeric(fitted(fit, what = "mu")[1]),
        0.55,
        tolerance = 0.15
    )
})

test_that("envelope simulation retains the fixed HSVASIQ quantile level", {
    set.seed(20260902)
    tau <- 0.25
    y <- rHSVASIQ(60, mu = 0.55, sigma = 0.25, tau = tau)

    with_global_tau(tau, {
        fit <- gamlss::gamlss(
            y ~ 1,
            sigma.formula = ~ 1,
            family = HSVASIQ(),
            control = gamlss::gamlss.control(n.cyc = 75, trace = FALSE)
        )
        fitted_mu <- as.numeric(fitted(fit, what = "mu"))
        fitted_sigma <- as.numeric(fitted(fit, what = "sigma"))

        set.seed(91)
        simulated <- .envelope_simulate_response(fit)
        set.seed(91)
        expected <- rHSVASIQ(
            length(y), fitted_mu, fitted_sigma, tau = tau
        )

        expect_identical(simulated, expected)
    })
})
hs_first_derivative <- function(f, x, h = 1e-6) {
    (f(x + h) - f(x - h)) / (2 * h)
}

hs_second_derivative <- function(f, x, h = 1e-4) {
    (f(x + h) - 2 * f(x) + f(x - h)) / h^2
}

hs_cross_derivative <- function(f, x, y, h = 1e-4) {
    (f(x + h, y + h) - f(x + h, y - h) -
         f(x - h, y + h) + f(x - h, y - h)) / (4 * h^2)
}
