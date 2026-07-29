test_that("NVASIM family moments agree with numerical integration", {
    mu <- 0.63
    sigma <- 0.27
    family <- NVASIM()

    numerical_mean <- integrate(
        function(x) x * dNVASIM(x, mu, sigma), 0, 1
    )$value
    numerical_variance <- integrate(
        function(x) (x - numerical_mean)^2 * dNVASIM(x, mu, sigma),
        0, 1
    )$value

    expect_equal(family$mean(mu), numerical_mean, tolerance = 1e-7)
    expect_equal(
        family$variance(mu, sigma),
        numerical_variance,
        tolerance = 1e-6
    )
})

test_that("NVASIQ uses global tau and returns correct moments", {
    mu <- 0.63
    sigma <- 0.27
    tau <- 0.20

    with_global_tau(tau, {
        family <- NVASIQ()
        numerical_mean <- integrate(
            function(x) x * dNVASIQ(x, mu, sigma, tau), 0, 1
        )$value
        numerical_variance <- integrate(
            function(x) (x - numerical_mean)^2 *
                dNVASIQ(x, mu, sigma, tau),
            0, 1
        )$value

        expect_equal(
            family$mean(mu, sigma),
            numerical_mean,
            tolerance = 1e-7
        )
        expect_equal(
            family$variance(mu, sigma),
            numerical_variance,
            tolerance = 1e-6
        )
        expect_identical(
            family$rqres[[1]][["tau"]],
            quote(tau)
        )
    })
})

test_that("LVASIQ uses global tau and returns correct moments", {
    mu <- 0.63
    sigma <- 0.27
    tau <- 0.20

    with_global_tau(tau, {
        family <- LVASIQ()
        numerical_mean <- integrate(
            function(p) qLVASIQ(p, mu, sigma, tau),
            lower = 0,
            upper = 1,
            subdivisions = 200L,
            rel.tol = 1e-8
        )$value
        numerical_variance <- integrate(
            function(p) {
                (qLVASIQ(p, mu, sigma, tau) - numerical_mean)^2
            },
            lower = 0,
            upper = 1,
            subdivisions = 200L,
            rel.tol = 1e-8
        )$value

        expect_equal(
            family$mean(mu, sigma),
            numerical_mean,
            tolerance = 1e-7
        )
        expect_equal(
            family$variance(mu, sigma),
            numerical_variance,
            tolerance = 1e-6
        )
        expect_identical(
            family$rqres[[1]][["tau"]],
            quote(tau)
        )
    })
})

test_that("family moment functions are vectorized", {
    mu <- c(0.30, 0.70)
    sigma <- c(0.20, 0.80)

    expect_length(NVASIM()$variance(mu, sigma), 2)
    expect_true(all(is.finite(NVASIM()$variance(mu, sigma))))
    with_global_tau(0.25, {
        expect_length(NVASIQ()$mean(mu, sigma), 2)
        expect_length(NVASIQ()$variance(mu, sigma), 2)
        expect_true(all(NVASIQ()$variance(mu, sigma) > 0))
        expect_length(LVASIQ()$mean(mu, sigma), 2)
        expect_length(LVASIQ()$variance(mu, sigma), 2)
        expect_true(all(LVASIQ()$variance(mu, sigma) > 0))
    })
})

test_that("NVASIQ requires a valid global quantile level", {
    expect_error(NVASIQ(tau = 0.25), "unused argument")
    with_global_tau(0, expect_error(NVASIQ(), "strictly between"))
    with_global_tau(1, expect_error(NVASIQ(), "strictly between"))
    with_global_tau(
        c(0.25, 0.75),
        expect_error(NVASIQ(), "single number")
    )
})

test_that("LVASIQ requires a valid global quantile level", {
    expect_error(LVASIQ(tau = 0.25), "unused argument")
    with_global_tau(0, expect_error(LVASIQ(), "strictly between"))
    with_global_tau(1, expect_error(LVASIQ(), "strictly between"))
    with_global_tau(
        c(0.25, 0.75),
        expect_error(LVASIQ(), "single number")
    )
})

test_that("normal-kernel families fit intercept-only GAMLSS models", {
    set.seed(123)
    y_mean <- rNVASIM(80, mu = 0.55, sigma = 0.25)
    y_quant <- rNVASIQ(80, mu = 0.55, sigma = 0.25, tau = 0.25)
    control <- gamlss::gamlss.control(n.cyc = 50, trace = FALSE)

    fit_mean <- gamlss::gamlss(
        y_mean ~ 1, family = NVASIM(), control = control
    )
    fit_quant <- with_global_tau(
        0.25,
        gamlss::gamlss(
            y_quant ~ 1, family = NVASIQ(), control = control
        )
    )

    expect_s3_class(fit_mean, "gamlss")
    expect_s3_class(fit_quant, "gamlss")
})

test_that("LVASIQ fits an intercept-only GAMLSS model", {
    set.seed(321)
    tau <- 0.25
    y <- rLVASIQ(100, mu = 0.55, sigma = 0.25, tau = tau)
    control <- gamlss::gamlss.control(n.cyc = 75, trace = FALSE)

    fit <- with_global_tau(
        tau,
        gamlss::gamlss(
            y ~ 1,
            family = LVASIQ(),
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
