first_derivative <- function(f, x, h = 1e-6) {
    (f(x + h) - f(x - h)) / (2 * h)
}

second_derivative <- function(f, x, h = 1e-4) {
    (f(x + h) - 2 * f(x) + f(x - h)) / h^2
}

cross_derivative <- function(f, x, y, h = 1e-4) {
    (f(x + h, y + h) - f(x + h, y - h) -
         f(x - h, y + h) + f(x - h, y - h)) / (4 * h^2)
}

test_that("NVASIM analytical derivatives agree with numerical derivatives", {
    y <- 0.37
    mu <- 0.69
    sigma <- 0.37
    family <- NVASIM()
    loglik <- function(m, s) dNVASIM(y, m, s, log = TRUE)

    expect_equal(
        family$dldm(y, mu, sigma),
        first_derivative(function(m) loglik(m, sigma), mu),
        tolerance = 1e-5
    )
    expect_equal(
        family$d2ldm2(y, mu, sigma),
        second_derivative(function(m) loglik(m, sigma), mu),
        tolerance = 1e-4
    )
    expect_equal(
        family$dldd(y, mu, sigma),
        first_derivative(function(s) loglik(mu, s), sigma),
        tolerance = 1e-5
    )
    expect_equal(
        family$d2ldd2(y, mu, sigma),
        second_derivative(function(s) loglik(mu, s), sigma),
        tolerance = 1e-4
    )
    expect_equal(
        family$d2ldmdd(y, mu, sigma),
        cross_derivative(loglik, mu, sigma),
        tolerance = 1e-4
    )
})

test_that("NVASIQ analytical derivatives agree with numerical derivatives", {
    y <- 0.37
    mu <- 0.69
    sigma <- 0.37
    tau <- 0.81
    loglik <- function(m, s) dNVASIQ(y, m, s, tau = tau, log = TRUE)

    with_global_tau(tau, {
        family <- NVASIQ()

        expect_equal(
            family$dldm(y, mu, sigma),
            first_derivative(function(m) loglik(m, sigma), mu),
            tolerance = 1e-5
        )
        expect_equal(
            family$d2ldm2(y, mu, sigma),
            second_derivative(function(m) loglik(m, sigma), mu),
            tolerance = 1e-4
        )
        expect_equal(
            family$dldd(y, mu, sigma),
            first_derivative(function(s) loglik(mu, s), sigma),
            tolerance = 1e-5
        )
        expect_equal(
            family$d2ldd2(y, mu, sigma),
            second_derivative(function(s) loglik(mu, s), sigma),
            tolerance = 1e-4
        )
        expect_equal(
            family$d2ldmdd(y, mu, sigma),
            cross_derivative(loglik, mu, sigma),
            tolerance = 1e-4
        )
    })
})

test_that("LVASIQ analytical derivatives agree with numerical derivatives", {
    y <- 0.37
    mu <- 0.69
    sigma <- 0.37
    tau <- 0.81
    loglik <- function(m, s) dLVASIQ(y, m, s, tau = tau, log = TRUE)

    with_global_tau(tau, {
        family <- LVASIQ()

        expect_equal(
            family$dldm(y, mu, sigma),
            first_derivative(function(m) loglik(m, sigma), mu),
            tolerance = 1e-5
        )
        expect_equal(
            family$d2ldm2(y, mu, sigma),
            second_derivative(function(m) loglik(m, sigma), mu),
            tolerance = 1e-4
        )
        expect_equal(
            family$dldd(y, mu, sigma),
            first_derivative(function(s) loglik(mu, s), sigma),
            tolerance = 1e-5
        )
        expect_equal(
            family$d2ldd2(y, mu, sigma),
            second_derivative(function(s) loglik(mu, s), sigma),
            tolerance = 1e-4
        )
        expect_equal(
            family$d2ldmdd(y, mu, sigma),
            cross_derivative(loglik, mu, sigma),
            tolerance = 1e-4
        )
    })
})
