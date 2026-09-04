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
    quantile <- 0.81
    loglik <- function(m, s) {
        dNVASIQ(y, m, s, quantile = quantile, log = TRUE)
    }
    family <- NVASIQ(quantile = quantile)

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

test_that("LVASIQ analytical derivatives agree with numerical derivatives", {
    y <- 0.37
    mu <- 0.69
    sigma <- 0.37
    quantile <- 0.81
    loglik <- function(m, s) {
        dLVASIQ(y, m, s, quantile = quantile, log = TRUE)
    }
    family <- LVASIQ(quantile = quantile)

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

test_that("ZANVASIM analytical derivatives agree with its log-likelihood", {
    y <- 0.37
    mu <- 0.69
    sigma <- 0.37
    nu <- 0.23
    responses <- c(0, y)
    probabilities <- c(nu, 1 - nu)
    family <- ZANVASIM()
    loglik <- function(response, m = mu, s = sigma, v = nu) {
        d0NVASIM(response, m, s, v, log = TRUE)
    }

    expect_equal(
        family$dldm(y, mu, sigma),
        first_derivative(function(m) loglik(y, m = m), mu),
        tolerance = 1e-5
    )
    expect_equal(
        family$d2ldm2(y, mu, sigma),
        second_derivative(function(m) loglik(y, m = m), mu),
        tolerance = 1e-4
    )
    expect_equal(
        family$dldd(y, mu, sigma),
        first_derivative(function(s) loglik(y, s = s), sigma),
        tolerance = 1e-5
    )
    expect_equal(
        family$d2ldd2(y, mu, sigma),
        second_derivative(function(s) loglik(y, s = s), sigma),
        tolerance = 1e-4
    )
    expect_equal(
        family$d2ldmdd(y, mu, sigma),
        cross_derivative(function(m, s) loglik(y, m, s), mu, sigma),
        tolerance = 1e-4
    )

    numerical_dldv <- vapply(
        responses,
        function(response) {
            first_derivative(function(v) loglik(response, v = v), nu)
        },
        numeric(1)
    )
    expect_equal(
        family$dldv(responses, nu), numerical_dldv, tolerance = 1e-5
    )

    expected_d2ldv2 <- sum(
        probabilities * vapply(
            responses,
            function(response) {
                second_derivative(function(v) loglik(response, v = v), nu)
            },
            numeric(1)
        )
    )
    expect_equal(family$d2ldv2(nu), expected_d2ldv2, tolerance = 1e-4)

    for (response in responses) {
        expect_equal(
            family$d2ldmdv(response),
            cross_derivative(
                function(m, v) loglik(response, m = m, v = v), mu, nu
            ),
            tolerance = 1e-4
        )
        expect_equal(
            family$d2ldddv(response),
            cross_derivative(
                function(s, v) loglik(response, s = s, v = v), sigma, nu
            ),
            tolerance = 1e-4
        )
    }

    expect_equal(family$dldm(0, mu, sigma), 0)
    expect_equal(family$d2ldm2(0, mu, sigma), 0)
    expect_equal(family$dldd(0, mu, sigma), 0)
    expect_equal(family$d2ldd2(0, mu, sigma), 0)
    expect_equal(family$d2ldmdd(0, mu, sigma), 0)
})

test_that("OANVASIM analytical derivatives agree with its log-likelihood", {
    y <- 0.37
    mu <- 0.69
    sigma <- 0.37
    nu <- 0.23
    responses <- c(y, 1)
    probabilities <- c(1 - nu, nu)
    family <- OANVASIM()
    loglik <- function(response, m = mu, s = sigma, v = nu) {
        d1NVASIM(response, m, s, v, log = TRUE)
    }

    expect_equal(
        family$dldm(y, mu, sigma),
        first_derivative(function(m) loglik(y, m = m), mu),
        tolerance = 1e-5
    )
    expect_equal(
        family$d2ldm2(y, mu, sigma),
        second_derivative(function(m) loglik(y, m = m), mu),
        tolerance = 1e-4
    )
    expect_equal(
        family$dldd(y, mu, sigma),
        first_derivative(function(s) loglik(y, s = s), sigma),
        tolerance = 1e-5
    )
    expect_equal(
        family$d2ldd2(y, mu, sigma),
        second_derivative(function(s) loglik(y, s = s), sigma),
        tolerance = 1e-4
    )
    expect_equal(
        family$d2ldmdd(y, mu, sigma),
        cross_derivative(function(m, s) loglik(y, m, s), mu, sigma),
        tolerance = 1e-4
    )

    numerical_dldv <- vapply(
        responses,
        function(response) {
            first_derivative(function(v) loglik(response, v = v), nu)
        },
        numeric(1)
    )
    expect_equal(
        family$dldv(responses, nu), numerical_dldv, tolerance = 1e-5
    )

    expected_d2ldv2 <- sum(
        probabilities * vapply(
            responses,
            function(response) {
                second_derivative(function(v) loglik(response, v = v), nu)
            },
            numeric(1)
        )
    )
    expect_equal(family$d2ldv2(nu), expected_d2ldv2, tolerance = 1e-4)

    for (response in responses) {
        expect_equal(
            family$d2ldmdv(response),
            cross_derivative(
                function(m, v) loglik(response, m = m, v = v), mu, nu
            ),
            tolerance = 1e-4
        )
        expect_equal(
            family$d2ldddv(response),
            cross_derivative(
                function(s, v) loglik(response, s = s, v = v), sigma, nu
            ),
            tolerance = 1e-4
        )
    }

    expect_equal(family$dldm(1, mu, sigma), 0)
    expect_equal(family$d2ldm2(1, mu, sigma), 0)
    expect_equal(family$dldd(1, mu, sigma), 0)
    expect_equal(family$d2ldd2(1, mu, sigma), 0)
    expect_equal(family$d2ldmdd(1, mu, sigma), 0)
})

test_that("ZOANVASIM analytical derivatives agree with its log-likelihood", {
    y <- 0.37
    mu <- 0.69
    sigma <- 0.37
    nu <- 0.35
    tau <- 0.55
    responses <- c(0, y, 1)
    probabilities <- c(nu, (1 - nu) * (1 - tau), (1 - nu) * tau)
    family <- ZOANVASIM()
    loglik <- function(response, m = mu, s = sigma,
                       v = nu, t = tau) {
        d01NVASIM(response, m, s, v, t, log = TRUE)
    }

    expect_equal(
        family$dldm(y, mu, sigma),
        first_derivative(function(m) loglik(y, m = m), mu),
        tolerance = 1e-5
    )
    expect_equal(
        family$d2ldm2(y, mu, sigma),
        second_derivative(function(m) loglik(y, m = m), mu),
        tolerance = 1e-4
    )
    expect_equal(
        family$dldd(y, mu, sigma),
        first_derivative(function(s) loglik(y, s = s), sigma),
        tolerance = 1e-5
    )
    expect_equal(
        family$d2ldd2(y, mu, sigma),
        second_derivative(function(s) loglik(y, s = s), sigma),
        tolerance = 1e-4
    )
    expect_equal(
        family$d2ldmdd(y, mu, sigma),
        cross_derivative(function(m, s) loglik(y, m, s), mu, sigma),
        tolerance = 1e-4
    )

    numerical_dldv <- vapply(
        responses,
        function(response) {
            first_derivative(function(v) loglik(response, v = v), nu)
        },
        numeric(1)
    )
    numerical_dldt <- vapply(
        responses,
        function(response) {
            first_derivative(function(t) loglik(response, t = t), tau)
        },
        numeric(1)
    )
    expect_equal(
        family$dldv(responses, nu), numerical_dldv, tolerance = 1e-5
    )
    expect_equal(
        family$dldt(responses, nu, tau), numerical_dldt, tolerance = 1e-5
    )

    expected_d2ldv2 <- sum(
        probabilities * vapply(
            responses,
            function(response) {
                second_derivative(function(v) loglik(response, v = v), nu)
            },
            numeric(1)
        )
    )
    expected_d2ldt2 <- sum(
        probabilities * vapply(
            responses,
            function(response) {
                second_derivative(function(t) loglik(response, t = t), tau)
            },
            numeric(1)
        )
    )
    expect_equal(family$d2ldv2(nu), expected_d2ldv2, tolerance = 1e-4)
    expect_equal(
        family$d2ldt2(nu, tau), expected_d2ldt2, tolerance = 1e-4
    )

    for (response in responses) {
        expect_equal(
            family$d2ldmdv(response),
            cross_derivative(
                function(m, v) loglik(response, m = m, v = v), mu, nu
            ),
            tolerance = 1e-4
        )
        expect_equal(
            family$d2ldmdt(response),
            cross_derivative(
                function(m, t) loglik(response, m = m, t = t), mu, tau
            ),
            tolerance = 1e-4
        )
        expect_equal(
            family$d2ldddv(response),
            cross_derivative(
                function(s, v) loglik(response, s = s, v = v), sigma, nu
            ),
            tolerance = 1e-4
        )
        expect_equal(
            family$d2ldddt(response),
            cross_derivative(
                function(s, t) loglik(response, s = s, t = t), sigma, tau
            ),
            tolerance = 1e-4
        )
        expect_equal(
            family$d2ldvdt(nu, tau),
            cross_derivative(
                function(v, t) loglik(response, v = v, t = t), nu, tau
            ),
            tolerance = 1e-4
        )
    }

    expect_equal(family$dldm(c(0, 1), mu, sigma), c(0, 0))
    expect_equal(family$d2ldm2(c(0, 1), mu, sigma), c(0, 0))
    expect_equal(family$dldd(c(0, 1), mu, sigma), c(0, 0))
    expect_equal(family$d2ldd2(c(0, 1), mu, sigma), c(0, 0))
    expect_equal(family$d2ldmdd(c(0, 1), mu, sigma), c(0, 0))
})
