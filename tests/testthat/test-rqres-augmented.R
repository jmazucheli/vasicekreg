evaluate_augmented_rqres <- function(family, y, mu, sigma, nu, tau = NULL) {
    evaluation_environment <- list2env(
        list(y = y, mu = mu, sigma = sigma, nu = nu),
        parent = parent.frame()
    )

    if (!is.null(tau)) {
        assign("tau", tau, envir = evaluation_environment)
    }

    eval(family$rqres, envir = evaluation_environment)
}

test_that("ZANVASIM randomized quantile residuals use the CDF jump at zero", {
    y <- c(0, 0.23, 0, 0.71)
    mu <- c(0.42, 0.48, 0.57, 0.63)
    sigma <- c(0.21, 0.26, 0.31, 0.36)
    nu <- c(0.12, 0.18, 0.27, 0.34)

    set.seed(20260813)
    residuals <- evaluate_augmented_rqres(
        ZANVASIM(), y, mu, sigma, nu
    )
    randomized_probabilities <- pnorm(residuals)
    at_zero <- y == 0

    expect_true(all(randomized_probabilities[at_zero] >= 0))
    expect_true(all(randomized_probabilities[at_zero] <= nu[at_zero]))
    expect_equal(
        randomized_probabilities[!at_zero],
        p0NVASIM(
            y[!at_zero], mu[!at_zero], sigma[!at_zero], nu[!at_zero]
        ),
        tolerance = 1e-12
    )

    set.seed(20260813)
    repeated_residuals <- evaluate_augmented_rqres(
        ZANVASIM(), y, mu, sigma, nu
    )
    expect_identical(residuals, repeated_residuals)
})

test_that("OANVASIM randomized quantile residuals use the CDF jump at one", {
    y <- c(0.19, 1, 0.56, 1)
    mu <- c(0.42, 0.48, 0.57, 0.63)
    sigma <- c(0.21, 0.26, 0.31, 0.36)
    nu <- c(0.12, 0.18, 0.27, 0.34)

    set.seed(20260813)
    residuals <- evaluate_augmented_rqres(
        OANVASIM(), y, mu, sigma, nu
    )
    randomized_probabilities <- pnorm(residuals)
    at_one <- y == 1

    expect_true(
        all(randomized_probabilities[at_one] >= 1 - nu[at_one])
    )
    expect_true(all(randomized_probabilities[at_one] <= 1))
    expect_equal(
        randomized_probabilities[!at_one],
        p1NVASIM(
            y[!at_one], mu[!at_one], sigma[!at_one], nu[!at_one]
        ),
        tolerance = 1e-12
    )

    set.seed(20260813)
    repeated_residuals <- evaluate_augmented_rqres(
        OANVASIM(), y, mu, sigma, nu
    )
    expect_identical(residuals, repeated_residuals)
})

test_that("ZOANVASIM residuals use both CDF jumps and the interior CDF", {
    y <- c(0, 0.23, 1, 0.71, 0, 1)
    mu <- c(0.42, 0.48, 0.52, 0.57, 0.61, 0.66)
    sigma <- c(0.21, 0.24, 0.27, 0.30, 0.33, 0.36)
    nu <- c(0.12, 0.16, 0.20, 0.24, 0.28, 0.32)
    tau <- c(0.11, 0.15, 0.19, 0.23, 0.27, 0.31)

    set.seed(20260813)
    residuals <- evaluate_augmented_rqres(
        ZOANVASIM(), y, mu, sigma, nu, tau
    )
    randomized_probabilities <- pnorm(residuals)
    at_zero <- y == 0
    at_one <- y == 1
    interior <- y > 0 & y < 1
    left_limit_at_one <- nu + (1 - nu) * (1 - tau)

    expect_true(all(randomized_probabilities[at_zero] >= 0))
    expect_true(all(randomized_probabilities[at_zero] <= nu[at_zero]))
    expect_true(
        all(
            randomized_probabilities[at_one] >=
                left_limit_at_one[at_one]
        )
    )
    expect_true(all(randomized_probabilities[at_one] <= 1))
    expect_equal(
        randomized_probabilities[interior],
        p01NVASIM(
            y[interior], mu[interior], sigma[interior],
            nu[interior], tau[interior]
        ),
        tolerance = 1e-12
    )

    set.seed(20260813)
    repeated_residuals <- evaluate_augmented_rqres(
        ZOANVASIM(), y, mu, sigma, nu, tau
    )
    expect_identical(residuals, repeated_residuals)
})

test_that("augmented-family quantile residuals are approximately normal", {
    set.seed(20260813)
    sample_size <- 10000L
    mu <- 0.58
    sigma <- 0.30
    nu <- 0.22
    tau <- 0.27

    generated_samples <- list(
        ZANVASIM = r0NVASIM(sample_size, mu, sigma, nu),
        OANVASIM = r1NVASIM(sample_size, mu, sigma, nu),
        ZOANVASIM = r01NVASIM(sample_size, mu, sigma, nu, tau)
    )
    families <- list(
        ZANVASIM = ZANVASIM(),
        OANVASIM = OANVASIM(),
        ZOANVASIM = ZOANVASIM()
    )

    residuals <- lapply(names(generated_samples), function(family_name) {
        family_tau <- if (family_name == "ZOANVASIM") tau else NULL
        evaluate_augmented_rqres(
            families[[family_name]],
            generated_samples[[family_name]],
            mu,
            sigma,
            nu,
            family_tau
        )
    })

    for (family_residuals in residuals) {
        expect_true(all(is.finite(family_residuals)))
        expect_equal(mean(family_residuals), 0, tolerance = 0.04)
        expect_equal(stats::sd(family_residuals), 1, tolerance = 0.04)
        expect_equal(
            as.numeric(stats::quantile(
                family_residuals, c(0.025, 0.50, 0.975)
            )),
            stats::qnorm(c(0.025, 0.50, 0.975)),
            tolerance = 0.10
        )
    }
})
