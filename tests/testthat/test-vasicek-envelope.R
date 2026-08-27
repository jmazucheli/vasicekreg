test_that("vasicek_envelope returns coherent quantile and Cox-Snell results", {
    skip_if_not_installed("gamlss")
    set.seed(20260827)
    dat <- data.frame(y = rNVASIM(80, mu = 0.58, sigma = 0.30))
    fit <- gamlss::gamlss(
        y ~ 1, sigma.formula = ~ 1, family = NVASIM(), data = dat,
        control = gamlss::gamlss.control(trace = FALSE)
    )

    envelope <- vasicek_envelope(
        fit, nsim = 5, seed = 20260827, data = dat, verbose = FALSE
    )

    expect_s3_class(envelope, "vasicek_envelope")
    expect_equal(envelope$nsim, 5L)
    expect_equal(envelope$failures, envelope$attempts - envelope$nsim)
    expect_named(envelope$results, c("quantile", "cox-snell"))
    expect_equal(dim(envelope$results$quantile$simulated), c(80L, 5L))
    expect_true(all(
        envelope$results$quantile$lower <=
            envelope$results$quantile$upper
    ))
    expect_equal(
        envelope$results[["cox-snell"]]$observed,
        sort(-pnorm(
            envelope$results$quantile$observed,
            lower.tail = FALSE,
            log.p = TRUE
        )),
        tolerance = 1e-12
    )
})

test_that("minmax envelope contains every simulated ordered residual", {
    skip_if_not_installed("gamlss")
    set.seed(20260827)
    dat <- data.frame(y = rNVASIM(60, mu = 0.52, sigma = 0.35))
    fit <- gamlss::gamlss(
        y ~ 1, sigma.formula = ~ 1, family = NVASIM(), data = dat,
        control = gamlss::gamlss.control(trace = FALSE)
    )
    envelope <- vasicek_envelope(
        fit, residual = "quantile", nsim = 4, envelope = "minmax",
        seed = 20260827, data = dat, verbose = FALSE
    )
    result <- envelope$results$quantile
    expect_equal(result$lower, apply(result$simulated, 1L, min))
    expect_equal(result$upper, apply(result$simulated, 1L, max))
})

test_that("seed is reproducible without changing the caller RNG state", {
    skip_if_not_installed("gamlss")
    set.seed(20260827)
    dat <- data.frame(y = rNVASIM(50, mu = 0.55, sigma = 0.28))
    fit <- gamlss::gamlss(
        y ~ 1, sigma.formula = ~ 1, family = NVASIM(), data = dat,
        control = gamlss::gamlss.control(trace = FALSE)
    )
    state_before <- .Random.seed
    first <- vasicek_envelope(
        fit, residual = "quantile", nsim = 3, seed = 17,
        data = dat, verbose = FALSE
    )
    expect_identical(.Random.seed, state_before)
    second <- vasicek_envelope(
        fit, residual = "quantile", nsim = 3, seed = 17,
        data = dat, verbose = FALSE
    )
    expect_identical(
        first$results$quantile$simulated,
        second$results$quantile$simulated
    )
})


test_that("Cox-Snell transformation remains finite in the upper tail", {
    residual <- 10
    cox_snell <- -pnorm(residual, lower.tail = FALSE, log.p = TRUE)

    expect_true(is.finite(cox_snell))
    expect_gt(cox_snell, 50)
})

test_that("automatic refitting rejects data with incompatible row count", {
    skip_if_not_installed("gamlss")
    set.seed(20260827)
    dat <- data.frame(y = rNVASIM(30, mu = 0.55, sigma = 0.30))
    fit <- gamlss::gamlss(
        y ~ 1, sigma.formula = ~ 1, family = NVASIM(), data = dat,
        control = gamlss::gamlss.control(trace = FALSE)
    )

    expect_error(
        .envelope_refit(fit, y = numeric(31), data = dat),
        "has 30 rows but the fitted model used 31 observations",
        fixed = TRUE
    )
})
