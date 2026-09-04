contains_symbol <- function(expression, name) {
    name %in% all.vars(expression, functions = FALSE)
}

test_that("quantile families are independent of global hyperparameters", {
    global_names <- c("tau", "quantile")
    existed <- vapply(
        global_names, exists, logical(1),
        envir = .GlobalEnv, inherits = FALSE
    )
    previous <- lapply(global_names[existed], get, envir = .GlobalEnv)
    names(previous) <- global_names[existed]
    on.exit({
        for (name in global_names) {
            if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
                rm(list = name, envir = .GlobalEnv)
            }
        }
        for (name in names(previous)) {
            assign(name, previous[[name]], envir = .GlobalEnv)
        }
    }, add = TRUE)

    assign("tau", 0.99, envir = .GlobalEnv)
    assign("quantile", 0.01, envir = .GlobalEnv)

    specifications <- list(
        NVASIQ = list(factory = NVASIQ, density = dNVASIQ),
        LVASIQ = list(factory = LVASIQ, density = dLVASIQ),
        HVASIQ = list(factory = HVASIQ, density = dHVASIQ)
    )
    for (specification in specifications) {
        lower <- specification$factory(quantile = 0.25)
        upper <- specification$factory(quantile = 0.75)

        expect_identical(lower$quantile, 0.25)
        expect_identical(upper$quantile, 0.75)
        expect_identical(lower$rqres[[1L]][["quantile"]], 0.25)
        expect_identical(upper$rqres[[1L]][["quantile"]], 0.75)

        evaluation_data <- list(y = 0.40, mu = 0.55, sigma = 0.30, w = 1)
        expect_equal(
            eval(
                body(lower$G.dev.incr),
                envir = evaluation_data,
                enclos = asNamespace("vasicekreg")
            ),
            -2 * specification$density(
                0.40, 0.55, 0.30, quantile = 0.25, log = TRUE
            )
        )
    }
})

test_that("GAMLSS family components contain no free quantile symbol", {
    component_names <- c(
        "dldm", "d2ldm2", "dldd", "d2ldmdd", "d2ldd2",
        "G.dev.incr", "mean", "variance"
    )

    for (factory in list(NVASIQ, LVASIQ, HVASIQ)) {
        family <- factory(quantile = 0.37)

        for (component_name in component_names) {
            expect_false(
                contains_symbol(body(family[[component_name]]), "quantile"),
                info = paste(family$family[1L], component_name)
            )
        }
        expect_false(
            contains_symbol(family$rqres, "quantile"),
            info = paste(family$family[1L], "rqres")
        )
        expect_false(
            contains_symbol(family$mu.initial, "quantile"),
            info = paste(family$family[1L], "mu.initial")
        )
    }
})

test_that("fixed quantile can be recovered after serialization", {
    family <- LVASIQ(quantile = 0.35)
    object <- list(rqres = family$rqres)

    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path), add = TRUE)
    saveRDS(object, path)
    restored <- readRDS(path)

    expect_identical(.envelope_fixed_quantile(restored), 0.35)
})

test_that("GAMLSS fits at different quantiles coexist safely", {
    skip_if_not_installed("gamlss")
    set.seed(20260903)
    data <- data.frame(
        y = rNVASIQ(
            80, mu = 0.55, sigma = 0.25, quantile = 0.25
        )
    )
    control <- gamlss::gamlss.control(n.cyc = 75, trace = FALSE)

    fit_025 <- gamlss::gamlss(
        y ~ 1, data = data,
        family = NVASIQ(quantile = 0.25), control = control
    )
    fit_075 <- gamlss::gamlss(
        y ~ 1, data = data,
        family = NVASIQ(quantile = 0.75), control = control
    )

    expect_identical(.envelope_fixed_quantile(fit_025), 0.25)
    expect_identical(.envelope_fixed_quantile(fit_075), 0.75)
    expect_equal(
        fit_025$residuals,
        stats::qnorm(pNVASIQ(
            data$y,
            fitted(fit_025, what = "mu"),
            fitted(fit_025, what = "sigma"),
            quantile = 0.25
        )),
        tolerance = 1e-10
    )
    expect_equal(
        fit_075$residuals,
        stats::qnorm(pNVASIQ(
            data$y,
            fitted(fit_075, what = "mu"),
            fitted(fit_075, what = "sigma"),
            quantile = 0.75
        )),
        tolerance = 1e-10
    )
})

test_that("automatic envelope refitting embeds a local quantile value", {
    skip_if_not_installed("gamlss")
    set.seed(20260903)
    data <- data.frame(
        y = rNVASIQ(
            60, mu = 0.55, sigma = 0.25, quantile = 0.25
        )
    )
    control <- gamlss::gamlss.control(n.cyc = 75, trace = FALSE)

    local_fit <- function(quantile_level) {
        local_formula <- y ~ 1
        gamlss::gamlss(
            formula = local_formula,
            sigma.formula = ~ 1,
            family = NVASIQ(quantile_level),
            data = data,
            control = control
        )
    }
    fit <- local_fit(0.25)

    set.seed(91)
    simulated_y <- .envelope_simulate_response(fit)
    refitted <- .envelope_refit(fit, simulated_y, data)

    expect_s3_class(refitted, "gamlss")
    expect_identical(.envelope_fixed_quantile(refitted), 0.25)
})
