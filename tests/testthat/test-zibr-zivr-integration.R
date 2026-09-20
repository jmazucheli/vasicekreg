test_that("longitudinal two-part functions are exported", {
    expect_true(is.function(zibr))
    expect_true(is.function(zivr))

    exports <- getNamespaceExports("vasicekreg")
    expect_true(all(c("zibr", "zivr") %in% exports))
})

test_that("the shared longitudinal engine is available internally", {
    engine <- getFromNamespace(
        ".fit_zero_inflated_engine",
        "vasicekreg"
    )
    expect_true(is.function(engine))
})

test_that("zibr S3 methods are registered", {
    generics <- c("print", "coef", "vcov", "logLik", "nobs", "BIC")

    for (generic in generics) {
        expect_true(is.function(
            utils::getS3method(generic, "zibr", optional = TRUE)
        ))
    }
})

test_that("zivr S3 methods are registered", {
    generics <- c("print", "coef", "vcov", "logLik", "nobs", "BIC")

    for (generic in generics) {
        expect_true(is.function(
            utils::getS3method(generic, "zivr", optional = TRUE)
        ))
    }
})

test_that("formula and legacy arguments remain present", {
    expect_true(all(
        c(
            "data", "y", "formula_bin", "formula_cont", "random",
            "logistic_cov", "beta_cov", "subject_ind", "time_ind"
        ) %in% names(formals(zibr))
    ))

    expect_true(all(
        c(
            "data", "y", "formula_bin", "formula_cont", "random",
            "logistic_cov", "vasicek_cov", "subject_ind", "time_ind"
        ) %in% names(formals(zivr))
    ))
})
