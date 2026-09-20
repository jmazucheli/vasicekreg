test_that("longitudinal two-part functions are exported", {
    expect_true(is.function(zab))
    expect_true(is.function(zavr))

    exports <- getNamespaceExports("vasicekreg")
    expect_true(all(c("zab", "zavr") %in% exports))
})

test_that("the shared longitudinal engine is available internally", {
    engine <- getFromNamespace(
        ".fit_zero_augmented_engine",
        "vasicekreg"
    )
    expect_true(is.function(engine))
})

test_that("zab S3 methods are registered", {
    generics <- c("print", "coef", "vcov", "logLik", "nobs", "BIC")

    for (generic in generics) {
        expect_true(is.function(
            utils::getS3method(generic, "zab", optional = TRUE)
        ))
    }
})

test_that("zavr S3 methods are registered", {
    generics <- c("print", "coef", "vcov", "logLik", "nobs", "BIC")

    for (generic in generics) {
        expect_true(is.function(
            utils::getS3method(generic, "zavr", optional = TRUE)
        ))
    }
})

test_that("formula and legacy arguments remain present", {
    expect_true(all(
        c(
            "data", "y", "formula_bin", "formula_cont", "random",
            "logistic_cov", "beta_cov", "subject_ind", "time_ind"
        ) %in% names(formals(zab))
    ))

    expect_true(all(
        c(
            "data", "y", "formula_bin", "formula_cont", "random",
            "logistic_cov", "vasicek_cov", "subject_ind", "time_ind"
        ) %in% names(formals(zavr))
    ))
})
