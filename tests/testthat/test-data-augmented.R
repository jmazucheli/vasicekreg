test_that("transport data retain their documented structure", {
    data("transport", package = "vasicekreg", envir = environment())

    expect_s3_class(transport, "data.frame")
    expect_equal(dim(transport), c(60L, 7L))
    expect_named(
        transport,
        c(
            "ntrips", "nbiked", "status", "gender", "parking",
            "distance", "propbiked"
        )
    )
    expect_equal(transport$propbiked, transport$nbiked / transport$ntrips)
    expect_equal(sum(transport$propbiked == 0), 24L)
    expect_equal(sum(transport$propbiked == 1), 0L)
})

test_that("trees data retain their documented structure", {
    data("trees", package = "vasicekreg", envir = environment())

    expect_s3_class(trees, "data.frame")
    expect_equal(dim(trees), c(26L, 7L))
    expect_named(
        trees,
        c(
            "planted", "survived", "pest", "fertilization", "precip",
            "wind", "prop"
        )
    )
    expect_equal(trees$prop, trees$survived / trees$planted)
    expect_equal(sum(trees$prop == 0), 0L)
    expect_equal(sum(trees$prop == 1), 6L)
})

test_that("aep data retain their documented structure", {
    data("aep", package = "vasicekreg", envir = environment())

    expect_s3_class(aep, "data.frame")
    expect_equal(dim(aep), c(1383L, 8L))
    expect_named(
        aep,
        c("los", "noinap", "loglos", "sex", "ward", "year", "age", "y")
    )
    expect_equal(aep$loglos, log(aep$los / 10))
    expect_equal(unname(aep$y[, 1L]), aep$noinap)
    expect_equal(unname(aep$y[, 2L]), aep$los - aep$noinap)

    inappropriate <- aep$noinap / aep$los
    expect_equal(sum(inappropriate == 0), 763L)
    expect_equal(sum(inappropriate > 0 & inappropriate < 1), 552L)
    expect_equal(sum(inappropriate == 1), 68L)
})
