grid_first_derivative <- function(f, x, h = 1e-6) {
    (f(x + h) - f(x - h)) / (2 * h)
}

grid_second_derivative <- function(f, x, h = 1e-4) {
    (f(x + h) - 2 * f(x) + f(x - h)) / h^2
}

grid_cross_derivative <- function(f, x, y, h = 1e-4) {
    (f(x + h, y + h) - f(x + h, y - h) -
         f(x - h, y + h) + f(x - h, y - h)) / (4 * h^2)
}

test_that("quantile-family analytical derivatives agree on a parameter grid", {
    parameter_grid <- expand.grid(
        y = c(0.20, 0.70),
        mu = c(0.35, 0.65),
        sigma = c(0.25, 0.60),
        quantile = c(0.20, 0.75)
    )
    specifications <- list(
        NVASIQ = list(family = NVASIQ, density = dNVASIQ),
        LVASIQ = list(family = LVASIQ, density = dLVASIQ),
        HSVASIQ = list(family = HSVASIQ, density = dHSVASIQ)
    )

    for (specification in specifications) {
        for (row in seq_len(nrow(parameter_grid))) {
            values <- parameter_grid[row, ]
            family <- specification$family(quantile = values$quantile)
            loglik <- function(mu, sigma) {
                specification$density(
                    values$y, mu, sigma,
                    quantile = values$quantile, log = TRUE
                )
            }
            expect_equal(
                family$dldm(values$y, values$mu, values$sigma),
                grid_first_derivative(
                    function(mu) loglik(mu, values$sigma), values$mu
                ),
                tolerance = 2e-5
            )
            expect_equal(
                family$dldd(values$y, values$mu, values$sigma),
                grid_first_derivative(
                    function(sigma) loglik(values$mu, sigma), values$sigma
                ),
                tolerance = 2e-5
            )
            expect_equal(
                family$d2ldm2(values$y, values$mu, values$sigma),
                grid_second_derivative(
                    function(mu) loglik(mu, values$sigma), values$mu
                ),
                tolerance = 3e-4
            )
            expect_equal(
                family$d2ldd2(values$y, values$mu, values$sigma),
                grid_second_derivative(
                    function(sigma) loglik(values$mu, sigma), values$sigma
                ),
                tolerance = 3e-4
            )
            expect_equal(
                family$d2ldmdd(values$y, values$mu, values$sigma),
                grid_cross_derivative(loglik, values$mu, values$sigma),
                tolerance = 3e-4
            )
        }
    }
})
