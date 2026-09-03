test_that("zero-augmented dpqr functions follow the mixture definition", {
    mu <- 0.61
    sigma <- 0.28
    nu <- 0.17
    x <- c(0.12, 0.44, 0.87)
    component_probability <- c(0.10, 0.50, 0.90)
    mixture_probability <- nu + (1 - nu) * component_probability

    expect_equal(d0NVASIM(0, mu, sigma, nu), nu)
    expect_equal(
        d0NVASIM(x, mu, sigma, nu),
        (1 - nu) * dNVASIM(x, mu, sigma)
    )
    expect_equal(p0NVASIM(0, mu, sigma, nu), nu)
    expect_equal(
        p0NVASIM(x, mu, sigma, nu),
        nu + (1 - nu) * pNVASIM(x, mu, sigma)
    )
    expect_equal(
        q0NVASIM(mixture_probability, mu, sigma, nu),
        qNVASIM(component_probability, mu, sigma),
        tolerance = 1e-12
    )
    expect_equal(q0NVASIM(c(0, nu), mu, sigma, nu), c(0, 0))
    expect_equal(
        d0NVASIM(0, mu, sigma, nu) +
            integrate(
                function(value) d0NVASIM(value, mu, sigma, nu),
                lower = 0,
                upper = 1,
                subdivisions = 500L,
                rel.tol = 1e-10
            )$value,
        1,
        tolerance = 1e-7
    )
})

test_that("one-augmented dpqr functions follow the mixture definition", {
    mu <- 0.61
    sigma <- 0.28
    nu <- 0.17
    x <- c(0.12, 0.44, 0.87)
    component_probability <- c(0.10, 0.50, 0.90)
    mixture_probability <- (1 - nu) * component_probability

    expect_equal(d1NVASIM(1, mu, sigma, nu), nu)
    expect_equal(
        d1NVASIM(x, mu, sigma, nu),
        (1 - nu) * dNVASIM(x, mu, sigma)
    )
    expect_equal(p1NVASIM(c(0, 1), mu, sigma, nu), c(0, 1))
    expect_equal(
        p1NVASIM(x, mu, sigma, nu),
        (1 - nu) * pNVASIM(x, mu, sigma)
    )
    expect_equal(
        q1NVASIM(mixture_probability, mu, sigma, nu),
        qNVASIM(component_probability, mu, sigma),
        tolerance = 1e-12
    )
    expect_equal(q1NVASIM(c(1 - nu, 1), mu, sigma, nu), c(1, 1))
    expect_equal(
        d1NVASIM(1, mu, sigma, nu) +
            integrate(
                function(value) d1NVASIM(value, mu, sigma, nu),
                lower = 0,
                upper = 1,
                subdivisions = 500L,
                rel.tol = 1e-10
            )$value,
        1,
        tolerance = 1e-7
    )
})

test_that("zero-and-one-augmented dpqr functions follow the mixture definition", {
    mu <- 0.61
    sigma <- 0.28
    nu <- 0.17
    tau <- 0.24
    continuous_probability <- (1 - nu) * (1 - tau)
    one_probability <- (1 - nu) * tau
    x <- c(0.12, 0.44, 0.87)
    component_probability <- c(0.10, 0.50, 0.90)
    mixture_probability <- nu + continuous_probability * component_probability

    expect_equal(d01NVASIM(0, mu, sigma, nu, tau), nu)
    expect_equal(d01NVASIM(1, mu, sigma, nu, tau), one_probability)
    expect_equal(p01NVASIM(c(0, 1), mu, sigma, nu, tau), c(nu, 1))
    expect_equal(
        d01NVASIM(x, mu, sigma, nu, tau),
        continuous_probability * dNVASIM(x, mu, sigma)
    )
    expect_equal(
        p01NVASIM(x, mu, sigma, nu, tau),
        nu + continuous_probability * pNVASIM(x, mu, sigma)
    )
    expect_equal(
        q01NVASIM(mixture_probability, mu, sigma, nu, tau),
        qNVASIM(component_probability, mu, sigma),
        tolerance = 1e-12
    )
    expect_equal(q01NVASIM(c(0, nu), mu, sigma, nu, tau), c(0, 0))
    expect_equal(
        q01NVASIM(c(nu + continuous_probability, 1), mu, sigma, nu, tau),
        c(1, 1)
    )
    expect_equal(
        d01NVASIM(0, mu, sigma, nu, tau) +
            d01NVASIM(1, mu, sigma, nu, tau) +
            integrate(
                function(value) d01NVASIM(value, mu, sigma, nu, tau),
                lower = 0,
                upper = 1,
                subdivisions = 500L,
                rel.tol = 1e-10
            )$value,
        1,
        tolerance = 1e-7
    )
})
