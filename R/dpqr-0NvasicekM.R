#' @name ZANVASIM
#' @aliases ZANVASIM d0NVASIM p0NVASIM q0NVASIM r0NVASIM
#'
#' @title Zero-adjusted N-Vasicek distribution with mean parameterization
#'
#' @description
#' Defines a zero-adjusted normal-kernel Vasicek distribution for responses
#' in \eqn{[0,1)}. The parameter \eqn{\nu} is the probability of a structural
#' zero. Conditional on a positive response, the distribution is
#' \code{NVASIM} with mean \eqn{\mu} and shape parameter \eqn{\sigma}.
#'
#' @param x Vector of values in \eqn{[0,1]} at which the density or
#'   probability mass is evaluated. The distribution has support
#'   \eqn{[0,1)}, and the returned value is zero at \eqn{x=1}.
#' @param q Vector of values in \eqn{[0,1]} at which the cumulative
#'   distribution function is evaluated.
#' @param p Vector of probabilities.
#' @param n Number of observations. If \code{length(n) > 1}, its length is
#'   taken to be the number required.
#' @param mu Mean of the positive Vasicek component, in \eqn{(0,1)}.
#' @param sigma Shape parameter of the positive Vasicek component, in
#'   \eqn{(0,1)}.
#' @param nu Probability of a structural zero, in \eqn{(0,1)}.
#' @param log Logical; if \code{TRUE}, log probabilities or log densities are
#'   returned.
#' @param lower.tail Logical; if \code{TRUE}, probabilities are
#'   \eqn{P(Y\leq y)}; otherwise, they are \eqn{P(Y>y)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are supplied or
#'   returned on the log scale.
#' @param mu.link Link function for \eqn{\mu}.
#' @param sigma.link Link function for \eqn{\sigma}.
#' @param nu.link Link function for \eqn{\nu}.
#'
#' @details
#' Let \eqn{Y_+\sim\mathrm{NVASIM}(\mu,\sigma)} and let
#' \eqn{0<\nu<1}. The zero-adjusted distribution is defined by
#' \deqn{P(Y=0)=\nu}
#' and
#' \deqn{f_Y(y)=(1-\nu)f_{Y_+}(y\mid\mu,\sigma),\quad 0<y<1.}
#' Its cumulative distribution function is
#' \deqn{F_Y(y)=\nu+(1-\nu)F_{Y_+}(y\mid\mu,\sigma),\quad 0<y<1.}
#' Consequently,
#' \deqn{E(Y)=(1-\nu)\mu}
#' and
#' \deqn{\mathrm{Var}(Y)=(1-\nu)\mathrm{Var}(Y_+)+
#' \nu(1-\nu)\mu^2.}
#'
#' Thus, \eqn{\mu} is the mean conditional on \eqn{Y>0}; it is not the
#' marginal mean when \eqn{\nu>0}. The marginal mean is
#' \eqn{(1-\nu)\mu}.
#'
#' @return
#' \code{ZANVASIM()} returns a \code{gamlss.family} object. The functions
#' \code{d0NVASIM()}, \code{p0NVASIM()}, \code{q0NVASIM()}, and
#' \code{r0NVASIM()} return density or probability mass values, cumulative
#' probabilities, quantiles, and random observations, respectively.
#'
#' @references
#' Mazucheli, J., Alves, B., Korkmaz, M. C., and Leiva, V. (2022).
#' Vasicek quantile and mean regression models for bounded data: New
#' formulation, mathematical derivations, and numerical applications.
#' \emph{Mathematics}, \bold{10}, 1389.
#' \doi{10.3390/math10091389}
#'
#' Ospina, R. and Ferrari, S. L. P. (2010). Inflated beta distributions.
#' \emph{Statistical Papers}, \bold{51}, 111--126.
#'
#' Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive models
#' for location, scale and shape. \emph{Applied Statistics}, \bold{54}(3),
#' 507--554.
#'
#' @seealso \code{\link[vasicekreg]{NVASIM}},
#'   \code{\link[gamlss.dist]{BEZI}}
#'
#' @examples
#' set.seed(123)
#' y <- r0NVASIM(1000, mu = 0.60, sigma = 0.30, nu = 0.20)
#' mean(y == 0)
#' mean(y)
#' (1 - 0.20) * 0.60
#'
#' library(gamlss)
#' fit <- gamlss(
#'   y ~ 1,
#'   sigma.formula = ~ 1,
#'   nu.formula = ~ 1,
#'   family = ZANVASIM(),
#'   control = gamlss.control(trace = FALSE)
#' )
#' fitted(fit, what = "mu")[1]
#' fitted(fit, what = "sigma")[1]
#' fitted(fit, what = "nu")[1]
#'
NULL

.recycle_0nvasim <- function(..., .n = NULL) {
    values <- list(...)
    sizes <- vapply(values, length, integer(1))
    target <- if (is.null(.n)) max(sizes) else .n

    if (any(sizes != 1L & sizes != target)) {
        stop("Arguments must have length 1 or a common length.", call. = FALSE)
    }

    lapply(values, rep_len, length.out = target)
}
#' @rdname ZANVASIM
#' @export
d0NVASIM <- function(x, mu = 0.5, sigma = 0.5, nu = 0.1,
                     log = FALSE) {
    .check_unit_interval(x, "x", closed = TRUE)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    .check_unit_interval(nu, "nu")
    .check_scalar_logical(log, "log")

    if (length(x) == 0L) {
        return(numeric(0))
    }

    args <- .recycle_0nvasim(x, mu, sigma, nu)
    x <- args[[1L]]
    mu <- args[[2L]]
    sigma <- args[[3L]]
    nu <- args[[4L]]

    ans <- rep(-Inf, length(x))
    is_zero <- x == 0
    is_positive <- x > 0 & x < 1

    ans[is_zero] <- log(nu[is_zero])
    if (any(is_positive)) {
        ans[is_positive] <- log1p(-nu[is_positive]) +
            dNVASIM(
                x[is_positive],
                mu[is_positive],
                sigma[is_positive],
                log = TRUE
            )
    }

    if (log) ans else exp(ans)
}

#' @rdname ZANVASIM
#' @export
p0NVASIM <- function(q, mu = 0.5, sigma = 0.5, nu = 0.1,
                     lower.tail = TRUE, log.p = FALSE) {
    .check_unit_interval(q, "q", closed = TRUE)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    .check_unit_interval(nu, "nu")
    .check_scalar_logical(lower.tail, "lower.tail")
    .check_scalar_logical(log.p, "log.p")

    if (length(q) == 0L) {
        return(numeric(0))
    }

    args <- .recycle_0nvasim(q, mu, sigma, nu)
    q <- args[[1L]]
    mu <- args[[2L]]
    sigma <- args[[3L]]
    nu <- args[[4L]]

    cdf <- nu
    is_positive <- q > 0 & q < 1
    if (any(is_positive)) {
        cdf[is_positive] <- nu[is_positive] +
            (1 - nu[is_positive]) *
            pNVASIM(
                q[is_positive],
                mu[is_positive],
                sigma[is_positive]
            )
    }
    cdf[q >= 1] <- 1

    if (lower.tail) {
        if (log.p) log(cdf) else cdf
    } else {
        if (log.p) log1p(-cdf) else 1 - cdf
    }
}

#' @rdname ZANVASIM
#' @export
q0NVASIM <- function(p, mu = 0.5, sigma = 0.5, nu = 0.1,
                     lower.tail = TRUE, log.p = FALSE) {
    .check_probability(p, log.p)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    .check_unit_interval(nu, "nu")
    .check_scalar_logical(lower.tail, "lower.tail")
    .check_scalar_logical(log.p, "log.p")

    if (length(p) == 0L) {
        return(numeric(0))
    }

    probability <- if (log.p) exp(p) else p
    if (!lower.tail) {
        probability <- if (log.p) -expm1(p) else 1 - probability
    }

    args <- .recycle_0nvasim(probability, mu, sigma, nu)
    probability <- args[[1L]]
    mu <- args[[2L]]
    sigma <- args[[3L]]
    nu <- args[[4L]]

    ans <- numeric(length(probability))
    is_positive <- probability > nu
    if (any(is_positive)) {
        adjusted <- (probability[is_positive] - nu[is_positive]) /
            (1 - nu[is_positive])
        ans[is_positive] <- qNVASIM(
            adjusted,
            mu[is_positive],
            sigma[is_positive]
        )
    }
    ans
}

#' @rdname ZANVASIM
#' @export
r0NVASIM <- function(n, mu = 0.5, sigma = 0.5, nu = 0.1) {
    n <- .n_random(n)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    .check_unit_interval(nu, "nu")

    if (n == 0L) {
        return(numeric(0))
    }

    args <- .recycle_0nvasim(mu, sigma, nu, .n = n)
    q0NVASIM(
        stats::runif(n),
        mu = args[[1L]],
        sigma = args[[2L]],
        nu = args[[3L]]
    )
}

#' @rdname ZANVASIM
#' @export
ZANVASIM <- function(mu.link = "logit", sigma.link = "logit",
                     nu.link = "logit") {
    mstats <- checklink(
        "mu.link", "ZANVASIM", substitute(mu.link),
        c("logit", "probit", "cloglog", "cauchit", "log", "own")
    )
    dstats <- checklink(
        "sigma.link", "ZANVASIM", substitute(sigma.link),
        c("logit", "probit", "cloglog", "cauchit", "log", "own")
    )
    vstats <- checklink(
        "nu.link", "ZANVASIM", substitute(nu.link),
        c("logit", "probit", "cloglog", "cauchit", "log", "own")
    )

    structure(
        list(
            family = c("ZANVASIM", "Zero-adjusted N-VasicekM"),
            parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE),
            nopar = 3,
            # GAMLSS classification: point mass at zero plus a continuous
            # component on (0, 1). This does not specify random effects.
            type = "Mixed",
            mu.link = as.character(substitute(mu.link)),
            sigma.link = as.character(substitute(sigma.link)),
            nu.link = as.character(substitute(nu.link)),
            mu.linkfun = mstats$linkfun,
            sigma.linkfun = dstats$linkfun,
            nu.linkfun = vstats$linkfun,
            mu.linkinv = mstats$linkinv,
            sigma.linkinv = dstats$linkinv,
            nu.linkinv = vstats$linkinv,
            mu.dr = mstats$mu.eta,
            sigma.dr = dstats$mu.eta,
            nu.dr = vstats$mu.eta,
            dldm = function(y, mu, sigma) {
                ans <- numeric(length(y))
                ind <- y > 0
                if (any(ind)) {
                    yy <- y[ind]
                    mm <- rep_len(mu, length(y))[ind]
                    ss <- rep_len(sigma, length(y))[ind]
                    zmu <- qnorm(mm)
                    ans[ind] <- (qnorm(yy) * sqrt(1 - ss) - zmu) /
                        (ss * dnorm(zmu))
                }
                ans
            },
            d2ldm2 = function(y, mu, sigma) {
                ans <- numeric(length(y))
                ind <- y > 0
                if (any(ind)) {
                    yy <- y[ind]
                    mm <- rep_len(mu, length(y))[ind]
                    ss <- rep_len(sigma, length(y))[ind]
                    zmu <- qnorm(mm)
                    invphi2 <- 1 / dnorm(zmu)^2
                    ans[ind] <- invphi2 / ss *
                        (-1 + zmu * (qnorm(yy) * sqrt(1 - ss) - zmu))
                }
                ans
            },
            dldd = function(y, mu, sigma) {
                ans <- numeric(length(y))
                ind <- y > 0
                if (any(ind)) {
                    yy <- y[ind]
                    mm <- rep_len(mu, length(y))[ind]
                    ss <- rep_len(sigma, length(y))[ind]
                    zy <- qnorm(yy)
                    a <- sqrt(1 - ss)
                    b <- zy * a - qnorm(mm)
                    ans[ind] <- -0.5 / (1 - ss) - 0.5 / ss +
                        0.5 * b * zy / (ss * a) +
                        0.5 * b^2 / ss^2
                }
                ans
            },
            d2ldd2 = function(y, mu, sigma) {
                ans <- numeric(length(y))
                ind <- y > 0
                if (any(ind)) {
                    yy <- y[ind]
                    mm <- rep_len(mu, length(y))[ind]
                    ss <- rep_len(sigma, length(y))[ind]
                    zy <- qnorm(yy)
                    a <- sqrt(1 - ss)
                    b <- zy * a - qnorm(mm)
                    ans[ind] <- -0.5 / (1 - ss)^2 + 0.5 / ss^2 -
                        0.25 * zy^2 / ((1 - ss) * ss) -
                        b * zy / (ss^2 * a) +
                        0.25 * b * zy / (ss * a * (1 - ss)) -
                        b^2 / ss^3
                }
                ans
            },
            d2ldmdd = function(y, mu, sigma) {
                ans <- numeric(length(y))
                ind <- y > 0
                if (any(ind)) {
                    yy <- y[ind]
                    mm <- rep_len(mu, length(y))[ind]
                    ss <- rep_len(sigma, length(y))[ind]
                    zy <- qnorm(yy)
                    a <- sqrt(1 - ss)
                    zmu <- qnorm(mm)
                    ans[ind] <- -0.5 * zy / (a * ss * dnorm(zmu)) -
                        (zy * a - zmu) / (ss^2 * dnorm(zmu))
                }
                ans
            },
            dldv = function(y, nu) {
                ifelse(y == 0, 1 / nu, -1 / (1 - nu))
            },
            d2ldv2 = function(nu) {
                -1 / (nu * (1 - nu))
            },
            d2ldmdv = function(y) numeric(length(y)),
            d2ldddv = function(y) numeric(length(y)),
            G.dev.incr = function(y, mu, sigma, nu, ...) {
                -2 * d0NVASIM(y, mu, sigma, nu, log = TRUE)
            },
            rqres = expression({
                uval <- ifelse(
                    y == 0,
                    stats::runif(length(y), min = 0, max = nu),
                    p0NVASIM(y, mu, sigma, nu)
                )
                stats::qnorm(uval)
            }),
            mu.initial = expression({
                positive_y <- y[y > 0]
                positive_mean <- if (length(positive_y)) {
                    mean(positive_y)
                } else {
                    0.5
                }
                mu <- (y + positive_mean) / 2
            }),
            sigma.initial = expression({
                sigma <- rep(0.5, length(y))
            }),
            nu.initial = expression({
                zero_fraction <- mean(y == 0)
                zero_fraction <- min(max(zero_fraction, 0.01), 0.99)
                nu <- rep(zero_fraction, length(y))
            }),
            mu.valid = function(mu) all(mu > 0 & mu < 1),
            sigma.valid = function(sigma) all(sigma > 0 & sigma < 1),
            nu.valid = function(nu) all(nu > 0 & nu < 1),
            y.valid = function(y) all(y >= 0 & y < 1),
            mean = function(mu, sigma, nu) {
                args <- .recycle_0nvasim(mu, sigma, nu)
                (1 - args[[3L]]) * args[[1L]]
            },
            variance = function(mu, sigma, nu) {
                args <- .recycle_0nvasim(mu, sigma, nu)
                mu <- args[[1L]]
                sigma <- args[[2L]]
                nu <- args[[3L]]
                positive_variance <- NVASIM()$variance(mu, sigma)
                (1 - nu) * positive_variance +
                    nu * (1 - nu) * mu^2
            }
        ),
        class = c("gamlss.family", "family")
    )
}
