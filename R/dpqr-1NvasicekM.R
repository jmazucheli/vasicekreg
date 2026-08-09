#' @name OANVASIM
#' @aliases OANVASIM d1NVASIM p1NVASIM q1NVASIM r1NVASIM dOANVASIM pOANVASIM qOANVASIM rOANVASIM
#'
#' @title One-adjusted N-Vasicek distribution with mean parameterization
#'
#' @description
#' Defines a one-adjusted normal-kernel Vasicek distribution for responses
#' in \eqn{(0,1]}. The parameter \eqn{\nu} is the probability at one.
#' Conditional on an observation in \eqn{(0,1)}, the distribution is
#' \code{NVASIM} with mean \eqn{\mu} and shape parameter \eqn{\sigma}.
#'
#' @param x Vector of values in \eqn{[0,1]} at which the density or
#'   probability mass is evaluated. The distribution has support
#'   \eqn{(0,1]}, and the returned value is zero at \eqn{x=0}.
#' @param q Vector of values in \eqn{[0,1]} at which the cumulative
#'   distribution function is evaluated.
#' @param p Vector of probabilities.
#' @param n Number of observations. If \code{length(n) > 1}, its length is
#'   taken to be the number required.
#' @param mu Mean of the continuous Vasicek component, in \eqn{(0,1)}.
#' @param sigma Shape parameter of the continuous Vasicek component, in
#'   \eqn{(0,1)}.
#' @param nu Probability at one, in \eqn{(0,1)}.
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
#' Let \eqn{Y_c\sim\mathrm{NVASIM}(\mu,\sigma)} and let
#' \eqn{0<\nu<1}. The BEOI-type one-adjusted distribution is defined by
#' \deqn{P(Y=1)=\nu}
#' and
#' \deqn{f_Y(y)=(1-\nu)f_{Y_c}(y\mid\mu,\sigma),\quad 0<y<1.}
#' Consequently,
#' \deqn{E(Y)=\nu+(1-\nu)\mu}
#' and
#' \deqn{\mathrm{Var}(Y)=(1-\nu)\mathrm{Var}(Y_c)+
#' \nu(1-\nu)(1-\mu)^2.}
#' Thus, \eqn{\mu=E(Y\mid 0<Y<1)} is the mean of the continuous component,
#' whereas \eqn{\nu+(1-\nu)\mu} is the marginal mean.
#'
#' @return
#' \code{OANVASIM()} returns a \code{gamlss.family} object. The functions
#' \code{d1NVASIM()}, \code{p1NVASIM()}, \code{q1NVASIM()}, and
#' \code{r1NVASIM()} return probability mass or density values, cumulative
#' probabilities, quantiles, and random observations, respectively.
#' \code{dOANVASIM()}, \code{pOANVASIM()}, \code{qOANVASIM()}, and
#' \code{rOANVASIM()} are equivalent names following the GAMLSS family-name
#' convention.
#'
#' @references
#' Ospina, R. and Ferrari, S. L. P. (2010). Inflated beta distributions.
#' \emph{Statistical Papers}, \bold{51}, 111--126.
#'
#' Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive models
#' for location, scale and shape. \emph{Applied Statistics}, \bold{54}(3),
#' 507--554.
#'
#' @seealso \code{\link[vasicekreg]{NVASIM}},
#'   \code{\link[vasicekreg]{ZANVASIM}},
#'   \code{\link[gamlss.dist]{BEOI}}
#'
#' @examples
#' set.seed(123)
#' y <- r1NVASIM(1000, mu = 0.60, sigma = 0.30, nu = 0.20)
#' mean(y == 1)
#' mean(y)
#' 0.20 + (1 - 0.20) * 0.60
#'
#' \dontrun{
#' library(gamlss)
#' fit <- gamlss(
#'   y ~ 1,
#'   sigma.formula = ~ 1,
#'   nu.formula = ~ 1,
#'   family = OANVASIM(),
#'   control = gamlss.control(trace = FALSE)
#' )
#' }
NULL

.recycle_1nvasim <- function(..., .n = NULL) {
    values <- list(...)
    sizes <- vapply(values, length, integer(1))
    target <- if (is.null(.n)) max(sizes) else .n

    if (any(sizes != 1L & sizes != target)) {
        stop("Arguments must have length 1 or a common length.", call. = FALSE)
    }

    lapply(values, rep_len, length.out = target)
}

#' @rdname OANVASIM
#' @export
d1NVASIM <- function(x, mu = 0.5, sigma = 0.5, nu = 0.1,
                     log = FALSE) {
    .check_unit_interval(x, "x", closed = TRUE)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    .check_unit_interval(nu, "nu")
    .check_scalar_logical(log, "log")

    if (length(x) == 0L) return(numeric(0))

    args <- .recycle_1nvasim(x, mu, sigma, nu)
    x <- args[[1L]]
    mu <- args[[2L]]
    sigma <- args[[3L]]
    nu <- args[[4L]]

    ans <- rep(-Inf, length(x))
    at_one <- x == 1
    interior <- x > 0 & x < 1
    ans[at_one] <- log(nu[at_one])
    if (any(interior)) {
        ans[interior] <- log1p(-nu[interior]) +
            dNVASIM(x[interior], mu[interior], sigma[interior], log = TRUE)
    }

    if (log) ans else exp(ans)
}

#' @rdname OANVASIM
#' @export
p1NVASIM <- function(q, mu = 0.5, sigma = 0.5, nu = 0.1,
                     lower.tail = TRUE, log.p = FALSE) {
    .check_unit_interval(q, "q", closed = TRUE)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    .check_unit_interval(nu, "nu")
    .check_scalar_logical(lower.tail, "lower.tail")
    .check_scalar_logical(log.p, "log.p")

    if (length(q) == 0L) return(numeric(0))

    args <- .recycle_1nvasim(q, mu, sigma, nu)
    q <- args[[1L]]
    mu <- args[[2L]]
    sigma <- args[[3L]]
    nu <- args[[4L]]

    cdf <- numeric(length(q))
    interior <- q > 0 & q < 1
    if (any(interior)) {
        cdf[interior] <- (1 - nu[interior]) *
            pNVASIM(q[interior], mu[interior], sigma[interior])
    }
    cdf[q >= 1] <- 1

    if (lower.tail) {
        if (log.p) log(cdf) else cdf
    } else {
        if (log.p) log1p(-cdf) else 1 - cdf
    }
}

#' @rdname OANVASIM
#' @export
q1NVASIM <- function(p, mu = 0.5, sigma = 0.5, nu = 0.1,
                     lower.tail = TRUE, log.p = FALSE) {
    .check_probability(p, log.p)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    .check_unit_interval(nu, "nu")
    .check_scalar_logical(lower.tail, "lower.tail")
    .check_scalar_logical(log.p, "log.p")

    if (length(p) == 0L) return(numeric(0))

    probability <- if (log.p) exp(p) else p
    if (!lower.tail) {
        probability <- if (log.p) -expm1(p) else 1 - probability
    }

    args <- .recycle_1nvasim(probability, mu, sigma, nu)
    probability <- args[[1L]]
    mu <- args[[2L]]
    sigma <- args[[3L]]
    nu <- args[[4L]]

    ans <- rep(1, length(probability))
    interior <- probability < 1 - nu
    if (any(interior)) {
        ans[interior] <- qNVASIM(
            probability[interior] / (1 - nu[interior]),
            mu[interior], sigma[interior]
        )
    }
    ans
}

#' @rdname OANVASIM
#' @export
r1NVASIM <- function(n, mu = 0.5, sigma = 0.5, nu = 0.1) {
    n <- .n_random(n)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    .check_unit_interval(nu, "nu")

    if (n == 0L) return(numeric(0))

    args <- .recycle_1nvasim(mu, sigma, nu, .n = n)
    q1NVASIM(stats::runif(n), args[[1L]], args[[2L]], args[[3L]])
}

# GAMLSS derives distribution-function names from family[[1]]. These
# wrappers preserve the established numbered API while satisfying that
# convention.
#' @rdname OANVASIM
#' @export
dOANVASIM <- function(x, mu = 0.5, sigma = 0.5, nu = 0.1,
                      log = FALSE) {
    d1NVASIM(x, mu, sigma, nu, log)
}

#' @rdname OANVASIM
#' @export
pOANVASIM <- function(q, mu = 0.5, sigma = 0.5, nu = 0.1,
                      lower.tail = TRUE, log.p = FALSE) {
    p1NVASIM(q, mu, sigma, nu, lower.tail, log.p)
}

#' @rdname OANVASIM
#' @export
qOANVASIM <- function(p, mu = 0.5, sigma = 0.5, nu = 0.1,
                      lower.tail = TRUE, log.p = FALSE) {
    q1NVASIM(p, mu, sigma, nu, lower.tail, log.p)
}

#' @rdname OANVASIM
#' @export
rOANVASIM <- function(n, mu = 0.5, sigma = 0.5, nu = 0.1) {
    r1NVASIM(n, mu, sigma, nu)
}

#' @rdname OANVASIM
#' @export
OANVASIM <- function(mu.link = "logit", sigma.link = "logit",
                     nu.link = "logit") {
    mstats <- checklink(
        "mu.link", "OANVASIM", substitute(mu.link),
        c("logit", "probit", "cloglog", "cauchit", "log", "own")
    )
    dstats <- checklink(
        "sigma.link", "OANVASIM", substitute(sigma.link),
        c("logit", "probit", "cloglog", "cauchit", "log", "own")
    )
    vstats <- checklink(
        "nu.link", "OANVASIM", substitute(nu.link),
        c("logit", "probit", "cloglog", "cauchit", "log", "own")
    )

    structure(
        list(
            family = c("OANVASIM", "One-adjusted N-VasicekM"),
            parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE),
            nopar = 3,
            # Point mass at one plus a continuous component; no random effects.
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
                ind <- y > 0 & y < 1
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
                ind <- y > 0 & y < 1
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
                ind <- y > 0 & y < 1
                if (any(ind)) {
                    yy <- y[ind]
                    mm <- rep_len(mu, length(y))[ind]
                    ss <- rep_len(sigma, length(y))[ind]
                    zy <- qnorm(yy)
                    a <- sqrt(1 - ss)
                    b <- zy * a - qnorm(mm)
                    ans[ind] <- -0.5 / (1 - ss) - 0.5 / ss +
                        0.5 * b * zy / (ss * a) + 0.5 * b^2 / ss^2
                }
                ans
            },
            d2ldd2 = function(y, mu, sigma) {
                ans <- numeric(length(y))
                ind <- y > 0 & y < 1
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
                ind <- y > 0 & y < 1
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
                ifelse(y == 1, 1 / nu, -1 / (1 - nu))
            },
            d2ldv2 = function(nu) {
                -1 / (nu * (1 - nu))
            },
            d2ldmdv = function(y) numeric(length(y)),
            d2ldddv = function(y) numeric(length(y)),
            G.dev.incr = function(y, mu, sigma, nu, ...) {
                -2 * d1NVASIM(y, mu, sigma, nu, log = TRUE)
            },
            rqres = expression({
                uval <- ifelse(
                    y == 1,
                    stats::runif(length(y), min = 1 - nu, max = 1),
                    p1NVASIM(y, mu, sigma, nu)
                )
                stats::qnorm(uval)
            }),
            mu.initial = expression({
                interior_y <- y[y > 0 & y < 1]
                interior_mean <- if (length(interior_y)) mean(interior_y) else 0.5
                mu <- (y + interior_mean) / 2
            }),
            sigma.initial = expression({
                sigma <- rep(0.5, length(y))
            }),
            nu.initial = expression({
                one_fraction <- min(max(mean(y == 1), 0.01), 0.99)
                nu <- rep(one_fraction, length(y))
            }),
            mu.valid = function(mu) all(mu > 0 & mu < 1),
            sigma.valid = function(sigma) all(sigma > 0 & sigma < 1),
            nu.valid = function(nu) all(nu > 0 & nu < 1),
            y.valid = function(y) all(y > 0 & y <= 1),
            mean = function(mu, sigma, nu) {
                args <- .recycle_1nvasim(mu, sigma, nu)
                args[[3L]] + (1 - args[[3L]]) * args[[1L]]
            },
            variance = function(mu, sigma, nu) {
                args <- .recycle_1nvasim(mu, sigma, nu)
                mu <- args[[1L]]
                sigma <- args[[2L]]
                nu <- args[[3L]]
                positive_variance <- NVASIM()$variance(mu, sigma)
                (1 - nu) * positive_variance +
                    nu * (1 - nu) * (1 - mu)^2
            }
        ),
        class = c("gamlss.family", "family")
    )
}
