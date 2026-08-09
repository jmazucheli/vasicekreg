#' @name ZOANVASIM
#' @aliases ZOANVASIM d01NVASIM p01NVASIM q01NVASIM r01NVASIM
#'
#' @title Zero-and-one-adjusted N-Vasicek distribution
#'
#' @description
#' Defines a normal-kernel Vasicek distribution augmented by point masses at
#' zero and one. Conditional on an observation in \eqn{(0,1)}, the continuous
#' component is \code{NVASIM} with mean \eqn{\mu} and shape parameter
#' \eqn{\sigma}.
#'
#' @param x Vector of values in \eqn{[0,1]} at which the density or
#'   probability mass is evaluated. The distribution has support
#'   \eqn{[0,1]}.
#' @param q Vector of values in \eqn{[0,1]} at which the cumulative
#'   distribution function is evaluated.
#' @param p Vector of probabilities.
#' @param n Number of observations. If \code{length(n) > 1}, its length is
#'   taken to be the number required.
#' @param mu Mean of the continuous Vasicek component, in \eqn{(0,1)}.
#' @param sigma Shape parameter of the continuous Vasicek component, in
#'   \eqn{(0,1)}.
#' @param nu Probability at zero, \eqn{\nu=P(Y=0)}, in \eqn{(0,1)}.
#' @param tau Conditional probability at one among nonzero observations,
#'   \eqn{\tau=P(Y=1\mid Y>0)}, in \eqn{(0,1)}. This parameter is unrelated
#'   to the fixed quantile level used by \code{NVASIQ()} and \code{LVASIQ()}.
#' @param log Logical; if \code{TRUE}, log probabilities or log densities are
#'   returned.
#' @param lower.tail Logical; if \code{TRUE}, probabilities are
#'   \eqn{P(Y\leq y)}; otherwise, they are \eqn{P(Y>y)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are supplied or
#'   returned on the log scale.
#' @param mu.link Link function for \eqn{\mu}.
#' @param sigma.link Link function for \eqn{\sigma}.
#' @param nu.link Link function for \eqn{\nu}.
#' @param tau.link Link function for \eqn{\tau}.
#'
#' @details
#' Let \eqn{Y_c\sim\mathrm{NVASIM}(\mu,\sigma)}. Write \eqn{p_0}, \eqn{p_1},
#' and \eqn{p_c} for the probabilities of zero, one, and the continuous
#' component. The sequential BEOI-type parameterization is
#' \deqn{\nu=P(Y=0),\qquad \tau=P(Y=1\mid Y>0).}
#' Hence,
#' \deqn{p_0=\nu,\qquad
#' p_1=(1-\nu)\tau,\qquad
#' p_c=(1-\nu)(1-\tau).}
#' The distribution is
#' \deqn{P(Y=0)=p_0,\qquad P(Y=1)=p_1}
#' and
#' \deqn{f_Y(y)=p_c f_{Y_c}(y\mid\mu,\sigma),\quad 0<y<1.}
#' Its marginal mean and variance are
#' \deqn{E(Y)=(1-\nu)\left[\tau+(1-\tau)\mu\right]}
#' and
#' \deqn{\mathrm{Var}(Y)=
#' (1-\nu)\left[(1-\tau)\left\{\mathrm{Var}(Y_c)+\mu^2\right\}
#' +\tau\right]
#' -\left\{(1-\nu)\left[\tau+(1-\tau)\mu\right]\right\}^2.}
#' Logit links for \eqn{\nu} and \eqn{\tau} guarantee valid probabilities.
#' If \eqn{\nu=0}, the model reduces to the BEOI-type one-adjusted model; if
#' \eqn{\tau=0}, it reduces to the zero-adjusted model.
#'
#' @return
#' \code{ZOANVASIM()} returns a four-parameter \code{gamlss.family} object.
#' The functions \code{d01NVASIM()}, \code{p01NVASIM()},
#' \code{q01NVASIM()}, and \code{r01NVASIM()} return probability mass or
#' density values, cumulative probabilities, quantiles, and random
#' observations, respectively.
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
#'   \code{\link[vasicekreg]{OANVASIM}},
#'   \code{\link[gamlss.dist]{BEOI}}
#'
#' @examples
#' set.seed(123)
#' y <- r01NVASIM(
#'   1000, mu = 0.60, sigma = 0.30, nu = 0.20, tau = 0.25
#' )
#' c(zero = mean(y == 0), one = mean(y == 1))
#' mean(y)
#' (1 - 0.20) * (0.25 + (1 - 0.25) * 0.60)
#'
#' \dontrun{
#' library(gamlss)
#' fit <- gamlss(
#'   y ~ 1,
#'   sigma.formula = ~ 1,
#'   nu.formula = ~ 1,
#'   tau.formula = ~ 1,
#'   family = ZOANVASIM(),
#'   control = gamlss.control(trace = FALSE)
#' )
#' }
NULL

.recycle_01nvasim <- function(..., .n = NULL) {
    values <- list(...)
    sizes <- vapply(values, length, integer(1))
    target <- if (is.null(.n)) max(sizes) else .n

    if (any(sizes != 1L & sizes != target)) {
        stop("Arguments must have length 1 or a common length.", call. = FALSE)
    }

    lapply(values, rep_len, length.out = target)
}
#' @rdname ZOANVASIM
#' @export
d01NVASIM <- function(x, mu = 0.5, sigma = 0.5,
                      nu = 0.1, tau = 0.1, log = FALSE) {
    .check_unit_interval(x, "x", closed = TRUE)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    .check_unit_interval(nu, "nu")
    .check_unit_interval(tau, "tau")
    .check_scalar_logical(log, "log")

    if (length(x) == 0L) return(numeric(0))

    args <- .recycle_01nvasim(x, mu, sigma, nu, tau)
    x <- args[[1L]]
    mu <- args[[2L]]
    sigma <- args[[3L]]
    nu <- args[[4L]]
    tau <- args[[5L]]
    ans <- rep(-Inf, length(x))
    at_zero <- x == 0
    at_one <- x == 1
    interior <- x > 0 & x < 1
    ans[at_zero] <- log(nu[at_zero])
    ans[at_one] <- log1p(-nu[at_one]) + log(tau[at_one])
    if (any(interior)) {
        ans[interior] <- log1p(-nu[interior]) +
            log1p(-tau[interior]) +
            dNVASIM(
                x[interior], mu[interior], sigma[interior], log = TRUE
            )
    }

    if (log) ans else exp(ans)
}

#' @rdname ZOANVASIM
#' @export
p01NVASIM <- function(q, mu = 0.5, sigma = 0.5,
                      nu = 0.1, tau = 0.1,
                      lower.tail = TRUE, log.p = FALSE) {
    .check_unit_interval(q, "q", closed = TRUE)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    .check_unit_interval(nu, "nu")
    .check_unit_interval(tau, "tau")
    .check_scalar_logical(lower.tail, "lower.tail")
    .check_scalar_logical(log.p, "log.p")

    if (length(q) == 0L) return(numeric(0))

    args <- .recycle_01nvasim(q, mu, sigma, nu, tau)
    q <- args[[1L]]
    mu <- args[[2L]]
    sigma <- args[[3L]]
    nu <- args[[4L]]
    tau <- args[[5L]]
    cdf <- nu
    interior <- q > 0 & q < 1
    if (any(interior)) {
        cdf[interior] <- nu[interior] +
            (1 - nu[interior]) * (1 - tau[interior]) *
            pNVASIM(q[interior], mu[interior], sigma[interior])
    }
    cdf[q >= 1] <- 1

    if (lower.tail) {
        if (log.p) log(cdf) else cdf
    } else {
        if (log.p) log1p(-cdf) else 1 - cdf
    }
}

#' @rdname ZOANVASIM
#' @export
q01NVASIM <- function(p, mu = 0.5, sigma = 0.5,
                      nu = 0.1, tau = 0.1,
                      lower.tail = TRUE, log.p = FALSE) {
    .check_probability(p, log.p)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    .check_unit_interval(nu, "nu")
    .check_unit_interval(tau, "tau")
    .check_scalar_logical(lower.tail, "lower.tail")
    .check_scalar_logical(log.p, "log.p")

    if (length(p) == 0L) return(numeric(0))

    probability <- if (log.p) exp(p) else p
    if (!lower.tail) {
        probability <- if (log.p) -expm1(p) else 1 - probability
    }

    args <- .recycle_01nvasim(probability, mu, sigma, nu, tau)
    probability <- args[[1L]]
    mu <- args[[2L]]
    sigma <- args[[3L]]
    nu <- args[[4L]]
    tau <- args[[5L]]
    p_zero <- nu
    p_continuous <- (1 - nu) * (1 - tau)
    p_below_one <- p_zero + p_continuous

    ans <- numeric(length(probability))
    at_one <- probability >= p_below_one
    ans[at_one] <- 1
    interior <- probability > p_zero & probability < p_below_one
    if (any(interior)) {
        adjusted <- (probability[interior] - nu[interior]) /
            p_continuous[interior]
        ans[interior] <- qNVASIM(
            adjusted, mu[interior], sigma[interior]
        )
    }
    ans
}

#' @rdname ZOANVASIM
#' @export
r01NVASIM <- function(n, mu = 0.5, sigma = 0.5,
                      nu = 0.1, tau = 0.1) {
    n <- .n_random(n)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    .check_unit_interval(nu, "nu")
    .check_unit_interval(tau, "tau")

    if (n == 0L) return(numeric(0))

    args <- .recycle_01nvasim(mu, sigma, nu, tau, .n = n)
    q01NVASIM(
        stats::runif(n), args[[1L]], args[[2L]], args[[3L]], args[[4L]]
    )
}

#' @rdname ZOANVASIM
#' @export
ZOANVASIM <- function(mu.link = "logit", sigma.link = "logit",
                      nu.link = "logit", tau.link = "logit") {
    mstats <- checklink(
        "mu.link", "ZOANVASIM", substitute(mu.link),
        c("logit", "probit", "cloglog", "cauchit", "log", "own")
    )
    dstats <- checklink(
        "sigma.link", "ZOANVASIM", substitute(sigma.link),
        c("logit", "probit", "cloglog", "cauchit", "log", "own")
    )
    vstats <- checklink(
        "nu.link", "ZOANVASIM", substitute(nu.link),
        c("logit", "probit", "cloglog", "cauchit", "log", "own")
    )
    tstats <- checklink(
        "tau.link", "ZOANVASIM", substitute(tau.link),
        c("logit", "probit", "cloglog", "cauchit", "log", "own")
    )

    structure(
        list(
            family = c("ZOANVASIM", "Zero-and-one-adjusted N-VasicekM"),
            parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE, tau = TRUE),
            nopar = 4,
            # Boundary masses plus a continuous component; no random effects.
            type = "Mixed",
            mu.link = as.character(substitute(mu.link)),
            sigma.link = as.character(substitute(sigma.link)),
            nu.link = as.character(substitute(nu.link)),
            tau.link = as.character(substitute(tau.link)),
            mu.linkfun = mstats$linkfun,
            sigma.linkfun = dstats$linkfun,
            nu.linkfun = vstats$linkfun,
            tau.linkfun = tstats$linkfun,
            mu.linkinv = mstats$linkinv,
            sigma.linkinv = dstats$linkinv,
            nu.linkinv = vstats$linkinv,
            tau.linkinv = tstats$linkinv,
            mu.dr = mstats$mu.eta,
            sigma.dr = dstats$mu.eta,
            nu.dr = vstats$mu.eta,
            tau.dr = tstats$mu.eta,
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
                ifelse(y == 0, 1 / nu, -1 / (1 - nu))
            },
            d2ldv2 = function(nu) {
                -1 / (nu * (1 - nu))
            },
            dldt = function(y, nu, tau) {
                ifelse(y == 0, 0, ifelse(y == 1, 1 / tau, -1 / (1 - tau)))
            },
            d2ldt2 = function(nu, tau) {
                -(1 - nu) / (tau * (1 - tau))
            },
            d2ldmdv = function(y) numeric(length(y)),
            d2ldmdt = function(y) numeric(length(y)),
            d2ldddv = function(y) numeric(length(y)),
            d2ldddt = function(y) numeric(length(y)),
            d2ldvdt = function(nu, tau) {
                rep(0, max(length(nu), length(tau)))
            },
            G.dev.incr = function(y, mu, sigma, nu, tau, ...) {
                -2 * d01NVASIM(y, mu, sigma, nu, tau, log = TRUE)
            },
            rqres = expression({
                pone.lower <- nu + (1 - nu) * (1 - tau)
                uval <- ifelse(
                    y == 0,
                    stats::runif(length(y), min = 0, max = nu),
                    ifelse(
                        y == 1,
                        stats::runif(length(y), min = pone.lower, max = 1),
                        p01NVASIM(y, mu, sigma, nu, tau)
                    )
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
                pzero <- min(max(mean(y == 0), 0.01), 0.99)
                nu <- rep(pzero, length(y))
            }),
            tau.initial = expression({
                nonzero <- sum(y > 0)
                pone.given.nonzero <- if (nonzero > 0) {
                    sum(y == 1) / nonzero
                } else {
                    0.5
                }
                pone.given.nonzero <- min(
                    max(pone.given.nonzero, 0.01), 0.99
                )
                tau <- rep(pone.given.nonzero, length(y))
            }),
            mu.valid = function(mu) all(mu > 0 & mu < 1),
            sigma.valid = function(sigma) all(sigma > 0 & sigma < 1),
            nu.valid = function(nu) all(nu > 0 & nu < 1),
            tau.valid = function(tau) all(tau > 0 & tau < 1),
            y.valid = function(y) all(y >= 0 & y <= 1),
            mean = function(mu, sigma, nu, tau) {
                args <- .recycle_01nvasim(mu, sigma, nu, tau)
                mu <- args[[1L]]
                nu <- args[[3L]]
                tau <- args[[4L]]
                (1 - nu) * (tau + (1 - tau) * mu)
            },
            variance = function(mu, sigma, nu, tau) {
                args <- .recycle_01nvasim(mu, sigma, nu, tau)
                mu <- args[[1L]]
                sigma <- args[[2L]]
                nu <- args[[3L]]
                tau <- args[[4L]]
                positive_variance <- NVASIM()$variance(mu, sigma)
                second_moment <- (1 - nu) * (
                    (1 - tau) * (positive_variance + mu^2) + tau
                )
                marginal_mean <- (1 - nu) * (
                    tau + (1 - tau) * mu
                )
                second_moment - marginal_mean^2
            }
        ),
        class = c("gamlss.family", "family")
    )
}
