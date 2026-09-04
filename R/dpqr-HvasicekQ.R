#' Hyperbolic-secant-kernel Vasicek-type quantile distribution
#'
#' @description
#' \code{HVASIQ()} defines the hyperbolic-secant-kernel Vasicek-type distribution
#' as a \code{gamlss.family} object for conditional quantile regression.
#' The functions \code{dHVASIQ()}, \code{pHVASIQ()}, \code{qHVASIQ()}, and
#' \code{rHVASIQ()} provide the density, distribution function, quantile
#' function, and random generation.
#'
#' @details
#' Let
#' \deqn{H(w)=\frac{2}{\pi}\arctan\{\exp(w)\},\qquad
#' h(w)=\frac{1}{\pi\cosh(w)},}
#' and
#' \deqn{q(u)=H^{-1}(u)=
#' \log\left\{\tan\left(\frac{\pi u}{2}\right)\right\}.}
#' For fixed \eqn{\tau\in(0,1)}, define
#' \deqn{a=\sqrt{\frac{1-\sigma}{\sigma}},\qquad
#' z=a\{q(x)-q(\mu)\}+q(\tau).}
#' The cumulative distribution function and density are
#' \deqn{F(x\mid\mu,\sigma,\tau)=H(z)}
#' and
#' \deqn{f(x\mid\mu,\sigma,\tau)=a\frac{h(z)}{h\{q(x)\}}.}
#' The quantile function is
#' \deqn{Q(p\mid\mu,\sigma,\tau)=
#' H\left[q(\mu)+\sqrt{\frac{\sigma}{1-\sigma}}
#' \{q(p)-q(\tau)\}\right].}
#' Consequently, \eqn{Q(\tau)=\mu}; \eqn{\mu} is exactly the conditional
#' \eqn{\tau}-th quantile and is not the conditional mean in general.
#'
#' The GAMLSS family uses analytical derivatives. For one observation, let
#' \deqn{d=q(y)-q(\mu),\quad z=ad+q(\tau),\quad
#' T=\tanh(z),\quad S=\operatorname{sech}^2(z),}
#' \deqn{b=\frac{1}{2\sigma(1-\sigma)},\qquad
#' g=\frac{\pi}{\sin(\pi\mu)}.}
#' If \eqn{\ell} is the individual log-likelihood contribution, then
#' \deqn{\frac{\partial\ell}{\partial\mu}=agT,}
#' \deqn{\frac{\partial\ell}{\partial\sigma}=b(adT-1),}
#' \deqn{\frac{\partial^2\ell}{\partial\mu^2}
#' =g^2\{-a\cos(\pi\mu)T-a^2S\},}
#' \deqn{\frac{\partial^2\ell}{\partial\mu\,\partial\sigma}
#' =-abg(T+adS),}
#' and
#' \deqn{\frac{\partial^2\ell}{\partial\sigma^2}
#' =b^2\{2(1-2\sigma)-(3-4\sigma)adT-a^2d^2S\}.}
#' The conditional mean and variance do not have elementary closed forms and
#' are evaluated by numerical integration of the quantile function.
#'
#' The fixed level is supplied through the \code{quantile} argument. It is
#' stored in the family definition and embedded as a numeric literal in the
#' family components used by GAMLSS; no global variable is required.
#'
#' @param x Vector of quantiles in \eqn{(0,1)}.
#' @param q Vector of values in \eqn{[0,1]} at which the cumulative
#'   distribution function is evaluated.
#' @param p Vector of probabilities in \eqn{[0,1]} on the probability scale.
#' @param n Number of observations.
#' @param mu Vector of conditional \eqn{\tau}-th quantiles,
#'   \eqn{0<\mu<1}.
#' @param sigma Vector of shape parameter values, \eqn{0<\sigma<1}.
#' @param quantile Fixed quantile level \eqn{\tau\in(0,1)} represented by
#'   \eqn{\mu}, used in the distribution functions and in \code{HVASIQ()}.
#' @param mu.link Link function for \eqn{\mu}.
#' @param sigma.link Link function for \eqn{\sigma}.
#' @param lower.tail Logical; if \code{TRUE}, probabilities are
#'   \eqn{P(X\leq x)}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @param log.p Logical; if \code{TRUE}, probabilities are supplied or
#'   returned on the logarithmic scale.
#'
#' @return
#' \code{dHVASIQ()} returns the density, \code{pHVASIQ()} the distribution
#' function, \code{qHVASIQ()} the quantile function, and \code{rHVASIQ()}
#' random deviates. \code{HVASIQ()} returns a \code{gamlss.family} object.
#'
#' @references
#' Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals.
#' \emph{Journal of Computational and Graphical Statistics}, \bold{5}(3),
#' 236--244.
#'
#' Fischer, M. J., Hui, A. and Hösle, S. (2017). wHS-type distributions
#' with application to finance. \emph{Journal of Statistics and Management
#' Systems}, \bold{20}(1), 67--89.
#' \doi{10.1080/09720510.2016.1190575}
#'
#' Mazucheli, J., Alves, B., Korkmaz, M. C. and Leiva, V. (2022).
#' Vasicek quantile and mean regression models for bounded data:
#' New formulation, mathematical derivations, and numerical applications.
#' \emph{Mathematics}, \bold{10}, 1389.
#' \doi{10.3390/math10091389}
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#'
#' @examples
#' set.seed(123)
#' y <- rHVASIQ(500, mu = 0.60, sigma = 0.30, quantile = 0.25)
#'
#' fit <- gamlss::gamlss(
#'     y ~ 1,
#'     sigma.formula = ~ 1,
#'     family = HVASIQ(quantile = 0.25),
#'     control = gamlss::gamlss.control(trace = FALSE)
#' )
#' fitted(fit, what = "mu")[1]
#'
#' @name HVASIQ
#' @aliases HVASIQ dHVASIQ pHVASIQ qHVASIQ rHVASIQ
#' @importFrom gamlss gamlss
#' @importFrom gamlss.dist checklink
#' @export
dHVASIQ <- function(x, mu, sigma, quantile = 0.5, log = FALSE) {
    .check_unit_interval(x, "x")
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    quantile <- .check_fixed_quantile(quantile)
    .check_scalar_logical(log, "log")
    cpp_dHVASIQ(x, mu, sigma, quantile, log)
}

#' @rdname HVASIQ
#' @export
pHVASIQ <- function(q, mu, sigma, quantile = 0.5,
                     lower.tail = TRUE, log.p = FALSE) {
    .check_unit_interval(q, "q", closed = TRUE)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    quantile <- .check_fixed_quantile(quantile)
    .check_scalar_logical(lower.tail, "lower.tail")
    .check_scalar_logical(log.p, "log.p")
    cpp_pHVASIQ(q, mu, sigma, quantile, lower.tail, log.p)
}

#' @rdname HVASIQ
#' @export
qHVASIQ <- function(p, mu, sigma, quantile = 0.5,
                     lower.tail = TRUE, log.p = FALSE) {
    .check_probability(p, log.p)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    quantile <- .check_fixed_quantile(quantile)
    .check_scalar_logical(lower.tail, "lower.tail")
    .check_scalar_logical(log.p, "log.p")
    cpp_qHVASIQ(p, mu, sigma, quantile, lower.tail, log.p)
}

#' @rdname HVASIQ
#' @export
rHVASIQ <- function(n, mu, sigma, quantile = 0.5) {
    n <- .n_random(n)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    quantile <- .check_fixed_quantile(quantile)
    cpp_rHVASIQ(n, mu, sigma, quantile)
}

#' @rdname HVASIQ
#' @export
HVASIQ <- function(quantile = 0.50, mu.link = "logit",
                    sigma.link = "logit") {
    quantile <- .check_fixed_quantile(quantile, "HVASIQ")

    mstats <- checklink(
        "mu.link", "HVASIQ", substitute(mu.link),
        c("logit", "probit", "cloglog", "cauchit", "log", "own")
    )
    dstats <- checklink(
        "sigma.link", "HVASIQ", substitute(sigma.link),
        c("logit", "probit", "cloglog", "cauchit", "log", "own")
    )

    structure(
        list(
            family = c("HVASIQ", "Hyperbolic-secant-kernel Vasicek-type quantile"),
            parameters = list(mu = TRUE, sigma = TRUE),
            nopar = 2,
            type = "Continuous",
            quantile = quantile,
            mu.link = as.character(substitute(mu.link)),
            sigma.link = as.character(substitute(sigma.link)),
            mu.linkfun = mstats$linkfun,
            sigma.linkfun = dstats$linkfun,
            mu.linkinv = mstats$linkinv,
            sigma.linkinv = dstats$linkinv,
            mu.dr = mstats$mu.eta,
            sigma.dr = dstats$mu.eta,
            dldm = eval(substitute(function(y, mu, sigma) {
                a <- sqrt((1 - sigma) / sigma)
                g <- pi / sin(pi * mu)
                d <- log(tan(pi * y / 2)) - log(tan(pi * mu / 2))
                z <- a * d + log(tan(pi * QUANTILE / 2))
                a * g * tanh(z)
            }, list(QUANTILE = quantile))),
            d2ldm2 = eval(substitute(function(y, mu, sigma) {
                a <- sqrt((1 - sigma) / sigma)
                g <- pi / sin(pi * mu)
                d <- log(tan(pi * y / 2)) - log(tan(pi * mu / 2))
                z <- a * d + log(tan(pi * QUANTILE / 2))
                hyperbolic_tangent <- tanh(z)
                squared_sech <- 1 / cosh(z)^2
                g^2 * (
                    -a * cos(pi * mu) * hyperbolic_tangent -
                        a^2 * squared_sech
                )
            }, list(QUANTILE = quantile))),
            dldd = eval(substitute(function(y, mu, sigma) {
                a <- sqrt((1 - sigma) / sigma)
                b <- 1 / (2 * sigma * (1 - sigma))
                d <- log(tan(pi * y / 2)) - log(tan(pi * mu / 2))
                z <- a * d + log(tan(pi * QUANTILE / 2))
                b * (a * d * tanh(z) - 1)
            }, list(QUANTILE = quantile))),
            d2ldmdd = eval(substitute(function(y, mu, sigma) {
                a <- sqrt((1 - sigma) / sigma)
                b <- 1 / (2 * sigma * (1 - sigma))
                g <- pi / sin(pi * mu)
                d <- log(tan(pi * y / 2)) - log(tan(pi * mu / 2))
                z <- a * d + log(tan(pi * QUANTILE / 2))
                hyperbolic_tangent <- tanh(z)
                squared_sech <- 1 / cosh(z)^2
                -a * b * g * (
                    hyperbolic_tangent + a * d * squared_sech
                )
            }, list(QUANTILE = quantile))),
            d2ldd2 = eval(substitute(function(y, mu, sigma) {
                a <- sqrt((1 - sigma) / sigma)
                b <- 1 / (2 * sigma * (1 - sigma))
                d <- log(tan(pi * y / 2)) - log(tan(pi * mu / 2))
                z <- a * d + log(tan(pi * QUANTILE / 2))
                hyperbolic_tangent <- tanh(z)
                squared_sech <- 1 / cosh(z)^2
                b^2 * (
                    2 * (1 - 2 * sigma) -
                        (3 - 4 * sigma) * a * d * hyperbolic_tangent -
                        a^2 * d^2 * squared_sech
                )
            }, list(QUANTILE = quantile))),
            G.dev.incr = eval(substitute(
                function(y, mu, sigma, w, ...) {
                    -2 * dHVASIQ(
                        y, mu, sigma, quantile = QUANTILE, log = TRUE
                    )
                },
                list(QUANTILE = quantile)
            )),
            rqres = as.expression(substitute(
                rqres(
                    pfun = "pHVASIQ", type = "Continuous", y = y,
                    mu = mu, sigma = sigma, quantile = QUANTILE
                ),
                list(QUANTILE = quantile)
            )),
            mu.initial = as.expression(substitute({
                mu0 <- as.numeric(
                    stats::quantile(y, probs = QUANTILE, names = FALSE)
                )
                mu <- (y + mu0) / 2
            }, list(QUANTILE = quantile))),
            sigma.initial = expression({
                sigma <- rep(0.5, length(y))
            }),
            mu.valid = function(mu) all(mu > 0 & mu < 1),
            sigma.valid = function(sigma) all(sigma > 0 & sigma < 1),
            y.valid = function(y) all(y > 0 & y < 1),
            mean = eval(substitute(function(mu, sigma) {
                n <- max(length(mu), length(sigma))
                mu <- rep_len(mu, n)
                sigma <- rep_len(sigma, n)

                vapply(seq_len(n), function(i) {
                    scale <- sqrt(sigma[i] / (1 - sigma[i]))
                    qfun <- function(p) {
                        linear_predictor <-
                            log(tan(pi * mu[i] / 2)) +
                            scale * (
                                log(tan(pi * p / 2)) -
                                log(tan(pi * QUANTILE / 2))
                            )
                        2 * atan(exp(linear_predictor)) / pi
                    }
                    stats::integrate(
                        qfun, lower = 0, upper = 1,
                        subdivisions = 200L, rel.tol = 1e-8
                    )$value
                }, numeric(1))
            }, list(QUANTILE = quantile))),
            variance = eval(substitute(function(mu, sigma) {
                n <- max(length(mu), length(sigma))
                mu <- rep_len(mu, n)
                sigma <- rep_len(sigma, n)

                vapply(seq_len(n), function(i) {
                    scale <- sqrt(sigma[i] / (1 - sigma[i]))
                    qfun <- function(p) {
                        linear_predictor <-
                            log(tan(pi * mu[i] / 2)) +
                            scale * (
                                log(tan(pi * p / 2)) -
                                log(tan(pi * QUANTILE / 2))
                            )
                        2 * atan(exp(linear_predictor)) / pi
                    }
                    expected <- stats::integrate(
                        qfun, lower = 0, upper = 1,
                        subdivisions = 200L, rel.tol = 1e-8
                    )$value
                    value <- stats::integrate(
                        function(p) (qfun(p) - expected)^2,
                        lower = 0, upper = 1,
                        subdivisions = 200L, rel.tol = 1e-8
                    )$value
                    max(value, 0)
                }, numeric(1))
            }, list(QUANTILE = quantile)))
        ),
        class = c("gamlss.family", "family")
    )
}
