#' Hyperbolic-secant-kernel Vasicek-type quantile distribution
#'
#' @description
#' \code{HSVASIQ()} defines the hyperbolic-secant-kernel Vasicek-type distribution
#' as a \code{gamlss.family} object for conditional quantile regression.
#' The functions \code{dHSVASIQ()}, \code{pHSVASIQ()}, \code{qHSVASIQ()}, and
#' \code{rHSVASIQ()} provide the density, distribution function, quantile
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
#' For GAMLSS fitting, \code{tau} must be defined as a scalar variable in the
#' global environment before \code{HSVASIQ()} is evaluated. It must retain the
#' same value when residuals or other post-fit quantities are computed.
#'
#' @param x,q Vector of quantiles in \eqn{(0,1)}.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param mu Vector of conditional \eqn{\tau}-th quantiles,
#'   \eqn{0<\mu<1}.
#' @param sigma Vector of shape parameter values, \eqn{0<\sigma<1}.
#' @param tau Scalar in \eqn{(0,1)} fixing the quantile represented by
#'   \eqn{\mu}. For \code{HSVASIQ()}, it must be defined globally.
#' @param mu.link Link function for \eqn{\mu}.
#' @param sigma.link Link function for \eqn{\sigma}.
#' @param lower.tail Logical; if \code{TRUE}, probabilities are
#'   \eqn{P(X\leq x)}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @param log.p Logical; if \code{TRUE}, probabilities are supplied or
#'   returned on the logarithmic scale.
#'
#' @return
#' \code{dHSVASIQ()} returns the density, \code{pHSVASIQ()} the distribution
#' function, \code{qHSVASIQ()} the quantile function, and \code{rHSVASIQ()}
#' random deviates. \code{HSVASIQ()} returns a \code{gamlss.family} object.
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
#' y <- rHSVASIQ(500, mu = 0.60, sigma = 0.30, tau = 0.25)
#'
#' tau <- 0.25
#' fit <- gamlss::gamlss(
#'     y ~ 1,
#'     sigma.formula = ~ 1,
#'     family = HSVASIQ(),
#'     control = gamlss::gamlss.control(trace = FALSE)
#' )
#' fitted(fit, what = "mu")[1]
#' rm(tau)
#'
#' @name HSVASIQ
#' @aliases HSVASIQ dHSVASIQ pHSVASIQ qHSVASIQ rHSVASIQ
#' @importFrom gamlss gamlss
#' @importFrom gamlss.dist checklink
#' @export
dHSVASIQ <- function(x, mu, sigma, tau = 0.5, log = FALSE) {
    .check_unit_interval(x, "x")
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    .check_unit_interval(tau, "tau")
    if (length(tau) != 1L) {
        stop("'tau' must be a single number in (0, 1).", call. = FALSE)
    }
    .check_scalar_logical(log, "log")
    cpp_dHSVASIQ(x, mu, sigma, tau, log)
}

#' @rdname HSVASIQ
#' @export
pHSVASIQ <- function(q, mu, sigma, tau = 0.5,
                     lower.tail = TRUE, log.p = FALSE) {
    .check_unit_interval(q, "q", closed = TRUE)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    .check_unit_interval(tau, "tau")
    if (length(tau) != 1L) {
        stop("'tau' must be a single number in (0, 1).", call. = FALSE)
    }
    .check_scalar_logical(lower.tail, "lower.tail")
    .check_scalar_logical(log.p, "log.p")
    cpp_pHSVASIQ(q, mu, sigma, tau, lower.tail, log.p)
}

#' @rdname HSVASIQ
#' @export
qHSVASIQ <- function(p, mu, sigma, tau = 0.5,
                     lower.tail = TRUE, log.p = FALSE) {
    .check_probability(p, log.p)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    .check_unit_interval(tau, "tau")
    if (length(tau) != 1L) {
        stop("'tau' must be a single number in (0, 1).", call. = FALSE)
    }
    .check_scalar_logical(lower.tail, "lower.tail")
    .check_scalar_logical(log.p, "log.p")
    cpp_qHSVASIQ(p, mu, sigma, tau, lower.tail, log.p)
}

#' @rdname HSVASIQ
#' @export
rHSVASIQ <- function(n, mu, sigma, tau = 0.5) {
    n <- .n_random(n)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    .check_unit_interval(tau, "tau")
    if (length(tau) != 1L) {
        stop("'tau' must be a single number in (0, 1).", call. = FALSE)
    }
    cpp_rHSVASIQ(n, mu, sigma, tau)
}

#' @rdname HSVASIQ
#' @export
HSVASIQ <- function(mu.link = "logit", sigma.link = "logit") {
    if (!exists("tau", envir = .GlobalEnv, inherits = FALSE)) {
        stop(
            "For HSVASIQ(), define a global scalar 'tau' in (0, 1).",
            call. = FALSE
        )
    }
    tau <- get("tau", envir = .GlobalEnv, inherits = FALSE)
    if (!is.numeric(tau) || length(tau) != 1L ||
        is.na(tau) || !is.finite(tau) || tau <= 0 || tau >= 1) {
        stop(
            paste(
                "For HSVASIQ(), global 'tau' must be a single number",
                "that is finite and strictly between 0 and 1."
            ),
            call. = FALSE
        )
    }
    tau <- as.numeric(tau)

    mstats <- checklink(
        "mu.link", "HSVASIQ", substitute(mu.link),
        c("logit", "probit", "cloglog", "cauchit", "log", "own")
    )
    dstats <- checklink(
        "sigma.link", "HSVASIQ", substitute(sigma.link),
        c("logit", "probit", "cloglog", "cauchit", "log", "own")
    )

    structure(
        list(
            family = c("HSVASIQ", "Hyperbolic-secant-kernel Vasicek-type quantile"),
            parameters = list(mu = TRUE, sigma = TRUE),
            nopar = 2,
            type = "Continuous",
            mu.link = as.character(substitute(mu.link)),
            sigma.link = as.character(substitute(sigma.link)),
            mu.linkfun = mstats$linkfun,
            sigma.linkfun = dstats$linkfun,
            mu.linkinv = mstats$linkinv,
            sigma.linkinv = dstats$linkinv,
            mu.dr = mstats$mu.eta,
            sigma.dr = dstats$mu.eta,
            dldm = function(y, mu, sigma) {
                a <- sqrt((1 - sigma) / sigma)
                g <- pi / sin(pi * mu)
                d <- log(tan(pi * y / 2)) - log(tan(pi * mu / 2))
                z <- a * d + log(tan(pi * tau / 2))
                a * g * tanh(z)
            },
            d2ldm2 = function(y, mu, sigma) {
                a <- sqrt((1 - sigma) / sigma)
                g <- pi / sin(pi * mu)
                d <- log(tan(pi * y / 2)) - log(tan(pi * mu / 2))
                z <- a * d + log(tan(pi * tau / 2))
                hyperbolic_tangent <- tanh(z)
                squared_sech <- 1 / cosh(z)^2
                g^2 * (
                    -a * cos(pi * mu) * hyperbolic_tangent -
                        a^2 * squared_sech
                )
            },
            dldd = function(y, mu, sigma) {
                a <- sqrt((1 - sigma) / sigma)
                b <- 1 / (2 * sigma * (1 - sigma))
                d <- log(tan(pi * y / 2)) - log(tan(pi * mu / 2))
                z <- a * d + log(tan(pi * tau / 2))
                b * (a * d * tanh(z) - 1)
            },
            d2ldmdd = function(y, mu, sigma) {
                a <- sqrt((1 - sigma) / sigma)
                b <- 1 / (2 * sigma * (1 - sigma))
                g <- pi / sin(pi * mu)
                d <- log(tan(pi * y / 2)) - log(tan(pi * mu / 2))
                z <- a * d + log(tan(pi * tau / 2))
                hyperbolic_tangent <- tanh(z)
                squared_sech <- 1 / cosh(z)^2
                -a * b * g * (
                    hyperbolic_tangent + a * d * squared_sech
                )
            },
            d2ldd2 = function(y, mu, sigma) {
                a <- sqrt((1 - sigma) / sigma)
                b <- 1 / (2 * sigma * (1 - sigma))
                d <- log(tan(pi * y / 2)) - log(tan(pi * mu / 2))
                z <- a * d + log(tan(pi * tau / 2))
                hyperbolic_tangent <- tanh(z)
                squared_sech <- 1 / cosh(z)^2
                b^2 * (
                    2 * (1 - 2 * sigma) -
                        (3 - 4 * sigma) * a * d * hyperbolic_tangent -
                        a^2 * d^2 * squared_sech
                )
            },
            G.dev.incr = function(y, mu, sigma, w, ...) {
                -2 * dHSVASIQ(y, mu, sigma, tau, log = TRUE)
            },
            rqres = expression(
                rqres(
                    pfun = "pHSVASIQ", type = "Continuous", y = y,
                    mu = mu, sigma = sigma, tau = tau
                )
            ),
            mu.initial = expression({
                mu0 <- as.numeric(
                    stats::quantile(y, probs = tau, names = FALSE)
                )
                mu <- (y + mu0) / 2
            }),
            sigma.initial = expression({
                sigma <- rep(0.5, length(y))
            }),
            mu.valid = function(mu) all(mu > 0 & mu < 1),
            sigma.valid = function(sigma) all(sigma > 0 & sigma < 1),
            y.valid = function(y) all(y > 0 & y < 1),
            mean = function(mu, sigma) {
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
                                log(tan(pi * tau / 2))
                            )
                        2 * atan(exp(linear_predictor)) / pi
                    }
                    stats::integrate(
                        qfun, lower = 0, upper = 1,
                        subdivisions = 200L, rel.tol = 1e-8
                    )$value
                }, numeric(1))
            },
            variance = function(mu, sigma) {
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
                                log(tan(pi * tau / 2))
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
            }
        ),
        class = c("gamlss.family", "family")
    )
}
