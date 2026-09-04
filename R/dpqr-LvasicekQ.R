#' @importFrom gamlss gamlss
#' @importFrom gamlss.dist checklink
#' @importFrom stats plogis qlogis runif
#'
#' @name LVASIQ
#' @aliases LVASIQ dLVASIQ pLVASIQ qLVASIQ rLVASIQ
#'
#' @title Logistic-kernel Vasicek-type distribution with quantile parameterization
#'
#' @description
#' The function \code{LVASIQ()} defines the logistic-kernel Vasicek-type
#' distribution as a \code{gamlss.family} object for conditional quantile
#' regression. The functions \code{dLVASIQ}, \code{pLVASIQ},
#' \code{qLVASIQ}, and \code{rLVASIQ} give the density, distribution
#' function, quantile function, and random generation. The parameter
#' \eqn{\mu} is the conditional \eqn{\tau}-th quantile
#' (\eqn{0<\mu<1}), \eqn{\sigma} is a shape parameter
#' (\eqn{0<\sigma<1}), and \eqn{\tau\in(0,1)} is fixed by the user.
#' The fixed level is supplied through the \code{quantile} argument.
#'
#' @details
#' Let
#' \eqn{\mathrm{logit}(u)=
#' \log\left(\frac{u}{1-u}\right)} and
#' \eqn{\Lambda(z)=\frac{1}{1+e^{-z}}},
#' with
#' \eqn{\lambda(z)=
#' \Lambda(z)\left[1-\Lambda(z)\right]}.
#' Define
#' \deqn{z =
#' \sqrt{\frac{1-\sigma}{\sigma}}
#' \left[\mathrm{logit}(x)-\mathrm{logit}(\mu)\right]
#' +\mathrm{logit}(\tau).}
#'
#' Cumulative distribution function
#' \deqn{F(x\mid\mu,\sigma,\tau)=\Lambda(z).}
#'
#' Probability density function
#' \deqn{f(x\mid\mu,\sigma,\tau)=
#' \sqrt{\frac{1-\sigma}{\sigma}}
#' \frac{\lambda(z)}{x(1-x)}.}
#'
#' Quantile function
#' \deqn{Q(p\mid\mu,\sigma,\tau)=
#' \Lambda\!\left\{
#' \mathrm{logit}(\mu)
#' +\sqrt{\frac{\sigma}{1-\sigma}}
#' \left[\mathrm{logit}(p)-\mathrm{logit}(\tau)\right]
#' \right\}.}
#'
#' By construction \eqn{Q(\tau)=\mu}, i.e. \eqn{\mu} is the \eqn{\tau}-th
#' quantile. Note that, unlike the normal-kernel Vasicek distribution, the
#' logistic kernel does \strong{not} yield a closed-form mean; in particular
#' \eqn{E(X)\neq\mu} in general.
#'
#' The GAMLSS family uses analytical derivatives. For one observation, let
#' \deqn{a=\sqrt{\frac{1-\sigma}{\sigma}},\qquad
#' d=\mathrm{logit}(y)-\mathrm{logit}(\mu),\qquad
#' P=\Lambda\left\{ad+\mathrm{logit}(\tau)\right\},}
#' and define
#' \deqn{V=P(1-P),\qquad
#' b=\frac{1}{2\sigma(1-\sigma)},\qquad
#' g=\frac{1}{\mu(1-\mu)}.}
#' If \eqn{\ell} denotes the individual log-likelihood contribution, the
#' first derivatives are
#' \deqn{\frac{\partial\ell}{\partial\mu}
#' =-ag(1-2P)}
#' and
#' \deqn{\frac{\partial\ell}{\partial\sigma}
#' =-b\left\{1+ad(1-2P)\right\}.}
#' The second and cross derivatives are
#' \deqn{\frac{\partial^2\ell}{\partial\mu^2}
#' =g^2\left\{
#' a(1-2\mu)(1-2P)-2a^2V
#' \right\},}
#' \deqn{\frac{\partial^2\ell}{\partial\mu\,\partial\sigma}
#' =abg\left\{
#' (1-2P)-2adV
#' \right\},}
#' and
#' \deqn{\frac{\partial^2\ell}{\partial\sigma^2}
#' =b^2\left\{
#' 2(1-2\sigma)
#' +ad(3-4\sigma)(1-2P)
#' -2a^2d^2V
#' \right\}.}
#' These expressions are evaluated directly by \code{LVASIQ()}; numerical
#' differentiation is not used. The \code{mean} and \code{variance}
#' components of the family object use numerical quadrature because the
#' corresponding moments do not have elementary closed forms.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#'
#' @references
#' Mazucheli, J., Alves, B., Korkmaz, M. C. and Leiva, V. (2022).
#' Vasicek quantile and mean regression models for bounded data:
#' New formulation, mathematical derivations, and numerical applications.
#' \emph{Mathematics}, \bold{10}, 1389.
#'
#' Vasicek, O. A. (2002). The distribution of loan portfolio value.
#' \emph{Risk}, \bold{15}(12), 160--162.
#'
#' Witzany, J. (2013). A note on the Vasicek's model with the logistic
#' distribution. \emph{Ekonomický časopis (Journal of Economics)},
#' \bold{61}(10), 1053--1066.
#'
#' @param x Vector of quantiles in \eqn{(0,1)}.
#' @param q Vector of values in \eqn{[0,1]} at which the cumulative
#'   distribution function is evaluated.
#' @param p Vector of probabilities in \eqn{[0,1]} on the probability scale.
#' @param n Number of observations.
#' @param mu Vector of \eqn{\tau}-quantiles, \eqn{0<\mu<1}.
#' @param sigma Vector of shape parameter values, \eqn{0<\sigma<1}.
#' @param quantile Fixed quantile level \eqn{\tau\in(0,1)} represented by
#'   \eqn{\mu}, used in the distribution functions and in \code{LVASIQ()}.
#' @param mu.link Link function for the \eqn{\mu} parameter.
#' @param sigma.link Link function for the \eqn{\sigma} parameter.
#' @param lower.tail Logical; if \code{TRUE}, probabilities are \eqn{P(X\le x)}.
#' @param log Logical; if \code{TRUE}, the log-density is returned.
#' @param log.p Logical; if \code{TRUE}, probabilities \code{p} are given
#' as \code{log(p)} or cumulative probabilities are returned on the log
#' scale, as appropriate.
#'
#' @return
#' \code{dLVASIQ} gives the density, \code{pLVASIQ} the distribution function,
#' \code{qLVASIQ} the quantile function, and \code{rLVASIQ} generates random
#' deviates. \code{LVASIQ()} returns a \code{gamlss.family} object.
#'
#' @note
#' The level supplied through \code{quantile} is stored in the family
#' definition and embedded as a numeric literal in the family components
#' used by GAMLSS; no global variable is required.
#'
#' @examples
#' set.seed(123)
#' x <- rLVASIQ(n = 1000, mu = 0.50, sigma = 0.25, quantile = 0.5)
#' S <- seq(min(x), max(x), length.out = 1000)
#'
#' hist(x, prob = TRUE, main = "Logistic-kernel Vasicek-type")
#' lines(S, dLVASIQ(x = S, mu = 0.50, sigma = 0.25, quantile = 0.5), col = 2)
#'
#' plot(ecdf(x))
#' lines(S, pLVASIQ(q = S, mu = 0.50, sigma = 0.25, quantile = 0.5), col = 2)
#'
#' data <- data.frame(
#'   y = rLVASIQ(n = 100, mu = 0.50, sigma = 0.25, quantile = 0.50)
#' )
#' fit <- gamlss::gamlss(
#'     y ~ 1,
#'     data = data,
#'     family = LVASIQ(
#'       quantile = 0.50, mu.link = "logit", sigma.link = "logit"
#'     )
#' )
#' fitted(fit, what = "mu")[1:5]
#'
NULL

##################################################
#' @rdname LVASIQ
#' @export
dLVASIQ <- function(x, mu, sigma, quantile = 0.5, log = FALSE) {
    .check_unit_interval(x, "x")
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    quantile <- .check_fixed_quantile(quantile)
    .check_scalar_logical(log, "log")
    cpp_dLVASIQ(x, mu, sigma, quantile, log)
}

##################################################
#' @rdname LVASIQ
#' @export
pLVASIQ <- function(q, mu, sigma, quantile = 0.5, lower.tail = TRUE, log.p = FALSE) {
    .check_unit_interval(q, "q", closed = TRUE)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    quantile <- .check_fixed_quantile(quantile)
    .check_scalar_logical(lower.tail, "lower.tail")
    .check_scalar_logical(log.p, "log.p")
    cpp_pLVASIQ(q, mu, sigma, quantile, lower.tail, log.p)
}

##################################################
#' @rdname LVASIQ
#' @export
qLVASIQ <- function(p, mu, sigma, quantile = 0.5, lower.tail = TRUE, log.p = FALSE) {
    .check_probability(p, log.p)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    quantile <- .check_fixed_quantile(quantile)
    .check_scalar_logical(lower.tail, "lower.tail")
    .check_scalar_logical(log.p, "log.p")
    cpp_qLVASIQ(p, mu, sigma, quantile, lower.tail, log.p)
}

##################################################
#' @rdname LVASIQ
#' @export
rLVASIQ <- function(n, mu, sigma, quantile = 0.5) {
    n <- .n_random(n)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    quantile <- .check_fixed_quantile(quantile)
    cpp_rLVASIQ(n, mu, sigma, quantile)
}

##################################################
#' @rdname LVASIQ
#' @export
LVASIQ <- function(quantile = 0.50, mu.link = "logit",
                   sigma.link = "logit") {
    quantile <- .check_fixed_quantile(quantile, "LVASIQ")

    mstats <- checklink(
        "mu.link", "LVASIQ", substitute(mu.link),
        c("logit", "probit", "cloglog", "cauchit", "log", "own")
    )
    dstats <- checklink(
        "sigma.link", "LVASIQ", substitute(sigma.link),
        c("logit", "probit", "cloglog", "cauchit", "log", "own")
    )

    structure(
        list(
            family = c("LVASIQ", "Logistic-kernel Vasicek-type quantile"),
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
                g <- 1 / (mu * (1 - mu))
                d <- stats::qlogis(y) - stats::qlogis(mu)
                z <- a * d + stats::qlogis(QUANTILE)
                prob <- stats::plogis(z)

                -a * g * (1 - 2 * prob)
            }, list(QUANTILE = quantile))),
            d2ldm2 = eval(substitute(function(y, mu, sigma) {
                a <- sqrt((1 - sigma) / sigma)
                g <- 1 / (mu * (1 - mu))
                d <- stats::qlogis(y) - stats::qlogis(mu)
                z <- a * d + stats::qlogis(QUANTILE)
                prob <- stats::plogis(z)
                kernel <- prob * (1 - prob)

                g^2 * (
                    a * (1 - 2 * mu) * (1 - 2 * prob) -
                        2 * a^2 * kernel
                )
            }, list(QUANTILE = quantile))),
            dldd = eval(substitute(function(y, mu, sigma) {
                a <- sqrt((1 - sigma) / sigma)
                b <- 1 / (2 * sigma * (1 - sigma))
                d <- stats::qlogis(y) - stats::qlogis(mu)
                z <- a * d + stats::qlogis(QUANTILE)
                prob <- stats::plogis(z)

                -b * (1 + a * d * (1 - 2 * prob))
            }, list(QUANTILE = quantile))),
            d2ldmdd = eval(substitute(function(y, mu, sigma) {
                a <- sqrt((1 - sigma) / sigma)
                b <- 1 / (2 * sigma * (1 - sigma))
                g <- 1 / (mu * (1 - mu))
                d <- stats::qlogis(y) - stats::qlogis(mu)
                z <- a * d + stats::qlogis(QUANTILE)
                prob <- stats::plogis(z)
                kernel <- prob * (1 - prob)

                a * b * g * (
                    (1 - 2 * prob) - 2 * a * d * kernel
                )
            }, list(QUANTILE = quantile))),
            d2ldd2 = eval(substitute(function(y, mu, sigma) {
                a <- sqrt((1 - sigma) / sigma)
                b <- 1 / (2 * sigma * (1 - sigma))
                d <- stats::qlogis(y) - stats::qlogis(mu)
                z <- a * d + stats::qlogis(QUANTILE)
                prob <- stats::plogis(z)
                kernel <- prob * (1 - prob)

                b^2 * (
                    2 * (1 - 2 * sigma) +
                        a * d * (3 - 4 * sigma) * (1 - 2 * prob) -
                        2 * a^2 * d^2 * kernel
                )
            }, list(QUANTILE = quantile))),
            G.dev.incr = eval(substitute(
                function(y, mu, sigma, w, ...) {
                    -2 * dLVASIQ(
                        y, mu, sigma, quantile = QUANTILE, log = TRUE
                    )
                },
                list(QUANTILE = quantile)
            )),
            rqres = as.expression(substitute(
                rqres(
                    pfun = "pLVASIQ", type = "Continuous", y = y,
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
                        stats::plogis(
                            stats::qlogis(mu[i]) +
                                scale * (
                                    stats::qlogis(p) - stats::qlogis(QUANTILE)
                                )
                        )
                    }
                    stats::integrate(
                        qfun,
                        lower = 0,
                        upper = 1,
                        subdivisions = 200L,
                        rel.tol = 1e-8
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
                        stats::plogis(
                            stats::qlogis(mu[i]) +
                                scale * (
                                    stats::qlogis(p) - stats::qlogis(QUANTILE)
                                )
                        )
                    }
                    expected <- stats::integrate(
                        qfun,
                        lower = 0,
                        upper = 1,
                        subdivisions = 200L,
                        rel.tol = 1e-8
                    )$value
                    variance <- stats::integrate(
                        function(p) (qfun(p) - expected)^2,
                        lower = 0,
                        upper = 1,
                        subdivisions = 200L,
                        rel.tol = 1e-8
                    )$value

                    max(variance, 0)
                }, numeric(1))
            }, list(QUANTILE = quantile)))
        ),
        class = c("gamlss.family", "family")
    )
}
