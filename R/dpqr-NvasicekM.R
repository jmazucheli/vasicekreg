#' @importFrom gamlss gamlss
#' @importFrom gamlss.dist checklink
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats dnorm qnorm runif
#'
#' @name NVASIM
#' @aliases dNVASIM pNVASIM qNVASIM rNVASIM
#'
#' @title Normal-kernel Vasicek-type distribution with mean parameterization
#'
#' @description
#' Defines the normal-kernel Vasicek distribution under a mean
#' parameterization for use as a \code{gamlss.family}. The parameter
#' \eqn{\mu} represents the mean of the distribution, with
#' \eqn{0 < \mu < 1}, and \eqn{\sigma} is a shape parameter.
#'
#' The density, distribution function, quantile function and random
#' number generation are provided by \code{dNVASIM()}, \code{pNVASIM()},
#' \code{qNVASIM()} and \code{rNVASIM()}, respectively.
#'
#' @author
#' Josmar Mazucheli \email{jmazucheli@gmail.com}
#' Bruna Alves \email{pg402900@uem.br}
#'
#' @references
#' Hastie, T. J. and Tibshirani, R. J. (1990).
#' \emph{Generalized Additive Models}. Chapman and Hall, London.
#'
#' Mazucheli, J., Alves, B., Korkmaz, M. C., and Leiva, V. (2022).
#' Vasicek quantile and mean regression models for bounded data:
#' New formulation, mathematical derivations, and numerical applications.
#' \emph{Mathematics}, \bold{10}, 1389.
#' \doi{10.3390/math10091389}
#'
#' Rigby, R. A. and Stasinopoulos, D. M. (2005).
#' Generalized additive models for location, scale and shape
#' (with discussion). \emph{Applied Statistics}, \bold{54}(3), 507--554.
#'
#' Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and
#' De Bastiani, F. (2019).
#' \emph{Distributions for Modeling Location, Scale, and Shape:
#' Using GAMLSS in R}. Chapman and Hall/CRC.
#'
#' Stasinopoulos, D. M. and Rigby, R. A. (2007).
#' Generalized additive models for location, scale and shape (GAMLSS)
#' in R. \emph{Journal of Statistical Software}, \bold{23}(7), 1--46.
#'
#' Stasinopoulos, D. M., Rigby, R. A., Heller, G., Voudouris, V.,
#' and De Bastiani, F. (2017).
#' \emph{Flexible Regression and Smoothing: Using GAMLSS in R}.
#' Chapman and Hall/CRC.
#'
#' Vasicek, O. A. (1987).
#' Probability of loss on loan portfolio. KMV Corporation.
#'
#' Vasicek, O. A. (2002).
#' The distribution of loan portfolio value.
#' \emph{Risk}, \bold{15}(12), 160--162.
#'
#' @param x Vector of quantiles in the interval \eqn{(0,1)}.
#' @param q Vector of values in \eqn{[0,1]} at which the cumulative
#'   distribution function is evaluated.
#' @param p Vector of probabilities in \eqn{[0,1]} on the probability scale.
#' @param n Number of observations. If \code{length(n) > 1}, the length
#' is taken to be the number required.
#' @param log Logical; if \code{TRUE}, the log-density is returned.
#' @param log.p Logical; if \code{TRUE}, probabilities \code{p} are given
#' as \code{log(p)} or cumulative probabilities are returned on the log
#' scale, as appropriate.
#' @param lower.tail Logical; if \code{TRUE},
#' probabilities \eqn{P(X \le x)} are returned.
#' @param mu.link Link function for the \eqn{\mu} parameter.
#' @param sigma.link Link function for the \eqn{\sigma} parameter.
#' @param mu Vector of mean values.
#' @param sigma Vector of shape parameter values.
#'
#' @return
#' \code{NVASIM()} returns a \code{gamlss.family} object.
#'
#' @note
#' In the \code{NVASIM()} parameterization, \eqn{\mu} corresponds to the
#' mean of the distribution and \eqn{\sigma} is a shape parameter.
#'
#' @seealso \code{\link[vasicekreg]{NVASIQ}}, \code{\link[mvtnorm]{pmvnorm}}
#'
#' @details
#' Probability density function 
#' \deqn{f(x\mid \mu ,\sigma )=\sqrt{\frac{1-\sigma }{\sigma }}\exp \left\{ \frac{1}{2}\left[ \Phi ^{-1}\left( x\right) ^{2}-\left( \frac{\Phi ^{-1}\left(  x\right)    \sqrt{1-\sigma }-\Phi ^{-1}\left( \mu \right) }{\sqrt{\sigma }}\right) ^{2}\right] \right\}}
#' 
#' Cumulative distribution function
#' \deqn{F(x\mid \mu ,\sigma )=\Phi \left( \frac{\Phi ^{-1}\left( x\right) \sqrt{1-\sigma }-\Phi ^{-1}\left( \mu \right) }{\sqrt{\sigma }}\right)}
#' 
#' Quantile function
#' \deqn{Q(p \mid \mu ,\sigma )=F^{-1}(p \mid \mu ,\sigma )=\Phi \left(\frac{\Phi ^{-1}\left(\mu\right) +\Phi ^{-1}\left( p \right) \sqrt{\sigma }}{\sqrt{1-\sigma }}\right) }
#' 
#' Expected value
#' \deqn{E(X) = \mu} 
#' 
#' Variance
#' \deqn{Var(X) = \Phi_2\left ( \Phi^{-1}(\mu),\Phi^{-1}(\mu),\sigma \right )-\mu^2}
#' where \eqn{(x, \mu, \sigma, p) \in (0,1)} and \eqn{\Phi_2(\cdot)} is the cumulative distribution function for the standard bivariate normal distribution with correlation \eqn{\sigma}.
#'
#' @examples
#'
#' set.seed(123)
#' x <- rNVASIM(n = 1000, mu = 0.50, sigma = 0.69)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], length.out = 1000)
#' 
#' hist(x, prob = TRUE, main = 'Vasicek')
#' lines(S, dNVASIM(x = S, mu = 0.50, sigma = 0.69), col = 2)
#' 
#' plot(ecdf(x))
#' lines(S, pNVASIM(q = S, mu = 0.50, sigma = 0.69), col = 2)
#' 
#' plot(quantile(x, probs = S), type = "l")
#' lines(qNVASIM(p = S, mu = 0.50, sigma = 0.69), col = 2)
#' 
#' library(gamlss)
#' set.seed(123)
#' data <- data.frame(y =  rNVASIM(n = 100, mu = 0.5, sigma = 0.69))
#' 
#' fit <- gamlss(y ~ 1, data = data, mu.link = 'logit', sigma.link = 'logit', family = NVASIM)
#' 1 /(1 + exp(-fit$mu.coefficients))
#' 1 /(1 + exp(-fit$sigma.coefficients))
#'
#' \dontrun{
#' library(gamlss)
#' set.seed(123)
#' 
#' n <- 1000
#' x <- rbinom(n, size = 1, prob = 0.5)
#' eta <- 0.5 + 1 * x;
#' mu <- 1 / (1 + exp(-eta));
#' sigma <- 0.5;
#' y <- rNVASIM(n, mu, sigma)
#' data <- data.frame(y, x)
#' 
#' fit <- gamlss(y ~ x, data = data, family = NVASIM, mu.link = 'logit', sigma.link = 'logit', 
#' control = gamlss.control(n.cyc = 200))
#' summary(fit)
#' }
##################################################
#' @rdname NVASIM
#' @export
dNVASIM <- function(x, mu, sigma, log = FALSE)
{
    .check_unit_interval(x, "x")
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    .check_scalar_logical(log, "log")
    cpp_dNVASIM(x, mu, sigma, log)
}
##################################################
#' @rdname NVASIM
#' @export
pNVASIM <- function(q, mu, sigma,
                   lower.tail = TRUE, log.p = FALSE)
{
    .check_unit_interval(q, "q", closed = TRUE)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    .check_scalar_logical(lower.tail, "lower.tail")
    .check_scalar_logical(log.p, "log.p")
    cpp_pNVASIM(q, mu, sigma,
                     lower.tail, log.p)
}
##################################################
#' @rdname NVASIM
#' @export
qNVASIM <- function(p, mu, sigma,
                   lower.tail = TRUE, log.p = FALSE)
{
    .check_probability(p, log.p)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    .check_scalar_logical(lower.tail, "lower.tail")
    .check_scalar_logical(log.p, "log.p")
    cpp_qNVASIM(p, mu, sigma,
                     lower.tail, log.p)
}
##################################################
#' @rdname NVASIM
#' @export
rNVASIM <- function(n, mu, sigma)
{
    n <- .n_random(n)
    .check_unit_interval(mu, "mu")
    .check_unit_interval(sigma, "sigma")
    cpp_qNVASIM(runif(n), mu, sigma, TRUE, FALSE)
}
##################################################
#' @rdname NVASIM
#' @export
#' 
NVASIM <- function (mu.link = "logit", sigma.link = "logit")
{
    mstats <- checklink("mu.link", "NVASIM", substitute(mu.link), c("logit", "probit", "cloglog", "cauchit", "log", "own"))
    dstats <- checklink("sigma.link", "NVASIM", substitute(sigma.link), c("logit", "probit", "cloglog", "cauchit", "log", "own"))
    structure(list(family = c("NVASIM", "Normal-kernel Vasicek-type mean"),
                   parameters = list(mu = TRUE, sigma = TRUE), nopar = 2, type = "Continuous", mu.link = as.character(substitute(mu.link)),
                   sigma.link = as.character(substitute(sigma.link)), mu.linkfun = mstats$linkfun, sigma.linkfun = dstats$linkfun, mu.linkinv = mstats$linkinv,
                   sigma.linkinv = dstats$linkinv, mu.dr = mstats$mu.eta, sigma.dr = dstats$mu.eta,
                   dldm = function(y, mu, sigma){
                       t2 <- sqrt(0.1e1 - sigma);
                       t4 <- qnorm(mu);
                       t8 <- 0.1e1 / dnorm(t4);
                       qnormx <- qnorm(y);
                       return(0.10e1 * (qnormx * t2 - t4) / sigma * t8);
                   },
                   d2ldm2 = function(y, mu, sigma){
                       t9 <- qnorm(mu);
                       t1 <- 0.1e1 / dnorm(t9);
                       t2 <- t1 * t1;
                       t3 <- 0.1e1 / sigma;
                       t7 <- sqrt(0.1e1 - sigma);
                       t12 <-  t9 * (t1 ^ 0.2e1);
                       qnormx <- qnorm(y);
                       return(-0.10e1 * t2 * t3 + 0.10e1 * (qnormx * t7 - t9) * t3 * t12);
                   },
                   dldd = function(y, mu, sigma){
                       qnormx <- qnorm(y);
                       t1 <- 0.1e1 - sigma;
                       t4 <- 0.1e1 / sigma;
                       t6 <- sqrt(t1);
                       t8 <- qnorm(mu);
                       t9 <- qnormx * t6 - t8;
                       t15 <- t9 * t9;
                       t16 <- sigma * sigma;
                       return(-0.5e0 / t1 - 0.5e0 * t4 + 0.5e0 * t9 * t4 * qnormx / t6 + 0.5e0 * t15 / t16);
                   },
                   d2ldd2 = function(y, mu, sigma){
                       qnormx <- qnorm(y);
                       t1 <- 0.1e1 - sigma;
                       t2 <- t1 * t1;
                       t5 <- sigma * sigma;
                       t6 <- 0.1e1 / t5;
                       t8 <- qnormx * qnormx;
                       t11 <- 0.1e1 / sigma;
                       t14 <- sqrt(t1);
                       t16 <- qnorm(mu);
                       t17 <- qnormx * t14 - t16;
                       t29 <- t17 * t17;
                       return(-0.5e0 / t2 + 0.5e0 * t6 - 0.2500000000e0 * t8 / t1 * t11 - 0.10e1 * t17 * t6 * qnormx / t14 + 
                                  0.2500000000e0 * t17 * t11 * qnormx / t14 / t1 - 0.10e1 * t29 / t5 / sigma);
                   },
                   d2ldmdd = function(y, mu, sigma){
                       t1 <- qnorm(mu);
                       t2 <- sqrt(0.1e1 - sigma);
                       t6 <- 0.1e1 / dnorm(t1)
                       t13 <- sigma * sigma;
                       qnormx <- qnorm(y);
                       return(-0.5000000000e0 * qnormx / t2 / sigma * t6 - 0.10e1 * (qnormx * t2 - t1) / t13 * t6);
                   },
                   G.dev.incr = function(y, mu, sigma, w, ...) -2 * dNVASIM(y, mu, sigma, log = TRUE),
                   rqres = expression(rqres(pfun = "pNVASIM", type = "Continuous", y = y, mu = mu, sigma = sigma)),
                   mu.initial = expression({mu <- (y + mean(y))/2}),
                   sigma.initial = expression({sigma <- rep(0.5, length(y))}),
                   mu.valid = function(mu) all(mu > 0 & mu < 1),
                   sigma.valid = function(sigma) all(sigma > 0 & sigma < 1),
                   y.valid = function(y) all(y > 0 & y < 1),
                   mean = function(mu) mu,
                   variance = function(mu, sigma) {
                       n <- max(length(mu), length(sigma))
                       mu <- rep_len(mu, n)
                       sigma <- rep_len(sigma, n)
                       vapply(seq_len(n), function(i) {
                           second_moment <- pmvnorm(
                               lower = c(-Inf, -Inf),
                               upper = rep(qnorm(mu[i]), 2),
                               mean = c(0, 0),
                               corr = matrix(c(1, sigma[i], sigma[i], 1),
                                             ncol = 2)
                           )
                           as.numeric(second_moment) - mu[i]^2
                       }, numeric(1))
                   }),
              class = c("gamlss.family", "family"))
}
