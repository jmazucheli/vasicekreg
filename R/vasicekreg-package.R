#' @title Overview of the vasicekreg package
#'
#' @description
#' The \pkg{vasicekreg} package provides distribution functions and GAMLSS
#' families for Vasicek-type distributions on the unit interval. The package
#' includes three parameterizations:
#' \itemize{
#'   \item \code{NVASIM}: normal kernel with mean parameterization, where
#'   \eqn{\mu=E(Y)}.
#'   \item \code{NVASIQ}: normal kernel with quantile parameterization, where
#'   \eqn{\mu=Q_Y(\tau)} for a fixed \eqn{\tau\in(0,1)}.
#'   \item \code{LVASIQ}: logistic kernel with quantile parameterization, where
#'   \eqn{\mu=Q_Y(\tau)} for a fixed \eqn{\tau\in(0,1)}.
#' }
#' In each case, \eqn{\sigma\in(0,1)} controls dispersion. The corresponding
#' \code{d}, \code{p}, \code{q}, and \code{r} functions provide the density,
#' cumulative distribution function, quantile function, and random number
#' generation, respectively.
#'
#' @details
#' \code{\link[vasicekreg]{bodyfat}}:
#' Body fat dataset.
#'
#' \code{\link[vasicekreg]{NVASIM}}:
#' Normal-kernel mean parameterization and GAMLSS family. In regression
#' models, covariates describe the conditional mean through \eqn{\mu}.
#'
#' \code{\link[vasicekreg]{NVASIQ}}:
#' Normal-kernel quantile parameterization and GAMLSS family. For a fixed
#' quantile level \eqn{\tau}, covariates describe the conditional
#' \eqn{\tau}-th quantile through \eqn{\mu}.
#'
#' \code{\link[vasicekreg]{LVASIQ}}:
#' Logistic-kernel quantile parameterization and GAMLSS family. For a fixed
#' quantile level \eqn{\tau}, covariates describe the conditional
#' \eqn{\tau}-th quantile through \eqn{\mu}. A logistic-kernel mean-regression
#' family is not provided because the mean has no closed-form expression and
#' does not equal \eqn{\mu} under this parameterization.
#'
#' The distribution functions \code{dNVASIM}, \code{pNVASIM},
#' \code{qNVASIM}, \code{dNVASIQ}, \code{pNVASIQ}, \code{qNVASIQ},
#' \code{dLVASIQ}, \code{pLVASIQ}, and \code{qLVASIQ} call compiled
#' \proglang{C++} routines through \pkg{Rcpp}. The functions
#' \code{rNVASIM} and \code{rNVASIQ} use inverse-transform generation with
#' their compiled quantile routines, whereas \code{rLVASIQ} calls a compiled
#' random-generation routine directly. Parameter validation, the GAMLSS
#' family definitions, and the analytical log-likelihood derivatives are
#' implemented in \proglang{R}. The mean and variance components of the
#' \code{LVASIQ()} family object are obtained by numerical quadrature because
#' these moments have no closed-form expressions.
#'
#' For the distribution functions associated with \code{NVASIQ} and
#' \code{LVASIQ}, \code{tau} is supplied as an argument. For GAMLSS fitting,
#' \code{NVASIQ()} and \code{LVASIQ()} require \code{tau} to be defined as a
#' scalar variable in the global environment. The same value must be retained
#' when residuals or other post-fit quantities are computed.
#'
#' @author
#' Josmar Mazucheli \email{jmazucheli@gmail.com}
#'
#' Bruna Alves \email{pg402900@uem.br}
#'
#' @useDynLib vasicekreg
#' @importFrom Rcpp sourceCpp
"_PACKAGE"

.onUnload <- function(libpath) {
    library.dynam.unload("vasicekreg", libpath)
}
