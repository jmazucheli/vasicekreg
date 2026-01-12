#' @title Overview of the vasicekreg package
#'
#' @description
#' The \pkg{vasicekreg} package implements probability density, cumulative
#' distribution, quantile, and random number generation functions for the
#' Vasicek distribution parameterized either by its mean or by its
#' \eqn{\tau}-th quantile, with \eqn{0 < \tau < 1}. In addition, two GAMLSS
#' frameworks for regression analysis are provided. Some functions are
#' written in \proglang{C++} using \pkg{Rcpp}.
#'
#' @details
#' \code{\link[vasicekreg]{bodyfat}}: Body fat dataset.
#'
#' \code{\link[vasicekreg]{VASIM}}: Mean modeling conditional or
#' unconditional on covariates.
#'
#' \code{\link[vasicekreg]{VASIQ}}: Quantile modeling conditional or
#' unconditional on covariates.
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
