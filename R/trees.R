#' Survival of young trees
#'
#' @description
#' Data from a study conducted by a parks and recreation department on the
#' two-year survival of young trees planted in 26 parks.
#'
#' @format A data frame with 26 observations and 7 variables:
#' \describe{
#'   \item{planted}{Number of trees planted.}
#'   \item{survived}{Number of trees that survived for two years.}
#'   \item{pest}{Frequency of pest-control treatment.}
#'   \item{fertilization}{Frequency of soil fertilization.}
#'   \item{precip}{Average annual precipitation, in inches.}
#'   \item{wind}{Average annual wind speed, in miles per hour.}
#'   \item{prop}{Proportion of planted trees that survived for two years,
#'     calculated as \code{survived / planted}.}
#' }
#'
#' @details
#' The data originate from a consulting study reported by Korosteleva (2019).
#' Menezes, Mazucheli, and Bourguignon (2021) analyzed them using an inflated
#' unit-Weibull quantile regression model. The object distributed here retains
#' the variable names and values supplied by the \pkg{uwquantreg} package.
#'
#' @source
#' Korosteleva, O. (2019). \emph{Advanced Regression Models with SAS and R}.
#' Boca Raton, FL: CRC Press.
#'
#' \url{https://github.com/AndrMenezes/uwquantreg}
#'
#' @references
#' Menezes, A. F. B., Mazucheli, J. and Bourguignon, M. (2021). A parametric
#' quantile regression approach for modeling zero- or one-inflated double
#' bounded data. \emph{Biometrical Journal}, \bold{63}(4), 841--858.
#' \doi{10.1002/bimj.202000126}
#'
#' Menezes, A. F. B. (2026). \pkg{uwquantreg}: unit-Weibull quantile
#' regression. R package version 0.1.0.
#' \url{https://github.com/AndrMenezes/uwquantreg}
#'
#' @examples
#' data("trees", package = "vasicekreg")
#' str(trees)
#' with(trees, all.equal(prop, survived / planted))
#'
#' @name trees
#' @aliases trees
#' @docType data
#' @keywords datasets
NULL
