#' Transportation to campus
#'
#' @description
#' Data from a stratified sample of 60 respondents concerning their mode of
#' transportation to campus. The sampling design oversampled respondents who
#' sometimes traveled to campus by bicycle.
#'
#' @format A data frame with 60 observations and 7 variables:
#' \describe{
#'   \item{ntrips}{Number of trips to campus during the preceding four weeks.}
#'   \item{nbiked}{Number of those trips made by bicycle.}
#'   \item{status}{Respondent's institutional status: faculty, staff, or
#'     student.}
#'   \item{gender}{Respondent's gender, recorded as \code{"F"} or \code{"M"}.}
#'   \item{parking}{Duration of the parking permit, in months.}
#'   \item{distance}{Distance to campus. The source documentation does not
#'     specify the measurement unit.}
#'   \item{propbiked}{Proportion of trips to campus made by bicycle, calculated
#'     as \code{nbiked / ntrips}.}
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
#' data("transport", package = "vasicekreg")
#' str(transport)
#' with(transport, all.equal(propbiked, nbiked / ntrips))
#'
#' @name transport
#' @aliases transport
#' @docType data
#' @keywords datasets
NULL
