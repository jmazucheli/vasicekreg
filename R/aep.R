#' Hospital-stay appropriateness data
#'
#' @description
#' Data on 1,383 patients admitted to Hospital del Mar, Barcelona, during 1988
#' and 1990. The data record the total length of stay and the number of days
#' classified as inappropriate.
#'
#' @format A data frame with 1,383 observations and 8 variables:
#' \describe{
#'   \item{los}{Total number of days spent in hospital.}
#'   \item{noinap}{Number of hospital-stay days classified as inappropriate.}
#'   \item{loglos}{Logarithm of length of stay divided by 10,
#'     \code{log(los / 10)}.}
#'   \item{sex}{Patient's sex: factor with levels \code{1} (male) and
#'     \code{2} (female).}
#'   \item{ward}{Hospital ward: factor with levels \code{1} (medical),
#'     \code{2} (surgical), and \code{3} (other).}
#'   \item{year}{Admission year: factor with levels \code{88} and \code{90}.}
#'   \item{age}{Patient's age minus 55 years.}
#'   \item{y}{Two-column matrix response whose first column is \code{noinap}
#'     and whose second column is \code{los - noinap}.}
#' }
#'
#' @details
#' Gange et al. (1996) modeled the number of inappropriate days conditional on
#' total length of stay using binomial and beta-binomial models. For augmented
#' continuous-response models, the patient-level proportion can be constructed
#' as \code{noinap / los}. The object distributed here retains the structure and
#' values supplied by the \pkg{gamlss.data} package.
#'
#' @source
#' Gange, S. J., Munoz, A., Saez, M. and Alonso, J. (1996). Use of the
#' beta-binomial distribution to model the effect of policy changes on
#' appropriateness of hospital stays. \emph{Applied Statistics}, \bold{45}(3),
#' 371--382.
#'
#' @references
#' Stasinopoulos, M. and Rigby, R. (2025). \pkg{gamlss.data}: Data for
#' Generalized Additive Models for Location Scale and Shape. R package version
#' 6.0-7. \url{https://CRAN.R-project.org/package=gamlss.data}
#'
#' @examples
#' data("aep", package = "vasicekreg")
#' str(aep)
#' inappropriate <- with(aep, noinap / los)
#' table(inappropriate == 0, inappropriate == 1)
#'
#' @name aep
#' @aliases aep
#' @docType data
#' @keywords datasets
NULL
