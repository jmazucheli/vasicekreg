#' @name bodyfat
#' @aliases bodyfat
#'
#' @title Percentage of Body Fat Dataset
#'
#' @description Percentage of body fat measurements from individuals
#' assisted in a public hospital in Curitiba, Paraná, Brazil.
#'
#' @format A data frame with 298 observations and 9 variables:
#'
#' \itemize{
#'   \item \code{ARMS}: arms fat percentage.
#'   \item \code{LEGS}: legs fat percentage.
#'   \item \code{BODY}: body fat percentage.
#'   \item \code{ANDROID}: android fat percentage.
#'   \item \code{GYNECOID}: gynoid fat percentage.
#'   \item \code{AGE}: age of individuals.
#'   \item \code{BMI}: body mass index.
#'   \item \code{SEX}: 1 for female and 2 for male.
#'   \item \code{IPAQ}: physical activity level according to IPAQ
#'   (0 = sedentary, 1 = insufficiently active, 2 = active).
#' }
#'
#' @author
#' Josmar Mazucheli \email{jmazucheli@gmail.com}
#'
#' Bruna Alves \email{pg402900@uem.br}
#'
#' @references
#' Mazucheli, J., Alves, B., Korkmaz, M. Ç., and Leiva, V. (2022). Vasicek quantile and mean regression models for bounded data: New formulation, mathematical derivations, and numerical applications. \emph{Mathematics}, \bold{10}, 1389.
#'
#' Mazucheli, J., Leiva, V., Alves, B., and Menezes, A. F. B. (2021). A new quantile regression for modeling bounded data under a unit Birnbaum-Saunders distribution with applications in medicine and politics. \emph{Symmetry}, \bold{13}(4), 1--21.
#'
#' Petterle, R. R., Bonat, W. H., Scarpin, C. T., Jonasson, T., and Borba, V. Z. C. (2020). Multivariate quasi-beta regression models for continuous bounded data. \emph{The International Journal of Biostatistics}, \bold{17}(1), 39--53.
#'
#' @examples
#' data(bodyfat, package = "vasicekreg")
#'
#' bodyfat$AGE <- bodyfat$AGE - 46.00
#' bodyfat$BMI <- bodyfat$BMI - 24.72
#' bodyfat$SEX <- as.factor(bodyfat$SEX)
#' bodyfat$IPAQ<- as.factor(bodyfat$IPAQ)
#'
#' library(gamlss)
#'
#' ## Mean regression model
#' fitmean <- gamlss(
#'   ARMS ~ AGE + BMI + SEX + IPAQ,
#'   data = bodyfat,
#'   family = NVASIM(mu.link = "logit", sigma.link = "logit")
#' )
#'
#' \dontrun{
#' tau_levels <- c(0.25, 0.50, 0.75)
#'
#' ## Quantile regression models with the normal kernel
#' fit_normal <- lapply(tau_levels, function(Tau) {
#'   tau <<- Tau
#'   gamlss(
#'     ARMS ~ AGE + BMI + SEX + IPAQ,
#'     data = bodyfat,
#'     family = NVASIQ(
#'       mu.link = "logit",
#'       sigma.link = "logit"
#'     )
#'   )
#' })
#'
#' ## Quantile regression models with the logistic kernel
#' fit_logistic <- lapply(tau_levels, function(Tau) {
#'   tau <<- Tau
#'   gamlss(
#'     ARMS ~ AGE + BMI + SEX + IPAQ,
#'     data = bodyfat,
#'     family = LVASIQ(
#'       mu.link = "logit",
#'       sigma.link = "logit"
#'     )
#'   )
#' })
#' 
#' ## Quantile regression models with the Hyperbolic-secant-kernel
#' fit_hsk <- lapply(tau_levels, function(Tau) {
#'   tau <<- Tau
#'   gamlss(
#'     ARMS ~ AGE + BMI + SEX + IPAQ,
#'     data = bodyfat,
#'     family = HSVASIQ(
#'       mu.link = "logit",
#'       sigma.link = "logit"
#'     )
#'   )
#' })
#'
#' lapply(fit_normal, summary)
#' lapply(fit_logistic, summary)
#' lapply(fit_hsk, summary)
#' }
"bodyfat"
