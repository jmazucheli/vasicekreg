#' @title Overview of the vasicekreg package
#'
#' @description
#' The \pkg{vasicekreg} package provides distribution functions and GAMLSS
#' families for Vasicek-type distributions on the unit interval. Four base
#' families are available:
#' \itemize{
#'   \item \code{NVASIM}: normal kernel with mean parameterization, where
#'   \eqn{\mu=E(Y)}.
#'   \item \code{NVASIQ}: normal kernel with quantile parameterization, where
#'   \eqn{\mu=Q_Y(\tau)} for a fixed \eqn{\tau\in(0,1)}.
#'   \item \code{LVASIQ}: logistic kernel with quantile parameterization, where
#'   \eqn{\mu=Q_Y(\tau)} for a fixed \eqn{\tau\in(0,1)}.
#'   \item \code{HVASIQ}: hyperbolic-secant kernel with quantile
#'   parameterization, where \eqn{\mu=Q_Y(\tau)} for a fixed
#'   \eqn{\tau\in(0,1)}.
#' }
#' For responses observed at the boundaries, the normal-kernel mean model is
#' also available as:
#' \itemize{
#'   \item \code{ZANVASIM}: point mass at zero and a continuous component on
#'   \eqn{(0,1)}.
#'   \item \code{OANVASIM}: point mass at one and a continuous component on
#'   \eqn{(0,1)}.
#'   \item \code{ZOANVASIM}: point masses at zero and one and a continuous
#'   component on \eqn{(0,1)}.
#' }
#' The documentation uses \emph{augmented} for these boundary mixtures and
#' \emph{Vasicek-type} for kernel-based constructions. The established family
#' names are retained for backward compatibility and follow familiar GAMLSS
#' abbreviations in which \code{ZA} and \code{OA} historically denote
#' zero- and one-adjusted families.
#' The shape parameter \eqn{\sigma\in(0,1)} controls dispersion in the continuous
#' Vasicek component. The corresponding \code{d}, \code{p}, \code{q}, and
#' \code{r} functions provide density or probability mass values, cumulative
#' probabilities, quantiles, and random observations, respectively.
#'
#' @details
#' Included datasets:
#' \itemize{
#'   \item \code{\link[vasicekreg]{bodyfat}}: body-fat proportions in
#'   \eqn{(0,1)} and demographic covariates for 298 individuals.
#'   \item \code{\link[vasicekreg]{aep}}: hospital-stay data from 1,383
#'   patients; \code{noinap / los} contains observations at zero and one.
#'   \item \code{\link[vasicekreg]{transport}}: bicycle-trip proportions for
#'   60 respondents, including observations at zero.
#'   \item \code{\link[vasicekreg]{trees}}: two-year tree-survival
#'   proportions for 26 parks, including observations at one.
#' }
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
#' \code{\link[vasicekreg]{HVASIQ}}:
#' Hyperbolic-secant-kernel quantile parameterization and GAMLSS family. For a
#' fixed quantile level \eqn{\tau}, covariates describe the conditional
#' \eqn{\tau}-th quantile through \eqn{\mu}. Its conditional mean and
#' variance are obtained by numerical quadrature and \eqn{\mu} must not be
#' interpreted as the mean.
#'
#' \code{\link[vasicekreg]{ZANVASIM}}:
#' Zero-augmented normal-kernel mean family. Here
#' \eqn{\nu=P(Y=0)}, \eqn{\mu=E(Y\mid Y>0)}, and the marginal mean is
#' \eqn{E(Y)=(1-\nu)\mu}.
#'
#' \code{\link[vasicekreg]{OANVASIM}}:
#' One-augmented normal-kernel mean family. Here
#' \eqn{\nu=P(Y=1)}, \eqn{\mu=E(Y\mid Y<1)}, and the marginal mean is
#' \eqn{E(Y)=\nu+(1-\nu)\mu}. The parameters \eqn{\mu} and \eqn{\nu}
#' therefore have the same interpretations as their counterparts in the
#' one-inflated beta family \code{\link[gamlss.dist]{BEOI}}. The shape
#' parameter \eqn{\sigma} is distribution-specific and should not be
#' compared directly between these families.
#'
#' \code{\link[vasicekreg]{ZOANVASIM}}:
#' Zero-and-one-augmented normal-kernel mean family. Here
#' \eqn{\nu=P(Y=0)}, \eqn{\tau=P(Y=1\mid Y>0)}, and
#' \eqn{\mu=E(Y\mid 0<Y<1)}. Consequently,
#' \eqn{P(Y=1)=(1-\nu)\tau} and
#' \eqn{E(Y)=(1-\nu)[\tau+(1-\tau)\mu]}.
#'
#' The distribution functions \code{dNVASIM}, \code{pNVASIM},
#' \code{qNVASIM}, \code{dNVASIQ}, \code{pNVASIQ}, \code{qNVASIQ},
#' \code{dLVASIQ}, \code{pLVASIQ}, \code{qLVASIQ},
#' \code{dHVASIQ}, \code{pHVASIQ}, and \code{qHVASIQ} call compiled
#' \proglang{C++} routines through \pkg{Rcpp}. The boundary-augmented
#' distribution functions are implemented in \proglang{R} and reuse the
#' compiled \code{NVASIM} functions for their continuous component.
#' Parameter validation, the GAMLSS family definitions, and all
#' log-likelihood derivatives are implemented in \proglang{R}. The mean and
#' variance components of the \code{LVASIQ()} and \code{HVASIQ()} family
#' objects are obtained by numerical quadrature because these moments have no
#' closed-form expressions.
#'
#' For the distribution functions and GAMLSS constructors associated with
#' \code{NVASIQ}, \code{LVASIQ}, and \code{HVASIQ}, the fixed quantile level
#' is supplied through the \code{quantile} argument. It is stored in the
#' family definition and embedded in the residual expression, so no global
#' variable is required. This fixed quantile level is distinct from the
#' parameter \code{tau} in \code{ZOANVASIM()}, which represents the
#' conditional probability at one among nonzero observations.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("betareg", quietly = TRUE)) {
#'     data("ReadingSkills", package = "betareg")
#'     ReadingSkills$dyslexia <- stats::relevel(
#'         factor(ReadingSkills$dyslexia), ref = "no"
#'     )
#'
#'     control <- gamlss::gamlss.control(n.cyc = 200, trace = FALSE)
#'
#'     ## In both models, mu = E(Y | Y < 1) and nu = P(Y = 1).
#'     fit_oanvasim <- gamlss::gamlss(
#'         accuracy1 ~ dyslexia * iq,
#'         sigma.formula = ~ dyslexia + iq,
#'         nu.formula = ~ 1,
#'         family = OANVASIM(),
#'         data = ReadingSkills,
#'         control = control
#'     )
#'
#'     fit_beoi <- gamlss::gamlss(
#'         accuracy1 ~ dyslexia * iq,
#'         sigma.formula = ~ dyslexia + iq,
#'         nu.formula = ~ 1,
#'         family = gamlss.dist::BEOI(),
#'         data = ReadingSkills,
#'         control = control
#'     )
#'
#'     n <- nrow(ReadingSkills)
#'     comparison <- data.frame(
#'         family = c("OANVASIM", "BEOI"),
#'         logLik = -c(
#'             fit_oanvasim$G.deviance,
#'             fit_beoi$G.deviance
#'         ) / 2,
#'         AIC = c(
#'             gamlss::GAIC(fit_oanvasim, k = 2),
#'             gamlss::GAIC(fit_beoi, k = 2)
#'         ),
#'         BIC = c(
#'             gamlss::GAIC(fit_oanvasim, k = log(n)),
#'             gamlss::GAIC(fit_beoi, k = log(n))
#'         )
#'     )
#'     comparison
#' }
#' }
#'
#' @author
#' Josmar Mazucheli \email{jmazucheli@gmail.com}
#'
#' Bruna Alves \email{pg402900@uem.br}
#'
#' @references
#' Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals.
#' \emph{Journal of Computational and Graphical Statistics}, \bold{5}(3),
#' 236--244. \doi{10.1080/10618600.1996.10474708}
#'
#' Fischer, M. J., Hui, A. and Hösle, S. (2017). wHS-type distributions
#' with application to finance. \emph{Journal of Statistics and Management
#' Systems}, \bold{20}(1), 67--89.
#' \doi{10.1080/09720510.2016.1190575}
#'
#' Mazucheli, J., Alves, B., Korkmaz, M. Ç., and Leiva, V. (2022).
#' Vasicek quantile and mean regression models for bounded data: New
#' formulation, mathematical derivations, and numerical applications.
#' \emph{Mathematics}, \bold{10}, 1389. \doi{10.3390/math10091389}
#'
#' Moral, R. A., Hinde, J. and Demetrio, C. G. B. (2017). Half-normal plots
#' and overdispersed models in R: The hnp package. \emph{Journal of
#' Statistical Software}, \bold{81}(10), 1--23.
#' \doi{10.18637/jss.v081.i10}
#'
#' Ospina, R. and Ferrari, S. L. P. (2010). Inflated beta distributions.
#' \emph{Statistical Papers}, \bold{51}, 111--126.
#' \doi{10.1007/s00362-008-0125-4}
#'
#' Ospina, R. and Ferrari, S. L. P. (2012). A general class of zero-or-one
#' inflated beta regression models. \emph{Computational Statistics & Data
#' Analysis}, \bold{56}(6), 1609--1623.
#' \doi{10.1016/j.csda.2011.10.005}
#'
#' Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive models
#' for location, scale and shape. \emph{Applied Statistics}, \bold{54}(3),
#' 507--554. \doi{10.1111/j.1467-9876.2005.00510.x}
#'
#' Vasicek, O. A. (2002). The distribution of loan portfolio value.
#' \emph{Risk}, \bold{15}(12), 160--162.
#'
#' Witzany, J. (2013). A note on the Vasicek's model with the logistic
#' distribution. \emph{Ekonomický časopis (Journal of Economics)},
#' \bold{61}(10), 1053--1066.
#'
#' Zhao, Y., Lee, A. H., Yau, K. K. W. and McLachlan, G. J. (2011).
#' Assessing the adequacy of Weibull survival models: A simulated envelope
#' approach. \emph{Journal of Applied Statistics}, \bold{38}(10),
#' 2089--2097. \doi{10.1080/02664763.2010.545115}
#'
#' @useDynLib vasicekreg
#' @importFrom Rcpp evalCpp
"_PACKAGE"

.onUnload <- function(libpath) {
    library.dynam.unload("vasicekreg", libpath)
}
