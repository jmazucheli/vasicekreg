# ======================================================================
# Zero-Augmented Vasicek model (formula interface)
# ======================================================================

#' Zero-Augmented Vasicek random-intercept model
#'
#' Fits a two-component model for responses in \eqn{[0,1)} with a point
#' mass at zero (discrete component, parameterized by \eqn{\gamma}) and a
#' Vasicek distribution for positive values (continuous component,
#' parameterized by \eqn{\beta}).
#'
#' @param data A data.frame.
#' @param y Character name of the response column.
#' @param formula_bin One-sided formula for the presence component.
#' @param formula_cont One-sided formula for the Vasicek component.
#' @param random Formula \code{~ 1 | Subject} for the random intercept.
#' @param logistic_cov Deprecated: character vector of presence covariates.
#' @param vasicek_cov Deprecated: character vector of Vasicek covariates.
#' @param subject_ind Deprecated: subject column name.
#' @param time_ind Character name of the time column.
#' @param component_wise_test Logical; compute component-wise LRTs?
#' @param quad_n Integer; number of Gauss-Hermite quadrature points.
#' @param verbose Logical; print progress?
#' @param joint_test NULL, TRUE, or FALSE; compute joint LRTs?
#' @param sd_lower Lower bound for the random-effect SD.
#' @param start Optional list of starting values.
#' @param control List of control parameters.
#' @param hessian Logical; compute Hessian-based standard errors?
#'
#' @details
#' Formula and legacy arguments are mutually exclusive for each model
#' component: use either \code{formula_bin} or \code{logistic_cov}, either
#' \code{formula_cont} or \code{vasicek_cov}, and either \code{random}
#' or \code{subject_ind}. Supplying both arguments in any pair is an error.
#'
#' For a response \eqn{Y_{it} \in [0,1)} observed on subject \eqn{i}
#' (\eqn{i = 1, \ldots, N}) at time \eqn{t} (\eqn{t = 1, \ldots, T}),
#' the model places a point mass at zero and a continuous component on
#' \eqn{(0,1)}:
#'
#' \deqn{
#'   Y_{it} = 0 \quad \mbox{with probability } 1 - p_{it},
#' }
#' \deqn{
#'   Y_{it} \sim \mathrm{NVASIM}\left(\mu_{it}, \sigma\right)
#'   \quad \mbox{with probability } p_{it},
#' }
#'
#' where \eqn{0 < p_{it} < 1}, \eqn{0 < \mu_{it} < 1}, and
#' \eqn{\sigma \in (0,1)}. The continuous component is the normal-kernel
#' Vasicek (mean-parameterized) distribution \code{NVASIM}, for which
#' \eqn{\mu_{it} = E(Y_{it} \mid Y_{it} > 0)}. Let \eqn{X_{it}} and
#' \eqn{Z_{it}} be the covariate vectors entering the discrete and
#' continuous components, respectively; they may share columns or be
#' disjoint. Both components are modeled on the logit scale:
#'
#' \deqn{
#'   \mathrm{logit}(p_{it}) = \log\left(\frac{p_{it}}{1 - p_{it}}\right)
#'   = a_i + \gamma_0 + X_{it}^\top \gamma,
#' }
#' \deqn{
#'   \mathrm{logit}(\mu_{it}) = \log\left(\frac{\mu_{it}}{1 - \mu_{it}}\right)
#'   = b_i + \beta_0 + Z_{it}^\top \beta,
#' }
#'
#' where \eqn{a_i} and \eqn{b_i} are subject-specific random intercepts
#' that induce correlation across repeated measurements on the same
#' subject,
#'
#' \deqn{
#'   a_i \sim N(0, \sigma_1^2), \qquad b_i \sim N(0, \sigma_2^2).
#' }
#'
#' The fixed-effect coefficients are denoted by \eqn{\gamma} (discrete
#' component) and \eqn{\beta} (continuous component); \eqn{\sigma} is the
#' shape parameter of the Vasicek component. The two random intercepts are
#' independent, and the components have no shared parameters. Consequently,
#' the marginal likelihood factorizes into discrete and continuous components.
#' The two components are optimized separately; this is equivalent to
#' maximizing their joint likelihood. Random effects are integrated out using
#' non-adaptive Gauss--Hermite quadrature.
#'
#' This function differs from \code{\link{zabr}} only in the continuous
#' component: \code{zabr} uses a Beta distribution with precision
#' \eqn{\phi}, whereas \code{zavr} uses the Vasicek \code{NVASIM}
#' distribution with shape \eqn{\sigma}. The two models share the same
#' discrete component and the same random-intercept structure.
#'
#' @section Fit statistics:
#' \code{fit_statistics} refers to the joint model and follows the
#' "Fit Statistics" table of SAS PROC NLMIXED:
#' \eqn{AIC = -2\ell + 2k}{AIC = -2 logL + 2k} and
#' \eqn{BIC = -2\ell + k\log(s)}{BIC = -2 logL + k log(s)}, where \eqn{\ell}{logL} is the maximized
#' log-likelihood, \eqn{k} the total number of parameters and \eqn{s} the
#' number of subjects. BIC uses subjects, not observations. \code{AIC()} and
#' \code{BIC()} return the same values as \code{fit_statistics}.
#'
#' \code{fit_statistics_bin} and \code{fit_statistics_cont} are informational.
#' The continuous component counts only observations with \eqn{Y > 0} and
#' subjects with at least one such observation, as if it were fitted alone to
#' the positive responses. Consequently, AIC is additive across components,
#' whereas BIC is additive only when every subject has at least one positive
#' response. To compare models, use \code{fit_statistics} or \code{BIC()};
#' do not sum component-wise values.
#'
#' \code{logLik()} sets the \code{nobs} attribute to the number of subjects,
#' so that \code{BIC()} follows PROC NLMIXED.
#'
#' The values coincide numerically with PROC NLMIXED only when the
#' log-likelihood coincides: the same model with independent random
#' intercepts and the same number of parameters, non-adaptive quadrature with
#' the same number of points (\code{NOAD} and \code{QPOINTS=} equal to
#' \code{quad_n}), and data sorted by subject.
#'
#' @return An object of class \code{"zavr"} with components including
#'   \code{logistic_est_table}, \code{vasicek_est_table},
#'   \code{shape_table}, \code{random_effects}, \code{loglikelihood},
#'   \code{joint_p}, \code{fit_statistics} (joint model),
#'   \code{fit_statistics_bin} and \code{fit_statistics_cont}
#'   (per component; see section \emph{Fit statistics}), and \code{vcov}.
#'   Fixed-effect names in \code{coef()} and \code{vcov()} use the prefix
#'   \code{gamma_} for the discrete component and \code{beta_} for the
#'   continuous component, including their intercepts.
#'
#' @seealso \code{\link{zabr}}, \code{\link{NVASIM}}
#'
#' @references
#' Chen, E. Z. and Li, H. (2016). A two-part mixed-effects model for
#' analyzing longitudinal microbiome compositional data.
#' \emph{Bioinformatics}, \bold{32}(17), 2611--2617.
#' \doi{10.1093/bioinformatics/btw308}
#'
#' Mazucheli, J., Alves, B., Korkmaz, M. Ç., and Leiva, V. (2022).
#' Vasicek quantile and mean regression models for bounded data: New
#' formulation, mathematical derivations, and numerical applications.
#' \emph{Mathematics}, \bold{10}, 1389.
#' \doi{10.3390/math10091389}
#'
#' @examples
#' \donttest{
#' data(please_microbiome)
#'
#' d <- subset(please_microbiome, Genus == "g__Bifidobacterium")
#' d$Subject <- factor(d$Subject)
#'
#' fit <- zavr(
#'   data         = d,
#'   y            = "Y",
#'   formula_bin  = ~ Baseline + Time + Treat,
#'   formula_cont = ~ Baseline + Time + Treat,
#'   random       = ~ 1 | Subject,
#'   time_ind     = "Time"
#' )
#' fit
#' }
#'
#' @export
zavr <- function(data, y, formula_bin = NULL, formula_cont = NULL, random = NULL,
                 logistic_cov = NULL, vasicek_cov = NULL, subject_ind = NULL, time_ind,
                 component_wise_test = TRUE, quad_n = 30, verbose = FALSE,
                 joint_test = NULL, sd_lower = 1e-5, start = NULL, control = list(), hessian = TRUE) {

  res <- .fit_zero_augmented_engine(
    data = data, y = y, formula_bin = formula_bin, formula_cont = formula_cont, random = random,
    logistic_cov = logistic_cov, positive_cov = vasicek_cov, subject_ind = subject_ind, time_ind = time_ind,
    component_wise_test = component_wise_test, quad_n = quad_n, verbose = verbose,
    joint_test = joint_test, sd_lower = sd_lower, start = start, control = control, hessian = hessian,
    required_pkgs = c("statmod", "numDeriv"),
    positive_name = "vasicek", shape_param_name = "shape",
    positive_density_fn = function(yy, mm, shape, sz) dNVASIM(x = yy, mu = mm, sigma = rep(shape, sz), log = TRUE),
    validate_shape_fn = function(t2) { shape <- stats::plogis(t2); if (!is.finite(shape) || shape <= 0 || shape >= 1) NULL else shape },
    natural_shape_fn = function(val, validate_only = FALSE, get_jacobian = FALSE) {
      if (validate_only) { if (val <= 0 || val >= 1) stop("Initial shape must be in (0,1)."); return(stats::qlogis(val)) }
      nat <- stats::plogis(val); if (get_jacobian) list(val = nat, jac = nat * (1 - nat)) else nat
    },
    start_seeds_fn = function(log_sd_start, q, Y) lapply(c(0.25, 0.5, 0.75), function(s) c(log_sd_start, stats::qlogis(s), stats::qlogis(mean(Y[Y > 0])), rep(0, q - 1L)))
  )

  out <- list(
    call = match.call(),
    logistic_est_table = res$fixed_table(res$cv_l, 1L, res$logistic_names, res$stat_l),
    logistic_s1_est = unname(res$cv_l$estimate[1L]),
    vasicek_est_table = res$fixed_table(res$cv_v, 2L, res$positive_names, res$stat_v),
    vasicek_s2_est = unname(res$cv_v$estimate[1L]), vasicek_shape_est = unname(res$cv_v$estimate[2L]),
    shape_table = data.frame(Estimate = unname(res$cv_v$estimate[2L]), SE = unname(res$cv_v$SE[2L]), row.names = "shape"),
    random_effects = res$random_table, loglikelihood = res$loglikelihood,
    joint_p = res$joint_p, joint_statistic = res$joint_statistic,
    fit_statistics_bin = res$stats_bin, fit_statistics_cont = res$stats_cont, fit_statistics = res$stats_total,
    estimates = res$estimates, vcov = res$vcov, diagnostics = res$diagnostics,
    quadrature = list(type = "nonadaptive Gauss-Hermite", points = quad_n),
    nobs = res$n, nsubjects = res$ns, sd_lower = sd_lower,
    versions = vapply(c("statmod", "numDeriv"), function(pkg) as.character(utils::packageVersion(pkg)), character(1)),
    optimization = list(logistic = res$full_l$opt, vasicek = res$full_v$opt),
    starts = list(logistic = res$full_l$candidates, vasicek = res$full_v$candidates)
  )
  class(out) <- "zavr"
  out
}

#' @export
print.zavr <- function(x, digits = 6, ...) {
  cat("Zero-augmented Vasicek random-intercept model\n")
  cat("Discrete component (gamma; Pvalue = LRT; Wald_Pvalue = normal reference):\n"); print(x$logistic_est_table, digits = digits)
  cat("\nContinuous component (beta):\n"); print(x$vasicek_est_table, digits = digits)
  cat("\nVasicek shape:\n"); print(x$shape_table, digits = digits)
  cat("\nRandom effects:\n"); print(x$random_effects, digits = digits)
  if (!is.null(x$joint_p)) { cat("\nJoint likelihood-ratio p-values (2 df):\n"); print(x$joint_p, digits = digits) }
  cat("\nFit statistics, joint model (BIC uses subjects, as in PROC NLMIXED):\n")
  print(x$fit_statistics, digits = digits, row.names = FALSE)
  cat("\nComponent-wise fit statistics (informational; AIC is additive,\n",
      "BIC is additive only if every subject has at least one Y > 0):\n", sep = "")
  print(rbind(cbind(Component = "presence", x$fit_statistics_bin),
              cbind(Component = "positive", x$fit_statistics_cont)),
        digits = digits, row.names = FALSE)
  cat("\nFull-model diagnostics:\n"); print(x$diagnostics[1:2, ], row.names = FALSE)
  invisible(x)
}

#' @export
coef.zavr <- function(object, ...) object$estimates

#' @export
vcov.zavr <- function(object, ...) object$vcov

#' @export
logLik.zavr <- function(object, ...) {
  # nobs = number of subjects, so that stats::BIC() reproduces PROC NLMIXED
  # and BIC(fit1, fit2, ...) compares several models in one table.
  structure(object$loglikelihood,
            df    = length(object$estimates),
            nobs  = object$nsubjects,
            class = "logLik")
}

#' @export
nobs.zavr <- function(object, ...) object$nobs
