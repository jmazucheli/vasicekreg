# ======================================================================
# Zero-Augmented Beta model (formula interface)
# ======================================================================

#' Zero-Augmented Beta random-intercept model
#'
#' Fits a two-component model for responses in \eqn{[0,1)} with a point
#' mass at zero (discrete component, parameterized by \eqn{\gamma}) and a
#' Beta distribution for positive values (continuous component,
#' parameterized by \eqn{\beta}).
#'
#' @param data A data.frame.
#' @param y Character name of the response column.
#' @param formula_bin One-sided formula for the presence component.
#' @param formula_cont One-sided formula for the Beta component.
#' @param random Formula \code{~ 1 | Subject} for the random intercept.
#' @param logistic_cov Deprecated: character vector of presence covariates.
#' @param beta_cov Deprecated: character vector of Beta covariates.
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
#' \code{formula_cont} or \code{beta_cov}, and either \code{random}
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
#'   Y_{it} \sim \mathrm{Beta}\left(\mu_{it}\phi,\; (1-\mu_{it})\phi\right)
#'   \quad \mbox{with probability } p_{it},
#' }
#'
#' where \eqn{0 < p_{it} < 1}, \eqn{0 < \mu_{it} < 1}, and \eqn{\phi > 0}.
#' Let \eqn{X_{it}} and \eqn{Z_{it}} be the covariate vectors entering the
#' discrete and continuous components, respectively; they may share
#' columns or be disjoint. Both components are modeled on the logit
#' scale:
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
#' component) and \eqn{\beta} (continuous component); \eqn{\phi} is the
#' Beta precision parameter. The two random intercepts are independent, and the components have no shared parameters. Consequently,
#' the marginal likelihood factorizes into discrete and continuous components.
#' The two components are optimized separately; this is equivalent to
#' maximizing their joint likelihood. Random effects are integrated out using
#' non-adaptive Gauss--Hermite quadrature.
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
#' @return An object of class \code{"zabr"} with components including
#'   \code{logistic_est_table}, \code{beta_est_table},
#'   \code{precision_table}, \code{random_effects}, \code{loglikelihood},
#'   \code{joint_p}, \code{fit_statistics} (joint model),
#'   \code{fit_statistics_bin} and \code{fit_statistics_cont}
#'   (per component; see section \emph{Fit statistics}), and \code{vcov}.
#'   Fixed-effect names in \code{coef()} and \code{vcov()} use the prefix
#'   \code{gamma_} for the discrete component and \code{beta_} for the
#'   continuous component, including their intercepts.
#'
#' @seealso \code{\link{zavr}}
#'
#' @references
#' Chen, E. Z. and Li, H. (2016). A two-part mixed-effects model for
#' analyzing longitudinal microbiome compositional data.
#' \emph{Bioinformatics}, \bold{32}(17), 2611--2617.
#' \doi{10.1093/bioinformatics/btw308}
#'
#' @examples
#' \donttest{
#' data(please_microbiome)
#'
#' d <- subset(please_microbiome, Genus == "g__Bifidobacterium")
#' d$Subject <- factor(d$Subject)
#'
#' fit <- zabr(
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
zabr <- function(data, y, formula_bin = NULL, formula_cont = NULL, random = NULL,
                 logistic_cov = NULL, beta_cov = NULL, subject_ind = NULL, time_ind,
                 component_wise_test = TRUE, quad_n = 30, verbose = FALSE,
                 joint_test = NULL, sd_lower = 1e-5, start = NULL, control = list(), hessian = TRUE) {

  res <- .fit_zero_augmented_engine(
    data = data, y = y, formula_bin = formula_bin, formula_cont = formula_cont, random = random,
    logistic_cov = logistic_cov, positive_cov = beta_cov, subject_ind = subject_ind, time_ind = time_ind,
    component_wise_test = component_wise_test, quad_n = quad_n, verbose = verbose,
    joint_test = joint_test, sd_lower = sd_lower, start = start, control = control, hessian = hessian,
    required_pkgs = c("statmod", "numDeriv"),
    positive_name = "beta", shape_param_name = "phi",
    positive_density_fn = function(yy, mm, phi, sz) stats::dbeta(yy, shape1 = mm * phi, shape2 = (1 - mm) * phi, log = TRUE),
    validate_shape_fn = function(t2) { phi <- exp(t2); if (!is.finite(phi) || phi <= 0) NULL else phi },
    natural_shape_fn = function(val, validate_only = FALSE, get_jacobian = FALSE) {
      if (validate_only) { if (val <= 0) stop("Initial precision must be > 0."); return(log(val)) }
      nat <- exp(val); if (get_jacobian) list(val = nat, jac = nat) else nat
    },
    start_seeds_fn = function(log_sd_start, q, Y) lapply(c(2, 5, 10), function(phi0) c(log_sd_start, log(phi0), stats::qlogis(mean(Y[Y > 0])), rep(0, q - 1L)))
  )

  out <- list(
    call = match.call(),
    logistic_est_table = res$fixed_table(res$cv_l, 1L, res$logistic_names, res$stat_l),
    logistic_s1_est = unname(res$cv_l$estimate[1L]),
    beta_est_table = res$fixed_table(res$cv_v, 2L, res$positive_names, res$stat_v),
    beta_s2_est = unname(res$cv_v$estimate[1L]), beta_v_est = unname(res$cv_v$estimate[2L]),
    precision_table = data.frame(Estimate = unname(res$cv_v$estimate[2L]), SE = unname(res$cv_v$SE[2L]), row.names = "phi"),
    random_effects = res$random_table, loglikelihood = res$loglikelihood,
    joint_p = res$joint_p, joint_statistic = res$joint_statistic,
    fit_statistics_bin = res$stats_bin, fit_statistics_cont = res$stats_cont, fit_statistics = res$stats_total,
    estimates = res$estimates, vcov = res$vcov, diagnostics = res$diagnostics,
    quadrature = list(type = "nonadaptive Gauss-Hermite", points = quad_n),
    nobs = res$n, nsubjects = res$ns, sd_lower = sd_lower,
    versions = vapply(c("statmod", "numDeriv"), function(pkg) as.character(utils::packageVersion(pkg)), character(1)),
    optimization = list(logistic = res$full_l$opt, beta = res$full_v$opt),
    starts = list(logistic = res$full_l$candidates, beta = res$full_v$candidates)
  )
  class(out) <- "zabr"
  out
}

#' @export
print.zabr <- function(x, digits = 6, ...) {
  cat("Zero-augmented Beta random-intercept model\n")
  cat("Discrete component (gamma; Pvalue = LRT; Wald_Pvalue = normal reference):\n"); print(x$logistic_est_table, digits = digits)
  cat("\nContinuous component (beta):\n"); print(x$beta_est_table, digits = digits)
  cat("\nBeta precision:\n"); print(x$precision_table, digits = digits)
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
coef.zabr <- function(object, ...) object$estimates

#' @export
vcov.zabr <- function(object, ...) object$vcov

#' @export
logLik.zabr <- function(object, ...) {
  # nobs = number of subjects, so that stats::BIC() reproduces PROC NLMIXED
  # and BIC(fit1, fit2, ...) compares several models in one table.
  structure(object$loglikelihood,
            df    = length(object$estimates),
            nobs  = object$nsubjects,
            class = "logLik")
}

#' @export
nobs.zabr <- function(object, ...) object$nobs
