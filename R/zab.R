# ======================================================================
# Zero-Augmented Beta model (formula interface)
# ======================================================================

#' Zero-Augmented Beta random-intercept model
#'
#' Fits a two-component model for responses in [0,1) with a point mass at
#' zero (presence component) and a Beta distribution for positive values.
#'
#' @param data A data.frame.
#' @param y Character name of the response column.
#' @param formula_bin One-sided formula for the presence component.
#' @param formula_cont One-sided formula for the Beta component.
#' @param random Formula ~ 1 | Subject for the random intercept.
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
#' @return An object of class "zab".
#' @export
zab <- function(data, y, formula_bin = NULL, formula_cont = NULL, random = NULL,
                 logistic_cov = NULL, beta_cov = NULL, subject_ind = NULL, time_ind,
                 component_wise_test = TRUE, quad_n = 30, verbose = FALSE,
                 joint_test = NULL, sd_lower = 1e-5, start = NULL, control = list(), hessian = TRUE) {

  res <- .fit_zero_augmented_engine(
    data = data, y = y, formula_bin = formula_bin, formula_cont = formula_cont, random = random,
    logistic_cov = logistic_cov, positive_cov = beta_cov, subject_ind = subject_ind, time_ind = time_ind,
    component_wise_test = component_wise_test, quad_n = quad_n, verbose = verbose,
    joint_test = joint_test, sd_lower = sd_lower, start = start, control = control, hessian = hessian,
    model_name = "zab", required_pkgs = c("statmod", "numDeriv"),
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
  class(out) <- "zab"
  out
}

#' @export
print.zab <- function(x, digits = 6, ...) {
  cat("Zero-inflated Beta random-intercept model\n")
  cat("Presence component (Pvalue = LRT; Wald_Pvalue = normal reference):\n"); print(x$logistic_est_table, digits = digits)
  cat("\nPositive-abundance component:\n"); print(x$beta_est_table, digits = digits)
  cat("\nBeta precision:\n"); print(x$precision_table, digits = digits)
  cat("\nRandom effects:\n"); print(x$random_effects, digits = digits)
  if (!is.null(x$joint_p)) { cat("\nJoint likelihood-ratio p-values (2 df):\n"); print(x$joint_p, digits = digits) }
  cat("\nFit statistics (joint -2LL; BIC follows PROC NLMIXED: ",
      "BIC = 2f + p*log(subjects)):\n", sep = "")
  print(x$fit_statistics, digits = digits, row.names = FALSE)
  cat("\nComponent-wise fit statistics (informational, NOT additive for BIC):\n")
  print(rbind(cbind(Component = "presence", x$fit_statistics_bin),
              cbind(Component = "positive", x$fit_statistics_cont)),
        digits = digits, row.names = FALSE)
  cat("\nFull-model diagnostics:\n"); print(x$diagnostics[1:2, ], row.names = FALSE)
  invisible(x)
}

#' @export
coef.zab <- function(object, ...) object$estimates

#' @export
vcov.zab <- function(object, ...) object$vcov

#' @export
logLik.zab <- function(object, ...) {
  structure(object$loglikelihood,
            df    = length(object$estimates),
            nobs  = object$nobs,
            class = "logLik")
}

#' @export
nobs.zab <- function(object, ...) object$nobs

#' @export
# BIC compatible with SAS/NLMIXED: uses log(number of subjects)
BIC.zab <- function(object, ...) {
  k <- length(object$estimates)
  -2 * as.numeric(stats::logLik(object)) + k * log(object$nsubjects)
}
