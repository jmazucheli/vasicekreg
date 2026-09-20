# ======================================================================
# Zero-Augmented Vasicek model (formula interface)
# ======================================================================

#' Zero-Augmented Vasicek random-intercept model
#'
#' Fits a two-component model for responses in [0,1) with a point mass at
#' zero (presence component) and a Vasicek distribution for positive values.
#'
#' @param data A data.frame.
#' @param y Character name of the response column.
#' @param formula_bin One-sided formula for the presence component.
#' @param formula_cont One-sided formula for the Vasicek component.
#' @param random Formula ~ 1 | Subject for the random intercept.
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
#' @return An object of class "zavr".
#' @export
zavr <- function(data, y, formula_bin = NULL, formula_cont = NULL, random = NULL,
                 logistic_cov = NULL, vasicek_cov = NULL, subject_ind = NULL, time_ind,
                 component_wise_test = TRUE, quad_n = 30, verbose = FALSE,
                 joint_test = NULL, sd_lower = 1e-5, start = NULL, control = list(), hessian = TRUE) {

  if (!requireNamespace("vasicekreg", quietly = TRUE))
    stop("Please install package 'vasicekreg'.")
  if (!"dNVASIM" %in% getNamespaceExports("vasicekreg"))
    stop("The vasicekreg package must export dNVASIM().")

  res <- .fit_zero_augmented_engine(
    data = data, y = y, formula_bin = formula_bin, formula_cont = formula_cont, random = random,
    logistic_cov = logistic_cov, positive_cov = vasicek_cov, subject_ind = subject_ind, time_ind = time_ind,
    component_wise_test = component_wise_test, quad_n = quad_n, verbose = verbose,
    joint_test = joint_test, sd_lower = sd_lower, start = start, control = control, hessian = hessian,
    model_name = "zavr", required_pkgs = c("vasicekreg", "statmod", "numDeriv"),
    positive_name = "vasicek", shape_param_name = "shape",
    positive_density_fn = function(yy, mm, shape, sz) vasicekreg::dNVASIM(x = yy, mu = mm, sigma = rep(shape, sz), log = TRUE),
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
    versions = vapply(c("vasicekreg", "statmod", "numDeriv"), function(pkg) as.character(utils::packageVersion(pkg)), character(1)),
    optimization = list(logistic = res$full_l$opt, vasicek = res$full_v$opt),
    starts = list(logistic = res$full_l$candidates, vasicek = res$full_v$candidates)
  )
  class(out) <- "zavr"
  out
}

#' @export
print.zavr <- function(x, digits = 6, ...) {
  cat("Zero-augmented Vasicek random-intercept model\n")
  cat("Presence component (Pvalue = LRT; Wald_Pvalue = normal reference):\n"); print(x$logistic_est_table, digits = digits)
  cat("\nPositive-abundance component:\n"); print(x$vasicek_est_table, digits = digits)
  cat("\nVasicek shape:\n"); print(x$shape_table, digits = digits)
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
coef.zavr <- function(object, ...) object$estimates

#' @export
vcov.zavr <- function(object, ...) object$vcov

#' @export
logLik.zavr <- function(object, ...) {
  structure(object$loglikelihood,
            df    = length(object$estimates),
            nobs  = object$nobs,
            class = "logLik")
}

#' @export
nobs.zavr <- function(object, ...) object$nobs

#' @export
# BIC compatible with SAS/NLMIXED: uses log(number of subjects)
BIC.zavr <- function(object, ...) {
  k <- length(object$estimates)
  -2 * as.numeric(stats::logLik(object)) + k * log(object$nsubjects)
}
