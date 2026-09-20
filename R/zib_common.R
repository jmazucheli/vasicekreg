# ======================================================================
# Shared computation engine (formula interface + SAS-aligned statistics)
# ======================================================================

#' Build a design matrix from a formula or character vector
#'
#' @param formula_or_names A one-sided formula (e.g., ~ x1 + x2) or a
#'   character vector of column names.
#' @param data A data.frame containing the variables.
#' @param label Character label for error messages.
#' @param legacy_arg_name Optional name of a deprecated argument.
#' @return A numeric matrix without the intercept column.
#' @keywords internal
#' @noRd
.build_design_matrix <- function(formula_or_names, data, label, legacy_arg_name = NULL) {
  if (is.null(formula_or_names) || (is.character(formula_or_names) && length(formula_or_names) == 0L))
    return(matrix(numeric(0), nrow = nrow(data), ncol = 0L))
  if (!is.null(legacy_arg_name))
    warning("Argument '", legacy_arg_name, "' is deprecated. Use the formula interface (e.g., formula_bin = ~ x1 + x2).", call. = FALSE)

  if (is.character(formula_or_names)) {
    formula_obj <- stats::reformulate(formula_or_names, intercept = TRUE)
  } else if (inherits(formula_or_names, "formula")) {
    if (length(formula_or_names) == 2L) formula_obj <- formula_or_names
    else stop(label, " must be a one-sided formula (e.g., ~ x1 + x2).")
  } else {
    stop(label, " must be a formula or a character vector.")
  }

  mm <- stats::model.matrix(formula_obj, data = data)
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  if (ncol(mm) > 0) {
    missing_cols <- setdiff(all.vars(formula_obj), names(data))
    if (length(missing_cols)) stop(label, ": column(s) not found: ", paste(missing_cols, collapse = ", "))
  }
  mm
}

#' Parse a random-effects formula of the form ~ 1 | Subject
#'
#' Only random intercepts are supported. Random slopes are rejected.
#'
#' @param random_formula A formula like ~ 1 | Subject.
#' @param legacy_arg Optional deprecated argument (subject_ind).
#' @return A character string with the subject variable name.
#' @keywords internal
#' @noRd
.parse_random_formula <- function(random_formula, legacy_arg = NULL) {
  if (!is.null(legacy_arg)) {
    warning("Argument '", legacy_arg, "' is deprecated. Use: random = ~ 1 | Subject", call. = FALSE)
    return(legacy_arg)
  }
  if (!inherits(random_formula, "formula")) stop("random must be a formula of the form ~ 1 | Subject.")
  rhs <- random_formula[[length(random_formula)]]
  if (length(rhs) < 3L || as.character(rhs[[1L]]) != "|")
    stop("random must be of the form ~ 1 | Subject.")
  # --- Validate LHS == 1 (only random intercepts are supported) ---
  lhs <- rhs[[2L]]
  if (!(is.numeric(lhs) && length(lhs) == 1L && isTRUE(all.equal(as.numeric(lhs), 1))))
    stop("Only random intercepts are supported. Use '~ 1 | Subject'.")
  # ---------------------------------------------------------------
  group_var <- as.character(rhs[[3L]])
  if (length(group_var) != 1L || is.na(group_var))
    stop("random must specify exactly one variable after '|'.")
  group_var
}

#' Fit a zero-inflated (or zero-augmented) random-intercept model
#'
#' Shared engine for two-component models: a presence (logistic) component
#' and a positive (Beta or Vasicek) component. Uses non-adaptive
#' Gauss-Hermite quadrature for the random intercept.
#'
#' @param data A data.frame.
#' @param y Character name of the response column (values in [0,1)).
#' @param formula_bin One-sided formula for the presence component.
#' @param formula_cont One-sided formula for the positive component.
#' @param random Formula ~ 1 | Subject specifying the random intercept.
#' @param logistic_cov Deprecated: character vector of presence covariates.
#' @param positive_cov Deprecated: character vector of positive covariates.
#' @param subject_ind Deprecated: subject column name.
#' @param time_ind Character name of the time column (used for duplicate checks).
#' @param component_wise_test Logical; compute component-wise LRTs?
#' @param quad_n Integer >= 2; number of Gauss-Hermite quadrature points.
#' @param verbose Logical; print progress?
#' @param joint_test NULL, TRUE, or FALSE; compute joint LRTs?
#' @param sd_lower Lower bound for the random-effect SD (> 0).
#' @param start Optional list of starting values.
#' @param control List of control parameters passed to nlminb.
#' @param hessian Logical; compute Hessian-based standard errors?
#' @param model_name Character; name of the model (for messages).
#' @param required_pkgs Character vector of required package names.
#' @param positive_density_fn Function(yy, mm, shape, sz) returning the
#'   log-density of the positive component.
#' @param positive_name Character; name of the positive component.
#' @param validate_shape_fn Function(t2) validating the shape parameter on
#'   the unconstrained scale; returns the natural-scale value or NULL.
#' @param natural_shape_fn Function(val, validate_only, get_jacobian)
#'   handling the shape transformation and its Jacobian.
#' @param start_seeds_fn Function(log_sd_start, q, Y) returning a list of
#'   starting-value vectors for the positive component.
#' @param shape_param_name Character; name of the shape parameter.
#' @return A list with all fitted quantities.
#' @keywords internal
#' @noRd
.fit_zero_inflated_engine <- function(data, y, formula_bin = NULL, formula_cont = NULL, random = NULL,
                                      logistic_cov = NULL, positive_cov = NULL, subject_ind = NULL, time_ind,
                                      component_wise_test, quad_n, verbose, joint_test, sd_lower,
                                      start, control, hessian, model_name, required_pkgs,
                                      positive_density_fn, positive_name, validate_shape_fn,
                                      natural_shape_fn, start_seeds_fn, shape_param_name = positive_name) {
  for (pkg in required_pkgs) if (!requireNamespace(pkg, quietly = TRUE)) stop("Please install package '", pkg, "'.")

  scalar_flag <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
  if (!all(vapply(list(component_wise_test, verbose, hessian), scalar_flag, logical(1)))) stop("Invalid logical argument.")
  if (!is.null(joint_test) && !scalar_flag(joint_test)) stop("joint_test must be NULL, TRUE, or FALSE.")
  if (!is.numeric(quad_n) || length(quad_n) != 1L || !is.finite(quad_n) || quad_n < 2 || quad_n != floor(quad_n)) stop("quad_n must be an integer >= 2.")
  if (!is.numeric(sd_lower) || length(sd_lower) != 1L || !is.finite(sd_lower) || sd_lower <= 0) stop("sd_lower must be > 0.")
  if (!is.data.frame(data)) stop("data must be a data.frame.")

  get_col <- function(nm, label) {
    if (!is.character(nm) || length(nm) != 1L || is.na(nm)) stop(label, " must be a column name (character).")
    if (!nm %in% names(data)) stop(label, ": column '", nm, "' not found.")
    data[[nm]]
  }

  Y <- get_col(y, "y")
  if (!is.null(random)) subject_col_name <- .parse_random_formula(random)
  else if (!is.null(subject_ind)) subject_col_name <- .parse_random_formula(NULL, subject_ind)
  else stop("Provide 'random' (formula) or 'subject_ind' (deprecated).")

  subject_col <- get_col(subject_col_name, "subject_ind")
  time_col <- get_col(time_ind, "time_ind")
  if (!is.numeric(Y)) stop("Column 'y' must be numeric.")

  if (!is.null(formula_bin)) logistic_mat <- .build_design_matrix(formula_bin, data, "formula_bin")
  else if (!is.null(logistic_cov)) logistic_mat <- .build_design_matrix(logistic_cov, data, "logistic_cov", "logistic_cov")
  else logistic_mat <- matrix(numeric(0), nrow = nrow(data), ncol = 0L)

  if (!is.null(formula_cont)) positive_mat <- .build_design_matrix(formula_cont, data, "formula_cont")
  else if (!is.null(positive_cov)) positive_mat <- .build_design_matrix(positive_cov, data, "positive_cov", paste0(positive_name, "_cov"))
  else positive_mat <- matrix(numeric(0), nrow = nrow(data), ncol = 0L)

  Y <- as.vector(Y); n <- length(Y)
  if (!n || any(!is.finite(Y)) || any(Y < 0 | Y >= 1)) stop("Y must be in [0,1).")
  if (!any(Y == 0) || !any(Y > 0)) stop("The model requires both zero and positive responses.")
  if (length(subject_col) != n || length(time_col) != n || anyNA(subject_col) || anyNA(time_col)) stop("Subject and Time cannot contain NA.")
  if (anyDuplicated(data.frame(subject_col, time_col))) stop("Duplicated Subject/Time pairs.")

  subject <- as.character(subject_col)
  group <- match(subject, unique(subject))
  ns <- length(unique(group))
  if (ns < 2L) stop("At least 2 subjects are required.")

  make_design <- function(x, label) {
    x <- as.matrix(x)
    if (!is.numeric(x) || nrow(x) != n || any(!is.finite(x))) stop(label, " must be a finite numeric matrix.")
    if (is.null(colnames(x)) && ncol(x) > 0L) colnames(x) <- paste0("var", seq_len(ncol(x)))
    if (anyDuplicated(colnames(x)) || "(Intercept)" %in% colnames(x)) stop(label, ": use unique names and omit the intercept.")
    a <- cbind("(Intercept)" = 1, x)
    if (qr(a)$rank != ncol(a)) stop(label, " is rank deficient.")
    a
  }

  X <- make_design(logistic_mat, "logistic_cov/formula_bin")
  Z <- make_design(positive_mat, "positive_cov/formula_cont")
  if (qr(Z[Y > 0, , drop = FALSE])$rank != ncol(Z)) stop("Positive-component design matrix is rank deficient.")

  same_design <- identical(colnames(X), colnames(Z)) && isTRUE(all.equal(unname(X), unname(Z), check.attributes = FALSE))
  if (is.null(joint_test)) joint_test <- same_design
  if (joint_test && !same_design) stop("Joint tests require identical covariates.")

  if (!is.list(control)) stop("control must be a list.")
  ctl <- utils::modifyList(list(iter.max = 2000L, eval.max = 5000L, rel.tol = 1e-10, trace = if (verbose) 1L else 0L), control)

  gh <- statmod::gauss.quad(quad_n, kind = "hermite")
  nodes <- sqrt(2) * gh$nodes
  log_weights <- log(gh$weights / sqrt(pi))
  softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
  logsumexp <- function(x) { m <- max(x); if (!is.finite(m)) return(m); m + log(sum(exp(x - m))) }
  integrate_subjects <- function(log_density) {
    sums <- rowsum(log_density, group, reorder = FALSE)
    value <- sum(apply(sweep(sums, 2L, log_weights, "+"), 1L, logsumexp))
    if (!is.finite(value)) Inf else -value
  }

  make_objective <- function(A, positive) {
    force(A); force(positive)
    function(t) {
      sd <- exp(t[1L]); offset <- if (positive) 2L else 1L
      if (!is.finite(sd)) return(Inf)
      eta <- outer(drop(A %*% t[-seq_len(offset)]), sd * nodes, "+")
      if (any(!is.finite(eta))) return(Inf)
      if (!positive) {
        ld <- -softplus(eta); ld[Y > 0, ] <- -softplus(-eta[Y > 0, , drop = FALSE])
      } else {
        shape_param <- validate_shape_fn(t[2L])
        if (is.null(shape_param)) return(Inf)
        ld <- matrix(0, nrow = n, ncol = quad_n)
        mu <- stats::plogis(eta[Y > 0, , drop = FALSE])
        yy <- rep(Y[Y > 0], times = quad_n); mm <- as.vector(mu)
        valid <- mm > 0 & mm < 1; logd <- rep(-Inf, length(mm))
        if (any(valid)) {
          vals <- positive_density_fn(yy[valid], mm[valid], shape_param, sum(valid))
          if (length(vals) != sum(valid) || anyNA(vals) || any(vals == Inf)) return(Inf)
          logd[valid] <- vals
        }
        ld[Y > 0, ] <- matrix(logd, nrow = sum(Y > 0), ncol = quad_n)
      }
      integrate_subjects(ld)
    }
  }

  fit_block <- function(A, positive, seeds) {
    fn <- make_objective(A, positive)
    lower <- c(log(sd_lower), rep(-Inf, length(seeds[[1L]]) - 1L))
    candidates <- lapply(seeds, function(s) tryCatch(stats::nlminb(s, fn, lower = lower, control = ctl), error = function(e) list(par = s, objective = Inf, convergence = 999L, message = conditionMessage(e))))
    values <- vapply(candidates, function(x) x$objective, numeric(1))
    if (!any(is.finite(values))) stop("No finite optimization result.")
    opt <- candidates[[which.min(values)]]
    list(opt = opt, fn = fn, boundary = opt$par[1L] <= log(sd_lower) + 1e-4, candidates = lapply(candidates, function(x) x[c("objective", "convergence", "message")]))
  }

  p <- ncol(X); q <- ncol(Z); log_sd_start <- log(max(1, sd_lower))
  seeds_l <- list(c(log_sd_start, rep(0, p)), c(log_sd_start, stats::qlogis(mean(Y > 0)), rep(0, p - 1L)))
  seeds_v <- start_seeds_fn(log_sd_start, q, Y)

  if (!is.null(start)) {
    if (!is.list(start) || (length(start) && is.null(names(start))) || any(!names(start) %in% c("logistic", positive_name))) stop("Invalid 'start'.")
    for (nm in c("logistic", positive_name)) {
      s <- start[[nm]]; if (is.null(s)) next
      expected <- if (nm == "logistic") p + 1L else q + 2L
      if (!is.numeric(s) || length(s) != expected || any(!is.finite(s)) || s[1L] < sd_lower) stop("Invalid starting vector: ", nm)
      s[1L] <- log(s[1L])
      if (nm == positive_name) { s[2L] <- natural_shape_fn(s[2L], validate_only = TRUE); seeds_v <- c(list(s), seeds_v) }
      else seeds_l <- c(list(s), seeds_l)
    }
  }

  full_l <- fit_block(X, FALSE, seeds_l)
  full_v <- fit_block(Z, TRUE, seeds_v)

  covariance_block <- function(fit, positive, names_coef) {
    t <- fit$opt$par; natural <- t; natural[1L] <- exp(t[1L]); jac <- rep(1, length(t)); jac[1L] <- natural[1L]
    if (positive) { res_shape <- natural_shape_fn(t[2L], get_jacobian = TRUE); natural[2L] <- res_shape$val; jac[2L] <- res_shape$jac }
    names(natural) <- if (positive) c("s2", shape_param_name, paste0("beta_", names_coef)) else c("s1", paste0("alpha_", names_coef))
    V <- matrix(NA_real_, length(t), length(t), dimnames = list(names(natural), names(natural)))
    status <- "Not requested"
    if (hessian) {
      status <- "Unavailable: nonconvergence or SD at lower bound"
      if (fit$opt$convergence == 0L && !fit$boundary) {
        ans <- tryCatch({ H <- numDeriv::hessian(fit$fn, t); H <- (H + t(H)) / 2; if (any(!is.finite(H))) stop("Nonfinite Hessian"); chol2inv(chol(H)) * outer(jac, jac) }, error = function(e) e)
        if (inherits(ans, "error")) status <- conditionMessage(ans) else { V[,] <- ans; status <- "OK" }
      }
    }
    list(estimate = natural, V = V, SE = sqrt(diag(V)), status = status)
  }

  cv_l <- covariance_block(full_l, FALSE, colnames(X))
  cv_v <- covariance_block(full_v, TRUE, colnames(Z))

  null_l <- vector("list", p); null_v <- vector("list", q)
  run_tests <- function(A, positive, full, seeds, indices) {
    offset <- if (positive) 2L else 1L; ans <- vector("list", ncol(A))
    for (j in indices) {
      if (verbose) message("LRT: ", if (positive) paste0(positive_name, " / ") else "presence / ", colnames(A)[j])
      reduced_seeds <- lapply(c(list(full$opt$par), seeds), function(s) s[-(offset + j)])
      ans[[j]] <- tryCatch(fit_block(A[, -j, drop = FALSE], positive, reduced_seeds), error = function(e) list(opt = list(objective = Inf, convergence = 999L, message = conditionMessage(e)), boundary = NA))
    }
    ans
  }

  idx_l <- if (component_wise_test) seq_len(p) else if (joint_test) seq_len(p)[-1L] else integer()
  idx_v <- if (component_wise_test) seq_len(q) else if (joint_test) seq_len(q)[-1L] else integer()
  if (length(idx_l)) null_l <- run_tests(X, FALSE, full_l, seeds_l, idx_l)
  if (length(idx_v)) null_v <- run_tests(Z, TRUE, full_v, seeds_v, idx_v)

  lrt_stat <- function(full, null) {
    if (is.null(null) || full$opt$convergence != 0L || null$opt$convergence != 0L) return(NA_real_)
    value <- 2 * (null$opt$objective - full$opt$objective)
    if (!is.finite(value) || value < -1e-6) return(NA_real_)
    max(0, value)
  }

  stat_l <- vapply(null_l, function(x) lrt_stat(full_l, x), numeric(1))
  stat_v <- vapply(null_v, function(x) lrt_stat(full_v, x), numeric(1))

  fixed_table <- function(cv, offset, labels, stat) {
    ind <- seq_along(labels) + offset; b <- unname(cv$estimate[ind]); se <- unname(cv$SE[ind])
    data.frame(Estimate = b, SE = se, Pvalue = if (component_wise_test) stats::pchisq(stat, 1, lower.tail = FALSE) else rep(NA_real_, length(labels)), Wald_Pvalue = 2 * stats::pnorm(-abs(b / se)), row.names = labels)
  }

  joint_stat <- if (joint_test) stats::setNames(stat_l[-1L] + stat_v[-1L], colnames(X)[-1L]) else NULL
  joint_p <- if (joint_test) stats::pchisq(joint_stat, 2, lower.tail = FALSE) else NULL

  estimates <- c(cv_l$estimate, cv_v$estimate); k <- length(estimates)
  V <- matrix(0, k, k, dimnames = list(names(estimates), names(estimates)))
  il <- seq_along(cv_l$estimate); iv <- length(il) + seq_along(cv_v$estimate)
  V[il, il] <- cv_l$V; V[iv, iv] <- cv_v$V

  random_table <- data.frame(Estimate = c(cv_l$estimate[1L], cv_l$estimate[1L]^2, cv_v$estimate[1L], cv_v$estimate[1L]^2),
                             SE = c(cv_l$SE[1L], 2 * cv_l$estimate[1L] * cv_l$SE[1L], cv_v$SE[1L], 2 * cv_v$estimate[1L] * cv_v$SE[1L]),
                             row.names = c("Presence_SD", "Presence_variance", "Positive_SD", "Positive_variance"))

  # ------------------------------------------------------------------
  # FIT STATISTICS (aligned with PROC NLMIXED, SAS/STAT 14.2)
  #   AIC  = 2f + 2p
  #   BIC  = 2f + p * log(s)
  #   where f = -logLik, p = #params, n = #obs, s = #subjects
  # ------------------------------------------------------------------
  minus2ll_bin   <- 2 * full_l$opt$objective
  minus2ll_cont  <- 2 * full_v$opt$objective
  minus2ll_total <- minus2ll_bin + minus2ll_cont

  # p and q already include the intercept.
  k_bin  <- p + 1L   # p fixed coefs + 1 random-effect SD
  k_cont <- q + 2L   # q fixed coefs + 1 shape + 1 random-effect SD
  k      <- k_bin + k_cont

  # --- Presence component (all n obs / ns subjects) ---
  n_bin  <- n
  ns_bin <- ns

  # --- Positive component (only obs > 0 and subjects with at least 1 such obs) ---
  pos_mask  <- Y > 0
  n_cont    <- sum(pos_mask)
  ns_cont   <- length(unique(subject[pos_mask]))

  # Joint statistics (what PROC NLMIXED reports)
  AIC_total  <- minus2ll_total + 2 * k
  BIC_total  <- minus2ll_total + k * log(ns)

  # Component-wise statistics (informational, NOT additive)
  AIC_bin    <- minus2ll_bin  + 2 * k_bin
  BIC_bin    <- minus2ll_bin  + k_bin  * log(ns_bin)

  AIC_cont   <- minus2ll_cont + 2 * k_cont
  BIC_cont   <- minus2ll_cont + k_cont * log(ns_cont)

  diag_row <- function(fit, label, hess = NA_character_) {
    data.frame(Fit = label, Convergence = fit$opt$convergence, Boundary = fit$boundary, LogLik = -fit$opt$objective, Hessian = hess, Message = paste(fit$opt$message, collapse = "; "))
  }

  diagnostics <- rbind(diag_row(full_l, "presence_full", cv_l$status), diag_row(full_v, paste0(positive_name, "_full"), cv_v$status))
  for (j in idx_l) diagnostics <- rbind(diagnostics, diag_row(null_l[[j]], paste0("presence_without_", colnames(X)[j])))
  for (j in idx_v) diagnostics <- rbind(diagnostics, diag_row(null_v[[j]], paste0(positive_name, "_without_", colnames(Z)[j])))

  if (any(diagnostics$Convergence != 0L)) warning("Some optimizations did not converge. Inspect $diagnostics.")
  if (any(diagnostics$Boundary, NA, na.rm = TRUE)) warning("An SD reached the lower bound; Wald/LRT inference requires caution.")
  if (hessian && (cv_l$status != "OK" || cv_v$status != "OK")) warning("Some standard errors are unavailable. Inspect $diagnostics.")
  if (anyNA(stat_l[idx_l]) || anyNA(stat_v[idx_v])) warning("Some LRTs are unavailable.")

  list(full_l = full_l, full_v = full_v, cv_l = cv_l, cv_v = cv_v, stat_l = stat_l, stat_v = stat_v, fixed_table = fixed_table,
       estimates = estimates, vcov = V, random_table = random_table,
       stats_bin = data.frame(N = n_bin, N_subjects = ns_bin, K = k_bin, Neg2LogLik = minus2ll_bin, AIC = AIC_bin, BIC = BIC_bin),
       stats_cont = data.frame(N = n_cont, N_subjects = ns_cont, K = k_cont, Neg2LogLik = minus2ll_cont, AIC = AIC_cont, BIC = BIC_cont),
       stats_total = data.frame(N_obs = n, N_subjects = ns, K = k, Neg2LogLik = minus2ll_total, AIC = AIC_total, BIC = BIC_total),
       loglikelihood = -minus2ll_total / 2, joint_p = joint_p, joint_statistic = joint_stat, diagnostics = diagnostics,
       n = n, ns = ns, k = k, sd_lower = sd_lower, quad_n = quad_n, idx_l = idx_l, idx_v = idx_v,
       null_l = null_l, null_v = null_v, logistic_names = colnames(X), positive_names = colnames(Z)
  )
}

