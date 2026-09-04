#' Simulated envelopes for residuals from Vasicek GAMLSS models
#'
#' @description
#' Constructs pointwise simulated envelopes for normalized randomized quantile
#' residuals and generalized Cox--Snell residuals from a Vasicek model fitted by
#' \code{gamlss()}. Each simulated response is generated under the fitted
#' model, the model is re-estimated, and the residuals from the re-estimated
#' model are ordered before the envelope is computed.
#'
#' @param object A fitted object of class \code{"gamlss"} using one of the
#'   Vasicek families supplied by \pkg{vasicekreg}.
#' @param residual Character vector selecting \code{"quantile"},
#'   \code{"cox-snell"}, or both.
#' @param nsim Number of successfully re-fitted simulated samples.
#' @param level Pointwise coverage probability used when
#'   \code{envelope = "quantile"}.
#' @param envelope Either \code{"quantile"}, for percentile limits, or
#'   \code{"minmax"}, for the minimum and maximum at each order position.
#' @param seed Optional integer seed. The previous random-number state is
#'   restored on exit.
#' @param data Optional data frame used in the original fit. It is normally
#'   recovered from \code{object$call$data}.
#' @param max.attempts Maximum number of simulation and re-fitting attempts.
#'   The default is five times \code{nsim}.
#' @param refit Optional function with arguments \code{object}, \code{y}, and
#'   \code{data}. It must return a re-fitted \code{"gamlss"} object. This is
#'   useful for transformed responses or nonstandard fitting calls.
#' @param simulate Optional function with the single argument \code{object}
#'   that returns one simulated response vector. By default, the random
#'   generator associated with the fitted Vasicek family is used.
#' @param verbose Logical; if \code{TRUE}, reports progress and failed fits.
#' @param x An object returned by \code{vasicek_envelope()}.
#' @param which Residual type to be plotted.
#' @param show.mean Logical; if \code{TRUE}, draws the pointwise mean of the
#'   ordered simulated residuals.
#' @param xlab,ylab,main Graphical labels.
#' @param envelope.col,mean.col,reference.col,point.col Graphical colors.
#' @param ... Further arguments. For the plot method, they are passed to
#'   \code{points()}.
#'
#' @details
#' Let \eqn{U_i} denote the probability integral transform used by the
#' normalized randomized quantile residual. The two residuals are
#' \deqn{r_i^q=\Phi^{-1}(U_i)}
#' and
#' \deqn{r_i^{CS}=-\log(1-U_i).}
#' For augmented families, \eqn{U_i} is randomized over the relevant jump of
#' the fitted distribution at zero or one, so the reported residuals inherit
#' whatever randomization \code{residuals()} applies for those families; the
#' Cox--Snell residual is then also randomized at boundary observations and
#' has an \eqn{\operatorname{Exp}(1)} reference distribution under a correctly
#' specified model. The Cox--Snell residual is obtained from the quantile
#' residual on the log survival scale, \eqn{-\log\{\Pr(Z>r_i^q)\}}, which is
#' numerically stable in the upper tail and preserves the ordering of the
#' quantile residuals.
#'
#' Each residual is displayed as a full quantile--quantile plot of the ordered
#' residuals against the corresponding theoretical quantiles (normal or
#' exponential), not as a half-normal plot. The envelope is pointwise, not
#' simultaneous: even under a correctly specified model a fraction of points
#' is expected to fall outside the band. With \code{envelope = "quantile"},
#' the two tail probabilities are equal and sum to \eqn{1-\mathrm{level}}.
#' With \code{envelope = "minmax"}, the limits are the minimum and maximum at
#' each order position, following the construction used by Zhao et al.
#'
#' The pointwise mean of the ordered simulated residuals is the reference
#' calibrated to the fitted model and to the sample size. The identity line
#' drawn by the plot method is a theoretical idealization and, for Cox--Snell
#' residuals, may separate from the mean curve in the upper tail at finite
#' \eqn{n}; in that region the mean curve is the more reliable reference.
#'
#' Failed or nonconverged re-fits are discarded and replaced until
#' \code{nsim} successful samples are obtained or \code{max.attempts} is
#' reached. The simulation uses fitted values for every distribution
#' parameter, so covariate-dependent parameters are retained. Automatic
#' simulation and refitting of \code{NVASIQ}, \code{LVASIQ}, and
#' \code{HVASIQ} recover the fixed level embedded in the fitted object; no
#' global quantile-level variable is consulted. Consequently, fitted models
#' at different quantile levels can be used in the same R session safely.
#'
#' Automatic
#' refitting assumes that \code{data} contains exactly the observations used
#' by the original fit; if rows were dropped for missing values or via
#' \code{subset}, supply a matching \code{data} or a custom \code{refit}.
#'
#' @return
#' An object of class \code{"vasicek_envelope"}. Its \code{results} component
#' contains, for each requested residual, the theoretical order statistics,
#' ordered observed residuals, pointwise lower, mean and upper curves, and the
#' matrix of ordered simulated residuals. The object also records the number
#' of successful simulations, attempts and failures.
#'
#' @references
#' Moral, R. A., Hinde, J., and Demetrio, C. G. B. (2017).
#' Half-normal plots and overdispersed models in R: The hnp package.
#' \emph{Journal of Statistical Software}, \bold{81}(10), 1--23.
#' \doi{10.18637/jss.v081.i10}
#'
#' Zhao, Y., Lee, A. H., Yau, K. K. W., and McLachlan, G. J. (2011).
#' Assessing the adequacy of Weibull survival models: A simulated envelope
#' approach. \emph{Journal of Applied Statistics}, \bold{38}, 2089--2097.
#' \doi{10.1080/02664763.2010.545115}
#'
#' @examples
#' \dontrun{
#' library(gamlss)
#' set.seed(123)
#' dat <- data.frame(y = rNVASIM(100, mu = 0.55, sigma = 0.30))
#' fit <- gamlss(y ~ 1, sigma.formula = ~ 1, family = NVASIM(),
#'               data = dat, trace = FALSE)
#'
#' env <- vasicek_envelope(fit, nsim = 200, seed = 123)
#' plot(env, which = "quantile")
#' plot(env, which = "cox-snell")
#' }
#'
#' @export
vasicek_envelope <- function(
    object,
    residual = c("quantile", "cox-snell"),
    nsim = 200L,
    level = 0.95,
    envelope = c("quantile", "minmax"),
    seed = NULL,
    data = NULL,
    max.attempts = 5L * nsim,
    refit = NULL,
    simulate = NULL,
    verbose = interactive()
) {
    if (!inherits(object, "gamlss")) {
        stop("'object' must inherit from class 'gamlss'.", call. = FALSE)
    }
    supported <- c(
        "NVASIM", "NVASIQ", "LVASIQ", "HVASIQ", "ZANVASIM", "OANVASIM",
        "ZOANVASIM"
    )
    family_name <- as.character(object$family[1L])
    if (!family_name %in% supported) {
        stop("The fitted family is not supported by vasicek_envelope().",
             call. = FALSE)
    }

    residual <- unique(match.arg(
        residual, c("quantile", "cox-snell"), several.ok = TRUE
    ))
    envelope <- match.arg(envelope)
    nsim <- .envelope_positive_integer(nsim, "nsim")
    max.attempts <- .envelope_positive_integer(max.attempts, "max.attempts")
    if (max.attempts < nsim) {
        stop("'max.attempts' must be at least 'nsim'.", call. = FALSE)
    }
    if (length(level) != 1L || !is.finite(level) ||
        level <= 0 || level >= 1) {
        stop("'level' must be a single number strictly between 0 and 1.",
             call. = FALSE)
    }
    if (!is.null(refit) && !is.function(refit)) {
        stop("'refit' must be NULL or a function.", call. = FALSE)
    }
    if (!is.null(simulate) && !is.function(simulate)) {
        stop("'simulate' must be NULL or a function.", call. = FALSE)
    }

    if (!is.null(seed)) {
        if (length(seed) != 1L || !is.finite(seed)) {
            stop("'seed' must be NULL or one finite number.", call. = FALSE)
        }
        had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
        on.exit({
            if (had_seed) {
                assign(".Random.seed", old_seed, envir = .GlobalEnv)
            } else if (exists(".Random.seed", envir = .GlobalEnv,
                              inherits = FALSE)) {
                rm(".Random.seed", envir = .GlobalEnv)
            }
        }, add = TRUE)
        set.seed(as.integer(seed))
    }

    observed_quantile <- as.numeric(stats::residuals(
        object, what = "z-scores"
    ))
    if (!length(observed_quantile) || any(!is.finite(observed_quantile))) {
        stop("The fitted object does not contain finite quantile residuals.",
             call. = FALSE)
    }
    n <- length(observed_quantile)

    if (is.null(data) && is.null(refit)) {
        data <- .envelope_model_data(object)
    }
    if (is.null(simulate)) {
        simulate <- .envelope_simulate_response
    }
    if (is.null(refit)) {
        refit <- .envelope_refit
    }

    simulated_quantile <- matrix(NA_real_, nrow = n, ncol = nsim)
    successes <- 0L
    attempts <- 0L
    failure_messages <- character(0)

    while (successes < nsim && attempts < max.attempts) {
        attempts <- attempts + 1L
        trial <- tryCatch({
            simulated_y <- as.numeric(simulate(object))
            if (length(simulated_y) != n) {
                stop("The simulated response has an unexpected length.")
            }
            fitted_trial <- refit(object, simulated_y, data)
            if (!inherits(fitted_trial, "gamlss")) {
                stop("The refit function did not return a 'gamlss' object.")
            }
            if (isFALSE(fitted_trial$converged)) {
                stop("The simulated model did not converge.")
            }
            trial_residuals <- as.numeric(stats::residuals(
                fitted_trial, what = "z-scores"
            ))
            if (length(trial_residuals) != n ||
                any(!is.finite(trial_residuals))) {
                stop("The simulated fit produced invalid residuals.")
            }
            sort(trial_residuals)
        }, error = function(e) e)

        if (inherits(trial, "error")) {
            failure_messages <- c(failure_messages, conditionMessage(trial))
            if (verbose) {
                message("Attempt ", attempts, " failed: ",
                        conditionMessage(trial))
            }
        } else {
            successes <- successes + 1L
            simulated_quantile[, successes] <- trial
            if (verbose && (successes == 1L || successes %% 10L == 0L ||
                            successes == nsim)) {
                message("Successful simulated fits: ", successes, "/", nsim)
            }
        }
    }

    if (successes < nsim) {
        stop(
            "Only ", successes, " successful simulated fits were obtained in ",
            attempts, " attempts. Increase 'max.attempts' or inspect the fit.",
            call. = FALSE
        )
    }

    observed <- list(quantile = observed_quantile)
    simulations <- list(quantile = simulated_quantile)
    # Cox--Snell residual r^{CS} = -log(1 - U) with U = Phi(r^q). It is computed
    # on the log survival scale, -log{ P(Z > r^q) }, rather than as
    # -log1p(-pnorm(r^q)). The direct form suffers from catastrophic
    # cancellation in the upper tail (pnorm saturates at 1), which caps the
    # residual; the log-scale form is stable and, being strictly increasing in
    # r^q, keeps the column-wise ordering of the sorted quantile residuals.
    observed[["cox-snell"]] <- -stats::pnorm(
        observed_quantile, lower.tail = FALSE, log.p = TRUE
    )
    simulations[["cox-snell"]] <- -stats::pnorm(
        simulated_quantile, lower.tail = FALSE, log.p = TRUE
    )

    probability_points <- stats::ppoints(n)
    results <- lapply(residual, function(kind) {
        simulated <- simulations[[kind]]
        if (envelope == "quantile") {
            alpha <- (1 - level) / 2
            limits <- apply(
                simulated, 1L, stats::quantile,
                probs = c(alpha, 1 - alpha), names = FALSE, type = 8
            )
            lower <- limits[1L, ]
            upper <- limits[2L, ]
        } else {
            lower <- apply(simulated, 1L, min)
            upper <- apply(simulated, 1L, max)
        }
        theoretical <- if (kind == "quantile") {
            stats::qnorm(probability_points)
        } else {
            stats::qexp(probability_points)
        }
        list(
            theoretical = theoretical,
            observed = sort(observed[[kind]]),
            lower = lower,
            mean = rowMeans(simulated),
            upper = upper,
            simulated = simulated
        )
    })
    names(results) <- residual

    structure(
        list(
            call = match.call(), family = family_name, residual = residual,
            nsim = nsim, level = level, envelope = envelope,
            attempts = attempts, failures = attempts - successes,
            failure.messages = unique(failure_messages), results = results
        ),
        class = "vasicek_envelope"
    )
}

.envelope_positive_integer <- function(x, name) {
    if (length(x) != 1L || !is.finite(x) || x < 1 || x != floor(x)) {
        stop("'", name, "' must be a positive integer.", call. = FALSE)
    }
    as.integer(x)
}

.envelope_model_data <- function(object) {
    data_call <- object$call$data
    if (is.null(data_call)) {
        stop("The original data cannot be recovered. Supply 'data' or 'refit'.",
             call. = FALSE)
    }
    formula_environment <- environment(stats::formula(object, what = "mu"))
    data <- tryCatch(
        eval(data_call, envir = formula_environment, enclos = parent.frame()),
        error = function(e) NULL
    )
    if (!is.data.frame(data)) {
        stop("The original data cannot be recovered. Supply 'data' or 'refit'.",
             call. = FALSE)
    }
    data
}

.envelope_simulate_response <- function(object) {
    family_name <- as.character(object$family[1L])
    parameter_map <- list(
        NVASIM = c("mu", "sigma"),
        NVASIQ = c("mu", "sigma"),
        LVASIQ = c("mu", "sigma"),
        HVASIQ = c("mu", "sigma"),
        ZANVASIM = c("mu", "sigma", "nu"),
        OANVASIM = c("mu", "sigma", "nu"),
        ZOANVASIM = c("mu", "sigma", "nu", "tau")
    )
    parameters <- parameter_map[[family_name]]
    random_function <- get(
        paste0("r", family_name), envir = asNamespace("vasicekreg"),
        inherits = FALSE
    )
    arguments <- lapply(parameters, function(parameter) {
        as.numeric(stats::fitted(object, what = parameter))
    })
    names(arguments) <- parameters
    arguments$n <- length(arguments[[1L]])
    quantile_families <- c("NVASIQ", "LVASIQ", "HVASIQ")
    if (family_name %in% quantile_families) {
        arguments$quantile <- .envelope_fixed_quantile(object)
    }
    arguments <- arguments[
        c(
            "n", parameters,
            if (family_name %in% quantile_families) "quantile"
        )
    ]
    do.call(random_function, arguments)
}

.envelope_fixed_quantile <- function(object) {
    rqres <- object$rqres
    if (is.expression(rqres) && length(rqres) == 1L) {
        rqres <- rqres[[1L]]
    }
    if (!is.call(rqres)) {
        stop(
            "The fixed quantile level cannot be recovered from the fitted model.",
            call. = FALSE
        )
    }
    quantile <- as.list(rqres)$quantile
    if (is.null(quantile)) {
        stop(
            "The fixed quantile level is absent from the fitted model. ",
            "Refit it with the current quantile-family implementation.",
            call. = FALSE
        )
    }
    .check_fixed_quantile(quantile)
}

.envelope_embed_quantile_call <- function(object) {
    family_name <- as.character(object$family[1L])
    quantile_families <- c("NVASIQ", "LVASIQ", "HVASIQ")
    if (!family_name %in% quantile_families) {
        return(object)
    }

    quantile <- .envelope_fixed_quantile(object)
    family_call <- object$call$family
    if (is.symbol(family_call) && identical(as.character(family_call), family_name)) {
        object$call$family <- as.call(list(
            as.name(family_name), quantile = quantile
        ))
        return(object)
    }
    if (is.character(family_call) && length(family_call) == 1L &&
        identical(family_call, family_name)) {
        object$call$family <- as.call(list(
            as.name(family_name), quantile = quantile
        ))
        return(object)
    }
    if (is.call(family_call)) {
        family_function <- get(
            family_name, envir = asNamespace("vasicekreg"), inherits = FALSE
        )
        matched_call <- tryCatch(
            match.call(
                definition = family_function,
                call = family_call,
                expand.dots = FALSE
            ),
            error = function(e) NULL
        )
        if (is.null(matched_call)) {
            stop(
                "The family call cannot be reconstructed automatically. ",
                "Supply a custom 'refit' function.",
                call. = FALSE
            )
        }
        matched_call$quantile <- quantile
        object$call$family <- matched_call
        return(object)
    }

    stop(
        "The family call cannot be reconstructed automatically. ",
        "Supply a custom 'refit' function.",
        call. = FALSE
    )
}

.envelope_refit <- function(object, y, data) {
    response <- stats::formula(object, what = "mu")[[2L]]
    if (!is.symbol(response)) {
        stop(
            "Automatic refitting requires an untransformed response name; ",
            "supply a custom 'refit' function.",
            call. = FALSE
        )
    }
    response_name <- as.character(response)
    if (!response_name %in% names(data)) {
        stop("The response variable is absent from 'data'.", call. = FALSE)
    }
    if (nrow(data) != length(y)) {
        stop(
            "'data' has ", nrow(data), " rows but the fitted model used ",
            length(y), " observations. This usually means rows were dropped ",
            "(missing values or 'subset') during the original fit. Supply a ",
            "'data' argument restricted to the complete cases used in the fit, ",
            "or a custom 'refit' function.",
            call. = FALSE
        )
    }
    bootstrap_data <- data
    bootstrap_data[[response_name]] <- y

    # Store formulas and the fixed quantile directly in the copied call. This
    # prevents update.gamlss() from depending on temporary symbols that may
    # have existed only inside the function that created the original fit.
    for (parameter in object$parameters) {
        formula_name <- if (identical(parameter, "mu")) {
            "formula"
        } else {
            paste0(parameter, ".formula")
        }
        object$call[[formula_name]] <- stats::formula(object, what = parameter)
    }
    object <- .envelope_embed_quantile_call(object)

    suppressWarnings(stats::update(
        object,
        data = bootstrap_data,
        control = object$control,
        trace = FALSE
    ))
}

#' @rdname vasicek_envelope
#' @export
print.vasicek_envelope <- function(x, ...) {
    envelope_label <- if (identical(x$envelope, "quantile")) {
        paste0(format(100 * x$level, trim = TRUE),
               "% pointwise quantile")
    } else {
        "pointwise minimum--maximum"
    }

    cat("Simulated residual envelope\n")
    cat("  Family: ", x$family, "\n", sep = "")
    cat("  Residuals: ", paste(x$residual, collapse = ", "), "\n", sep = "")
    cat("  Envelope: ", envelope_label, "\n", sep = "")
    cat("  Successful simulations: ", x$nsim, "\n", sep = "")
    cat("  Attempts: ", x$attempts, "\n", sep = "")
    cat("  Failed fits: ", x$failures, "\n", sep = "")
    invisible(x)
}

#' @rdname vasicek_envelope
#' @export
plot.vasicek_envelope <- function(
    x,
    which = x$residual[1L],
    show.mean = TRUE,
    xlab = NULL,
    ylab = NULL,
    main = "",
    envelope.col = "grey85",
    mean.col = "blue",
    reference.col = "red",
    point.col = "black",
    ...
) {
    which <- match.arg(which, x$residual)
    result <- x$results[[which]]
    if (is.null(xlab)) {
        xlab <- if (which == "quantile") {
            "Standard normal quantiles"
        } else {
            "Standard exponential quantiles"
        }
    }
    if (is.null(ylab)) {
        ylab <- if (which == "quantile") {
            "Ordered randomized quantile residuals"
        } else {
            "Ordered Cox-Snell residuals"
        }
    }
    limits <- range(
        result$observed, result$lower, result$upper, result$theoretical,
        finite = TRUE
    )
    graphics::plot(
        result$theoretical, result$observed, type = "n",
        xlab = xlab, ylab = ylab, main = main, ylim = limits
    )
    graphics::polygon(
        c(result$theoretical, rev(result$theoretical)),
        c(result$lower, rev(result$upper)),
        col = envelope.col, border = NA
    )
    graphics::lines(result$theoretical, result$lower, col = "grey55")
    graphics::lines(result$theoretical, result$upper, col = "grey55")
    graphics::abline(a = 0, b = 1, col = reference.col, lwd = 1.2)
    if (isTRUE(show.mean)) {
        graphics::lines(
            result$theoretical, result$mean, col = mean.col, lwd = 1.5
        )
    }
    graphics::points(
        result$theoretical, result$observed, col = point.col, ...
    )
    invisible(x)
}
