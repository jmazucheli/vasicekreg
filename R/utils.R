.check_unit_interval <- function(x, name, closed = FALSE) {
    if (!is.numeric(x) || anyNA(x)) {
        stop(sprintf("'%s' must be a numeric vector without missing values.", name),
             call. = FALSE)
    }

    valid <- if (closed) {
        x >= 0 & x <= 1
    } else {
        x > 0 & x < 1
    }

    if (!all(valid)) {
        interval <- if (closed) "[0, 1]" else "(0, 1)"
        stop(sprintf("'%s' must contain values in %s.", name, interval),
             call. = FALSE)
    }

    invisible(TRUE)
}

.check_probability <- function(p, log.p = FALSE) {
    if (!is.numeric(p) || anyNA(p)) {
        stop("'p' must be a numeric vector without missing values.",
             call. = FALSE)
    }

    valid <- if (isTRUE(log.p)) {
        p <= 0
    } else {
        p >= 0 & p <= 1
    }

    if (!all(valid)) {
        scale <- if (isTRUE(log.p)) "(-Inf, 0]" else "[0, 1]"
        stop(sprintf("'p' must contain values in %s.", scale),
             call. = FALSE)
    }

    invisible(TRUE)
}

.check_scalar_logical <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
        stop(sprintf("'%s' must be TRUE or FALSE.", name), call. = FALSE)
    }
    invisible(TRUE)
}

.n_random <- function(n) {
    if (length(n) > 1L) {
        return(length(n))
    }

    if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
        !is.finite(n) || n < 0 || n != floor(n)) {
        stop("'n' must be a non-negative integer or a vector.", call. = FALSE)
    }

    as.integer(n)
}
