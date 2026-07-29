with_global_tau <- function(value, code) {
    had_tau <- exists("tau", envir = .GlobalEnv, inherits = FALSE)
    if (had_tau) {
        old_tau <- get("tau", envir = .GlobalEnv, inherits = FALSE)
    }

    on.exit({
        if (had_tau) {
            assign("tau", old_tau, envir = .GlobalEnv)
        } else if (exists("tau", envir = .GlobalEnv, inherits = FALSE)) {
            rm("tau", envir = .GlobalEnv)
        }
    }, add = TRUE)

    assign("tau", value, envir = .GlobalEnv)
    force(code)
}
