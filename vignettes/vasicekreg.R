## ----setup, include=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  fig.height = 4.8,
  fig.width = 7.0,
  message = FALSE,
  warning = FALSE
)

library(gamlss)
library(vasicekreg)

control <- gamlss.control(n.cyc = 200, trace = FALSE)
options(width = 90)
set.seed(2026)

## ----helpers----------------------------------------------------------------------------
model_fit_table <- function(models, n) {
  rows <- lapply(names(models), function(label) {
    object <- models[[label]]
    data.frame(
      model = label,
      family = as.character(object$family[1L]),
      parameters = object$df.fit,
      logLik = -object$G.deviance / 2,
      AIC = gamlss::GAIC(object, k = 2),
      BIC = gamlss::GAIC(object, k = log(n)),
      converged = isTRUE(object$converged),
      row.names = NULL
    )
  })
  do.call(rbind, rows)
}

coefficient_table <- function(object, parameters) {
  rows <- lapply(parameters, function(parameter) {
    estimate <- stats::coef(object, what = parameter)
    data.frame(
      parameter = parameter,
      term = names(estimate),
      estimate = unname(estimate),
      row.names = NULL
    )
  })
  do.call(rbind, rows)
}

fitted_summary <- function(values) {
  rows <- lapply(names(values), function(quantity) {
    x <- values[[quantity]]
    data.frame(
      quantity = quantity,
      minimum = min(x),
      first_quartile = unname(stats::quantile(x, 0.25)),
      median = stats::median(x),
      mean = mean(x),
      third_quartile = unname(stats::quantile(x, 0.75)),
      maximum = max(x),
      row.names = NULL
    )
  })
  do.call(rbind, rows)
}

response_profile <- function(x) {
  c(
    observations = length(x),
    zero = sum(x == 0),
    interior = sum(x > 0 & x < 1),
    one = sum(x == 1)
  )
}

## ----load-data--------------------------------------------------------------------------
data("bodyfat", package = "vasicekreg")
data("transport", package = "vasicekreg")
data("trees", package = "vasicekreg")
data("aep", package = "vasicekreg")

bodyfat_responses <- c("ARMS", "LEGS", "BODY", "ANDROID", "GYNECOID")
aep$inappropriate <- with(aep, noinap / los)

profiles <- rbind(
  do.call(rbind, lapply(bodyfat_responses, function(response) {
    data.frame(
      data = "bodyfat",
      response = response,
      t(response_profile(bodyfat[[response]])),
      row.names = NULL
    )
  })),
  data.frame(
    data = "transport",
    response = "propbiked",
    t(response_profile(transport$propbiked)),
    row.names = NULL
  ),
  data.frame(
    data = "trees",
    response = "prop",
    t(response_profile(trees$prop)),
    row.names = NULL
  ),
  data.frame(
    data = "aep",
    response = "noinap / los",
    t(response_profile(aep$inappropriate)),
    row.names = NULL
  )
)

knitr::kable(profiles, caption = "Observed support of the bounded responses.")

## ----bodyfat-prepare--------------------------------------------------------------------
bodyfat_analysis <- within(bodyfat, {
  AGE_centered <- AGE - mean(AGE)
  BMI_centered <- BMI - mean(BMI)
  SEX <- factor(SEX, levels = c(1, 2), labels = c("female", "male"))
  IPAQ <- factor(
    IPAQ,
    levels = c(0, 1, 2),
    labels = c("sedentary", "insufficiently_active", "active")
  )
})

bodyfat_ranges <- data.frame(
  response = bodyfat_responses,
  minimum = vapply(bodyfat_analysis[bodyfat_responses], min, numeric(1)),
  maximum = vapply(bodyfat_analysis[bodyfat_responses], max, numeric(1)),
  row.names = NULL
)
knitr::kable(bodyfat_ranges, digits = 3,
             caption = "Ranges of the five body-fat proportions.")

## ----bodyfat-mean-fits------------------------------------------------------------------
fit_bodyfat_mean <- setNames(
  lapply(bodyfat_responses, function(response) {
    mu_formula <- stats::reformulate(
      c("AGE_centered", "BMI_centered", "SEX", "IPAQ"),
      response = response
    )
    gamlss(
      formula = mu_formula,
      sigma.formula = ~ 1,
      family = NVASIM(),
      data = bodyfat_analysis,
      control = control
    )
  }),
  bodyfat_responses
)

bodyfat_mean_statistics <- model_fit_table(
  fit_bodyfat_mean,
  n = nrow(bodyfat_analysis)
)
knitr::kable(
  bodyfat_mean_statistics,
  digits = 3,
  caption = "Normal-kernel Vasicek mean regressions for the body-fat responses."
)

## ----bodyfat-mean-coefficients----------------------------------------------------------
knitr::kable(
  coefficient_table(fit_bodyfat_mean[["ARMS"]], c("mu", "sigma")),
  digits = 4,
  caption = "Coefficient estimates for the ARMS mean-regression model."
)

## ----bodyfat-quantile-fits--------------------------------------------------------------
quantile_level <- 0.50

fam_arms_nq <- NVASIQ(quantile = quantile_level)
fam_arms_lq <- LVASIQ(quantile = quantile_level)
fam_arms_hq <- HVASIQ(quantile = quantile_level)

fit_arms_nq <- gamlss(
  ARMS ~ AGE_centered + BMI_centered + SEX + IPAQ,
  sigma.formula = ~ 1,
  family = fam_arms_nq,
  data = bodyfat_analysis,
  control = control
)

fit_arms_lq <- gamlss(
  ARMS ~ AGE_centered + BMI_centered + SEX + IPAQ,
  sigma.formula = ~ 1,
  family = fam_arms_lq,
  data = bodyfat_analysis,
  control = control
)

fit_arms_hq <- gamlss(
  ARMS ~ AGE_centered + BMI_centered + SEX + IPAQ,
  sigma.formula = ~ 1,
  family = fam_arms_hq,
  data = bodyfat_analysis,
  control = control
)

arms_models <- c(
  list(NVASIM_mean = fit_bodyfat_mean[["ARMS"]]),
  list(
    NVASIQ_median = fit_arms_nq,
    LVASIQ_median = fit_arms_lq,
    HVASIQ_median = fit_arms_hq
  )
)

knitr::kable(
  model_fit_table(arms_models, n = nrow(bodyfat_analysis)),
  digits = 3,
  caption = "Likelihood-based summaries for the ARMS models."
)

## ----bodyfat-fitted-summary-------------------------------------------------------------
knitr::kable(
  fitted_summary(list(
    NVASIM_conditional_mean = fitted(arms_models$NVASIM_mean, what = "mu"),
    NVASIQ_conditional_median = fitted(arms_models$NVASIQ_median, what = "mu"),
    LVASIQ_conditional_median = fitted(arms_models$LVASIQ_median, what = "mu"),
    HVASIQ_conditional_median = fitted(arms_models$HVASIQ_median, what = "mu")
  )),
  digits = 4,
  caption = "Summaries of fitted means and medians for ARMS."
)

## ----transport-fit----------------------------------------------------------------------
transport_analysis <- within(transport, {
  gender <- stats::relevel(factor(gender), ref = "F")
  status <- stats::relevel(factor(status), ref = "faculty")
  parking_centered <- parking - mean(parking)
  distance_centered <- distance - mean(distance)
})

fit_transport <- gamlss(
  propbiked ~ gender + parking_centered + status,
  sigma.formula = ~ 1,
  nu.formula = ~ gender + distance_centered,
  family = ZANVASIM(),
  data = transport_analysis,
  control = control
)

knitr::kable(
  model_fit_table(
    list(ZANVASIM = fit_transport),
    n = nrow(transport_analysis)
  ),
  digits = 3,
  caption = "Likelihood-based summary for the transport model."
)

knitr::kable(
  coefficient_table(fit_transport, c("mu", "sigma", "nu")),
  digits = 4,
  caption = "Coefficient estimates for the zero-augmented transport model."
)

## ----transport-fitted-------------------------------------------------------------------
transport_mu <- fitted(fit_transport, what = "mu")
transport_nu <- fitted(fit_transport, what = "nu")
transport_marginal_mean <- (1 - transport_nu) * transport_mu

knitr::kable(
  fitted_summary(list(
    positive_component_mean = transport_mu,
    probability_zero = transport_nu,
    marginal_mean = transport_marginal_mean
  )),
  digits = 4,
  caption = "Fitted quantities from the transport model."
)

## ----trees-fit--------------------------------------------------------------------------
trees_analysis <- within(trees, {
  precip_centered <- precip - mean(precip)
  wind_centered <- wind - mean(wind)
})

fit_trees <- gamlss(
  prop ~ pest + fertilization + precip_centered + wind_centered,
  sigma.formula = ~ 1,
  nu.formula = ~ wind_centered,
  family = OANVASIM(),
  data = trees_analysis,
  control = control
)

knitr::kable(
  model_fit_table(
    list(OANVASIM = fit_trees),
    n = nrow(trees_analysis)
  ),
  digits = 3,
  caption = "Likelihood-based summary for the tree-survival model."
)

knitr::kable(
  coefficient_table(fit_trees, c("mu", "sigma", "nu")),
  digits = 4,
  caption = "Coefficient estimates for the one-augmented tree-survival model."
)

## ----trees-fitted-----------------------------------------------------------------------
trees_mu <- fitted(fit_trees, what = "mu")
trees_nu <- fitted(fit_trees, what = "nu")
trees_marginal_mean <- trees_nu + (1 - trees_nu) * trees_mu

knitr::kable(
  fitted_summary(list(
    continuous_component_mean = trees_mu,
    probability_one = trees_nu,
    marginal_mean = trees_marginal_mean
  )),
  digits = 4,
  caption = "Fitted quantities from the tree-survival model."
)

## ----aep-fit----------------------------------------------------------------------------
aep_analysis <- within(aep, {
  sex <- stats::relevel(factor(sex), ref = "1")
  ward <- stats::relevel(factor(ward), ref = "1")
  year <- stats::relevel(factor(year), ref = "88")
})

fit_aep <- gamlss(
  inappropriate ~ sex + ward + year + age + loglos,
  sigma.formula = ~ loglos,
  nu.formula = ~ loglos,
  tau.formula = ~ loglos,
  family = ZOANVASIM(),
  data = aep_analysis,
  control = control
)

knitr::kable(
  model_fit_table(
    list(ZOANVASIM = fit_aep),
    n = nrow(aep_analysis)
  ),
  digits = 3,
  caption = "Likelihood-based summary for the hospital-stay model."
)

knitr::kable(
  coefficient_table(fit_aep, c("mu", "sigma", "nu", "tau")),
  digits = 4,
  caption = "Coefficient estimates for the zero-and-one-augmented hospital-stay model."
)

## ----aep-fitted-------------------------------------------------------------------------
aep_mu <- fitted(fit_aep, what = "mu")
aep_sigma <- fitted(fit_aep, what = "sigma")
aep_nu <- fitted(fit_aep, what = "nu")
aep_tau <- fitted(fit_aep, what = "tau")

aep_probability_zero <- aep_nu
aep_probability_one <- (1 - aep_nu) * aep_tau
aep_probability_continuous <- (1 - aep_nu) * (1 - aep_tau)
aep_marginal_mean <- (1 - aep_nu) * (
  aep_tau + (1 - aep_tau) * aep_mu
)

stopifnot(all.equal(
  aep_probability_zero + aep_probability_one + aep_probability_continuous,
  rep(1, nrow(aep_analysis)),
  tolerance = 1e-8
))

knitr::kable(
  fitted_summary(list(
    continuous_component_mean = aep_mu,
    shape = aep_sigma,
    probability_zero = aep_probability_zero,
    probability_one = aep_probability_one,
    probability_continuous = aep_probability_continuous,
    marginal_mean = aep_marginal_mean
  )),
  digits = 4,
  caption = "Fitted quantities from the hospital-stay model."
)

## ----aep-envelope, eval=FALSE-----------------------------------------------------------
# envelope_aep <- vasicek_envelope(
#   object = fit_aep,
#   residual = c("quantile", "cox-snell"),
#   nsim = 500,
#   level = 0.95,
#   envelope = "quantile",
#   seed = 2026,
#   data = aep_analysis
# )
# 
# old_par <- graphics::par(no.readonly = TRUE)
# graphics::par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
# plot(envelope_aep, which = "quantile", pch = 19, cex = 0.55)
# plot(envelope_aep, which = "cox-snell", pch = 19, cex = 0.55)
# graphics::par(old_par)

## ----session-info-----------------------------------------------------------------------
sessionInfo()

