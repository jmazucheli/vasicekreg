# Fitting Vasicek-Type Regression Models to the Data in vasicekreg

## Scope

The `vasicekreg` package provides distribution functions and GAMLSS
families for responses bounded by the unit interval. This vignette is
deliberately application-oriented: its main purpose is to fit models to
all four data sets distributed with the package and to explain what the
fitted parameters mean. Detailed derivations of the densities and a
systematic study of random generation are left to a separate vignette.

The package is built on the Vasicek construction ([Vasicek
2002](#ref-Vasicek2002); [Mazucheli et al.
2022](#ref-MazucheliEtAl2022)) and on the Generalized Additive Models
for Location, Scale and Shape (GAMLSS) framework ([Rigby and
Stasinopoulos 2005](#ref-RigbyStasinopoulos2005); [Stasinopoulos and
Rigby 2007](#ref-StasinopoulosRigby2007)). The logistic and
hyperbolic-secant kernels are motivated by Witzany
([2013](#ref-Witzany2013)) and Fischer et al.
([2017](#ref-FischerEtAl2017)), respectively.

The available regression families are summarized below.

| Family | Kernel | Support | Parameterization | Meaning of `mu` |
|:---|:---|:---|:---|:---|
| `NVASIM` | Standard normal | $`(0,1)`$ | Mean | $`E(Y)`$ |
| `NVASIQ` | Standard normal | $`(0,1)`$ | Fixed quantile | $`Q_Y(\mathtt{quantile})`$ |
| `LVASIQ` | Logistic | $`(0,1)`$ | Fixed quantile | $`Q_Y(\mathtt{quantile})`$ |
| `HVASIQ` | Hyperbolic secant | $`(0,1)`$ | Fixed quantile | $`Q_Y(\mathtt{quantile})`$ |
| `ZANVASIM` | Standard normal | $`[0,1)`$ | Zero-augmented mean | $`E(Y\mid Y>0)`$ |
| `OANVASIM` | Standard normal | $`(0,1]`$ | One-augmented mean | $`E(Y\mid Y<1)`$ |
| `ZOANVASIM` | Standard normal | $`[0,1]`$ | Zero-and-one-augmented mean | $`E(Y\mid 0<Y<1)`$ |

In every family, `sigma` is a shape parameter in $`(0,1)`$. The default
link is the logit for every distributional parameter. Thus, unless a
different link is requested, a coefficient is an effect on a logit
scale; it is not an additive change in a mean, quantile, or probability
on the response scale.

For `NVASIQ`, `LVASIQ`, and `HVASIQ`, the target probability is supplied
directly to the family constructor. For example,
`NVASIQ(quantile = 0.25)` defines a model in which `mu` is the
conditional first quartile. No global variable is required. The argument
`quantile` in these three families is unrelated to `tau` in `ZOANVASIM`,
where $`\tau=P(Y=1\mid Y>0)`$ is an estimated distributional parameter.

## Conditional and marginal interpretation

The distinction between the continuous-component mean and the marginal
mean is essential for the augmented models. If `mu`, `nu`, and `tau`
denote fitted values on the response scale, then the mixture
parameterizations used here follow the general treatment of
boundary-inflated models in Ospina and Ferrari
([2010](#ref-OspinaFerrari2010)) and Ospina and Ferrari
([2012](#ref-OspinaFerrari2012)).

| Family | Boundary probabilities | Marginal mean |
|:---|:---|:---|
| `ZANVASIM` | $`P(Y=0)=\nu`$ | $`(1-\nu)\mu`$ |
| `OANVASIM` | $`P(Y=1)=\nu`$ | $`\nu+(1-\nu)\mu`$ |
| `ZOANVASIM` | $`P(Y=0)=\nu`$, $`P(Y=1)=(1-\nu)\tau`$ | $`(1-\nu)\{\tau+(1-\tau)\mu\}`$ |

Consequently, a coefficient in the `mu` predictor of an augmented model
describes the conditional mean in the open interval. Its effect on the
marginal mean also depends on the fitted boundary probabilities.

## Helpers used in the examples

The following small functions keep the output compact. All reported
values are calculated when the vignette is built; no numerical results
are hard-coded.

``` r

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
```

## Overview of the included data

``` r

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
```

| data      | response     | observations | zero | interior | one |
|:----------|:-------------|-------------:|-----:|---------:|----:|
| bodyfat   | ARMS         |          298 |    0 |      298 |   0 |
| bodyfat   | LEGS         |          298 |    0 |      298 |   0 |
| bodyfat   | BODY         |          298 |    0 |      298 |   0 |
| bodyfat   | ANDROID      |          298 |    0 |      298 |   0 |
| bodyfat   | GYNECOID     |          298 |    0 |      298 |   0 |
| transport | propbiked    |           60 |   24 |       36 |   0 |
| trees     | prop         |           26 |    0 |       20 |   6 |
| aep       | noinap / los |         1383 |  763 |      552 |  68 |

Observed support of the bounded responses. {.table}

This empirical support determines the admissible families. The five
`bodyfat` responses lie strictly inside the interval and can be fitted
by any of the four base families. `transport`, `trees`, and `aep`
require boundary models because their responses contain zero, one, or
both, respectively. Replacing an exact boundary value by an arbitrary
small offset is unnecessary and would change the observed data.

## Body-fat proportions: mean and quantile regression

The `bodyfat` data contain five body-fat proportions measured on 298
individuals, together with age, body mass index (BMI), sex, and physical
activity. A value such as `ARMS = 0.163` represents 16.3%, so the
responses are already on the unit-interval scale and must not be divided
by 100. The data source and previous analyses are described by Petterle
et al. ([2020](#ref-PetterleEtAl2020)), Mazucheli et al.
([2021](#ref-MazucheliEtAl2021)), and Mazucheli et al.
([2022](#ref-MazucheliEtAl2022)).

The five responses are repeated measurements on the same individuals.
The models below analyze each response separately and therefore do not
estimate their cross-response dependence. A joint multivariate
interpretation would require a different model.

We center the continuous covariates to give the intercept a more useful
interpretation and explicitly label the factor levels.

``` r

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
```

| response | minimum | maximum |
|:---------|--------:|--------:|
| ARMS     |   0.042 |   0.547 |
| LEGS     |   0.068 |   0.574 |
| BODY     |   0.067 |   0.538 |
| ANDROID  |   0.072 |   0.580 |
| GYNECOID |   0.105 |   0.593 |

Ranges of the five body-fat proportions. {.table}

### Mean regression for all five responses

For each body-fat response, the same linear predictor is used for the
conditional mean and `sigma` is held constant. Keeping the formula
common facilitates comparisons of coefficient patterns across anatomical
regions, but AIC and BIC values from different responses should not be
interpreted as a competition among the responses.

``` r

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
```

| model    | family | parameters |  logLik |      AIC |      BIC | converged |
|:---------|:-------|-----------:|--------:|---------:|---------:|:----------|
| ARMS     | NVASIM |          7 | 455.611 | -897.221 | -871.342 | TRUE      |
| LEGS     | NVASIM |          7 | 424.419 | -834.838 | -808.958 | TRUE      |
| BODY     | NVASIM |          7 | 437.690 | -861.379 | -835.500 | TRUE      |
| ANDROID  | NVASIM |          7 | 375.644 | -737.288 | -711.408 | TRUE      |
| GYNECOID | NVASIM |          7 | 440.223 | -866.445 | -840.565 | TRUE      |

Normal-kernel Vasicek mean regressions for the body-fat responses.
{.table}

For these `NVASIM` fits, `fitted(object, what = "mu")` returns the
fitted conditional mean on the response scale. The coefficients
themselves are on the logit scale. The complete coefficient table for
any response can be obtained as follows.

``` r

knitr::kable(
  coefficient_table(fit_bodyfat_mean[["ARMS"]], c("mu", "sigma")),
  digits = 4,
  caption = "Coefficient estimates for the ARMS mean-regression model."
)
```

| parameter | term                      | estimate |
|:----------|:--------------------------|---------:|
| mu        | (Intercept)               |  -0.4684 |
| mu        | AGE_centered              |   0.0044 |
| mu        | BMI_centered              |   0.0877 |
| mu        | SEXmale                   |  -0.9049 |
| mu        | IPAQinsufficiently_active |  -0.1154 |
| mu        | IPAQactive                |  -0.2514 |
| sigma     | (Intercept)               |  -3.4739 |

Coefficient estimates for the ARMS mean-regression model. {.table}

### Comparing the three quantile kernels

We next fit median regressions to `ARMS` with identical predictors and
the normal, logistic, and hyperbolic-secant kernels. Here, `mu` is the
fitted conditional median, not the conditional mean. The same
construction can be used for another fixed level by changing
`quantile_level`.

``` r

quantile_level <- 0.50

fit_arms_nq <- gamlss(
  ARMS ~ AGE_centered + BMI_centered + SEX + IPAQ,
  sigma.formula = ~ 1,
  family = NVASIQ(quantile = quantile_level),
  data = bodyfat_analysis,
  control = control
)

fit_arms_lq <- gamlss(
  ARMS ~ AGE_centered + BMI_centered + SEX + IPAQ,
  sigma.formula = ~ 1,
  family = LVASIQ(quantile = quantile_level),
  data = bodyfat_analysis,
  control = control
)

fit_arms_hq <- gamlss(
  ARMS ~ AGE_centered + BMI_centered + SEX + IPAQ,
  sigma.formula = ~ 1,
  family = HVASIQ(quantile = quantile_level),
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
```

| model         | family | parameters |  logLik |      AIC |      BIC | converged |
|:--------------|:-------|-----------:|--------:|---------:|---------:|:----------|
| NVASIM_mean   | NVASIM |          7 | 455.611 | -897.221 | -871.342 | TRUE      |
| NVASIQ_median | NVASIQ |          7 | 455.495 | -896.989 | -871.109 | TRUE      |
| LVASIQ_median | LVASIQ |          7 | 451.201 | -888.402 | -862.522 | TRUE      |
| HVASIQ_median | HVASIQ |          7 | 446.193 | -878.387 | -852.507 | TRUE      |

Likelihood-based summaries for the ARMS models. {.table}

Because all four models use the same response observations, their
maximized likelihoods, AICs, and BICs can be compared. Such a comparison
concerns the complete fitted distributions. It does not make `mu`
directly comparable between the mean model and the median models.

``` r

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
```

| quantity | minimum | first_quartile | median | mean | third_quartile | maximum |
|:---|---:|---:|---:|---:|---:|---:|
| NVASIM_conditional_mean | 0.0942 | 0.1843 | 0.2502 | 0.2655 | 0.3406 | 0.5125 |
| NVASIQ_conditional_median | 0.0911 | 0.1807 | 0.2469 | 0.2627 | 0.3383 | 0.5128 |
| LVASIQ_conditional_median | 0.0913 | 0.1804 | 0.2457 | 0.2633 | 0.3398 | 0.5145 |
| HVASIQ_conditional_median | 0.0921 | 0.1807 | 0.2445 | 0.2634 | 0.3384 | 0.5147 |

Summaries of fitted means and medians for ARMS. {.table}

## Bicycle-trip proportions: zero augmentation

The `transport` data contain 60 respondents from a stratified
transportation study. The response is the proportion of trips to campus
made by bicycle. The data originate from the consulting study reported
by Korosteleva ([2019](#ref-Korosteleva2019)) and were subsequently
analyzed by Menezes et al. ([2021](#ref-MenezesEtAl2021)); the object
distributed by `vasicekreg` retains the values supplied by `uwquantreg`
([Menezes 2026](#ref-uwquantreg)).

Because `propbiked` contains exact zeros but no ones, `ZANVASIM` is the
appropriate augmented family. We use the covariate structure from the
previous analysis: gender, parking-permit duration, and institutional
status in the positive-component mean; gender and distance in the zero
probability. Centering the continuous predictors changes the intercepts
but not the fitted values or slopes.

``` r

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
```

| model    | family   | parameters |  logLik |    AIC |    BIC | converged |
|:---------|:---------|-----------:|--------:|-------:|-------:|:----------|
| ZANVASIM | ZANVASIM |          9 | -15.106 | 48.212 | 67.061 | TRUE      |

Likelihood-based summary for the transport model. {.table}

``` r


knitr::kable(
  coefficient_table(fit_transport, c("mu", "sigma", "nu")),
  digits = 4,
  caption = "Coefficient estimates for the zero-augmented transport model."
)
```

| parameter | term              | estimate |
|:----------|:------------------|---------:|
| mu        | (Intercept)       |   0.5379 |
| mu        | genderM           |   0.4728 |
| mu        | parking_centered  |  -0.1895 |
| mu        | statusstaff       |  -1.1928 |
| mu        | statusstudent     |  -1.1238 |
| sigma     | (Intercept)       |  -0.8408 |
| nu        | (Intercept)       |  -0.8846 |
| nu        | genderM           |   1.7728 |
| nu        | distance_centered |   0.3396 |

Coefficient estimates for the zero-augmented transport model. {.table}

Here, `mu` is the fitted mean bicycle-trip proportion among positive
responses and `nu` is the fitted probability of no bicycle trips. The
fitted marginal mean combines both components.

``` r

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
```

| quantity | minimum | first_quartile | median | mean | third_quartile | maximum |
|:---|---:|---:|---:|---:|---:|---:|
| positive_component_mean | 0.2290 | 0.4327 | 0.4980 | 0.5161 | 0.6142 | 0.7532 |
| probability_zero | 0.0217 | 0.1134 | 0.3364 | 0.4000 | 0.6776 | 1.0000 |
| marginal_mean | 0.0000 | 0.1580 | 0.3016 | 0.3014 | 0.4706 | 0.6213 |

Fitted quantities from the transport model. {.table}

Although `ntrips` and `nbiked` are available, this analysis treats each
respondent’s proportion as one mixed continuous–discrete response. It is
not a binomial model for `nbiked` conditional on `ntrips`, and
respondents with larger denominators do not automatically receive larger
weights.

## Tree-survival proportions: one augmentation

The `trees` data record two-year survival proportions in 26 parks. Their
provenance is the same as that of `transport` ([Korosteleva
2019](#ref-Korosteleva2019); [Menezes et al.
2021](#ref-MenezesEtAl2021); [Menezes 2026](#ref-uwquantreg)). The
response contains exact ones but no zeros, so we use `OANVASIM`.

The continuous-component mean depends on pest-control frequency,
fertilization frequency, precipitation, and wind speed. The probability
of complete survival depends on wind speed. Precipitation and wind are
centered to make the intercepts refer to their sample-average values.

``` r

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
```

| model    | family   | parameters | logLik |    AIC |    BIC | converged |
|:---------|:---------|-----------:|-------:|-------:|-------:|:----------|
| OANVASIM | OANVASIM |          8 |  0.539 | 14.923 | 24.987 | TRUE      |

Likelihood-based summary for the tree-survival model. {.table}

``` r


knitr::kable(
  coefficient_table(fit_trees, c("mu", "sigma", "nu")),
  digits = 4,
  caption = "Coefficient estimates for the one-augmented tree-survival model."
)
```

| parameter | term            | estimate |
|:----------|:----------------|---------:|
| mu        | (Intercept)     |  -1.2055 |
| mu        | pest            |   0.3280 |
| mu        | fertilization   |   0.9807 |
| mu        | precip_centered |  -0.0893 |
| mu        | wind_centered   |  -0.2141 |
| sigma     | (Intercept)     |  -1.2919 |
| nu        | (Intercept)     |  -1.5453 |
| nu        | wind_centered   |  -0.4671 |

Coefficient estimates for the one-augmented tree-survival model.
{.table}

For this model, `mu` is the fitted mean survival proportion conditional
on a value below one and `nu` is the probability of complete survival.
Therefore, the marginal mean is `nu + (1 - nu) * mu`.

``` r

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
```

| quantity | minimum | first_quartile | median | mean | third_quartile | maximum |
|:---|---:|---:|---:|---:|---:|---:|
| continuous_component_mean | 0.2553 | 0.5242 | 0.6711 | 0.6410 | 0.7872 | 0.9363 |
| probability_one | 0.0233 | 0.0684 | 0.2292 | 0.2308 | 0.3478 | 0.5928 |
| marginal_mean | 0.3348 | 0.5454 | 0.7607 | 0.7130 | 0.8782 | 0.9520 |

Fitted quantities from the tree-survival model. {.table}

The small sample of 26 parks warrants caution, particularly because
several distributional parameters are estimated. As in the transport
example, the model treats each park-level proportion as one response; it
is not a binomial model that uses `planted` as a number of trials.

## Inappropriate hospital-stay proportions: two-boundary augmentation

The `aep` data contain 1,383 patients admitted to Hospital del Mar in
Barcelona in 1988 and 1990. For each patient, `noinap` is the number of
days classified as inappropriate and `los` is the total length of stay.
The data were studied by Gange et al. ([1996](#ref-GangeEtAl1996)) and
the object in `vasicekreg` retains the structure supplied by
`gamlss.data` ([Stasinopoulos and Rigby 2025](#ref-gamlssdata)). We
define

``` math
Y_i=\mathrm{noinap}_i / \mathrm{los}_i.
```

The response includes both zero and one and is therefore fitted with
`ZOANVASIM`. In its sequential boundary parameterization,

``` math
P(Y_i=0)=\nu_i,\qquad
P(Y_i=1)=(1-\nu_i)\tau_i,\qquad
P(0<Y_i<1)=(1-\nu_i)(1-\tau_i).
```

Sex, ward, admission year, centered age, and length of stay enter the
continuous-component mean. Length of stay also enters the shape and both
boundary components. The variable `age` supplied in the data is already
age minus 55 years, and `loglos` is $`\log(\mathtt{los}/10)`$.

``` r

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
```

| model     | family    | parameters |   logLik |      AIC |      BIC | converged |
|:----------|:----------|-----------:|---------:|---------:|---------:|:----------|
| ZOANVASIM | ZOANVASIM |         13 | -841.493 | 1708.986 | 1777.002 | TRUE      |

Likelihood-based summary for the hospital-stay model. {.table}

``` r


knitr::kable(
  coefficient_table(fit_aep, c("mu", "sigma", "nu", "tau")),
  digits = 4,
  caption = "Coefficient estimates for the zero-and-one-augmented hospital-stay model."
)
```

| parameter | term        | estimate |
|:----------|:------------|---------:|
| mu        | (Intercept) |   0.0003 |
| mu        | sex2        |  -0.0025 |
| mu        | ward2       |  -0.4006 |
| mu        | ward3       |  -0.4206 |
| mu        | year90      |  -0.3826 |
| mu        | age         |   0.0041 |
| mu        | loglos      |   0.0525 |
| sigma     | (Intercept) |  -1.0861 |
| sigma     | loglos      |   0.6381 |
| nu        | (Intercept) |  -0.1827 |
| nu        | loglos      |  -1.1574 |
| tau       | (Intercept) |  -2.2294 |
| tau       | loglos      |  -1.1406 |

Coefficient estimates for the zero-and-one-augmented hospital-stay
model. {.table}

The four fitted components must be interpreted jointly. In particular,
`tau` is conditional on a nonzero response; it is not the marginal
probability of one. The following calculations recover the three
component probabilities and the marginal fitted mean.

``` r

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
```

| quantity | minimum | first_quartile | median | mean | third_quartile | maximum |
|:---|---:|---:|---:|---:|---:|---:|
| continuous_component_mean | 0.2613 | 0.3427 | 0.3991 | 0.4017 | 0.4510 | 0.5539 |
| shape | 0.0721 | 0.1354 | 0.2119 | 0.2255 | 0.2950 | 0.6050 |
| probability_zero | 0.0509 | 0.3607 | 0.5573 | 0.5517 | 0.7704 | 0.9229 |
| probability_one | 0.0068 | 0.0437 | 0.0557 | 0.0521 | 0.0633 | 0.0688 |
| probability_continuous | 0.0310 | 0.1611 | 0.3811 | 0.3962 | 0.5956 | 0.9423 |
| marginal_mean | 0.0542 | 0.1284 | 0.2182 | 0.2205 | 0.2953 | 0.5179 |

Fitted quantities from the hospital-stay model. {.table}

Gange et al. ([1996](#ref-GangeEtAl1996)) modeled the number of
inappropriate days conditional on length of stay using binomial and
beta-binomial models. The present analysis has a different sampling
formulation: the patient is the observational unit, and the
patient-level proportion is modeled by a distribution with two boundary
masses and a continuous interior component. Thus, the approaches should
not be described as the same likelihood with a different continuous
kernel.

## Reading coefficients and comparing models

For all fitted models, `summary(object)` supplies the usual GAMLSS
coefficient tables. A disciplined interpretation proceeds in three
stages:

1.  identify the modeled component (`mu`, `sigma`, `nu`, or `tau`);
2.  interpret its coefficient on the chosen link scale; and
3.  transform predictions to the response scale with
    [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and, for an
    augmented family, combine the components using the appropriate
    marginal mean formula.

Likelihood-based criteria compare complete fitted distributions only
when the response and observations are the same. A smaller AIC or BIC
does not by itself establish adequate residual behavior, and
coefficients belonging to different parameterizations should not be
equated merely because they share the name `mu`.

## Residual diagnostics and simulated envelopes

GAMLSS uses normalized randomized quantile residuals ([Dunn and Smyth
1996](#ref-DunnSmyth1996)). For an augmented distribution, the
probability integral transform is randomized over the fitted CDF jump at
the observed boundary. Repeated residual calculations can therefore
differ at observations equal to zero or one.

The package function
[`vasicek_envelope()`](https://jmazucheli.github.io/vasicekreg/reference/vasicek_envelope.md)
implements parametric-bootstrap pointwise envelopes for these residuals
and for generalized Cox–Snell residuals ([Cox and Snell
1968](#ref-CoxSnell1968)). Each accepted bootstrap sample is simulated
from the fitted model, the model is re-estimated, and the ordered
residuals are recalculated. The construction is related to
simulated-envelope diagnostics described by Atkinson
([1985](#ref-Atkinson1985)), Moral et al. ([2017](#ref-MoralEtAl2017)),
and Zhao et al. ([2011](#ref-ZhaoEtAl2011)).

Because a publication-quality envelope requires hundreds of model
re-estimations, the code is shown but not executed while the vignette is
built. It can be applied to any of the fitted objects above.

``` r

envelope_aep <- vasicek_envelope(
  object = fit_aep,
  residual = c("quantile", "cox-snell"),
  nsim = 500,
  level = 0.95,
  envelope = "quantile",
  seed = 2026,
  data = aep_analysis
)

old_par <- graphics::par(no.readonly = TRUE)
graphics::par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
plot(envelope_aep, which = "quantile", pch = 19, cex = 0.55)
plot(envelope_aep, which = "cox-snell", pch = 19, cex = 0.55)
graphics::par(old_par)
```

These are full quantile–quantile plots. The gray region is a pointwise,
not simultaneous, envelope; the red identity line is the theoretical
reference and the blue curve is the pointwise mean of the ordered
bootstrap residuals. In finite samples, particularly in the upper tail
of the Cox–Snell plot, the simulated mean can be the more informative
reference.

## Practical workflow

For a new bounded response, the following sequence is recommended:

1.  verify whether the observed support is $`(0,1)`$, $`[0,1)`$,
    $`(0,1]`$, or $`[0,1]`$;
2.  choose a compatible family rather than transforming exact
    boundaries;
3.  state whether `mu` is a mean or a fixed quantile;
4.  specify predictors separately for every scientifically relevant
    distributional component;
5.  inspect convergence and coefficient estimates;
6.  calculate fitted quantities on the response scale, including the
    marginal mean for augmented models; and
7.  assess residual behavior, preferably with simulated envelopes when
    the final model is selected.

The examples in this vignette are reproducible templates, not automatic
model-selection prescriptions. Covariate structures should ultimately
follow the scientific question, the sampling design, and the information
available in each data set.

## Session information

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] parallel  splines   stats     graphics  grDevices utils     datasets  methods  
#> [9] base     
#> 
#> other attached packages:
#> [1] vasicekreg_1.2.0  gamlss_5.5-0      nlme_3.1-169      gamlss.dist_6.1-1
#> [5] gamlss.data_6.0-7
#> 
#> loaded via a namespace (and not attached):
#>  [1] cli_3.6.6         knitr_1.51        rlang_1.3.0       xfun_0.60        
#>  [5] otel_0.2.0        textshaping_1.0.5 jsonlite_2.0.0    htmltools_0.5.9  
#>  [9] ragg_1.5.2        sass_0.4.10       rmarkdown_2.32    grid_4.6.1       
#> [13] evaluate_1.0.5    jquerylib_0.1.4   MASS_7.3-65       fastmap_1.2.0    
#> [17] mvtnorm_1.4-2     yaml_2.3.12       lifecycle_1.0.5   compiler_4.6.1   
#> [21] fs_2.1.0          Rcpp_1.1.2        lattice_0.22-9    systemfonts_1.3.2
#> [25] digest_0.6.39     R6_2.6.1          Matrix_1.7-5      bslib_0.12.0     
#> [29] tools_4.6.1       survival_3.8-6    pkgdown_2.2.1     cachem_1.1.0     
#> [33] desc_1.4.3
```

## References

Atkinson, Anthony C. 1985. *Plots, Transformations and Regression: An
Introduction to Graphical Methods of Diagnostic Regression Analysis*.
Oxford University Press.

Cox, David R., and E. J. Snell. 1968. “A General Definition of
Residuals.” *Journal of the Royal Statistical Society: Series B
(Methodological)* 30 (2): 248–75.
<https://doi.org/10.1111/j.2517-6161.1968.tb00724.x>.

Dunn, Peter K., and Gordon K. Smyth. 1996. “Randomized Quantile
Residuals.” *Journal of Computational and Graphical Statistics* 5 (3):
236–44. <https://doi.org/10.1080/10618600.1996.10474708>.

Fischer, Matthias J., Alexander Hui, and Stefan Hösle. 2017. “wHS-Type
Distributions with Application to Finance.” *Journal of Statistics and
Management Systems* 20 (1): 67–89.
<https://doi.org/10.1080/09720510.2016.1190575>.

Gange, Stephen J., Alvaro Munoz, Marc Saez, and Jordi Alonso. 1996. “Use
of the Beta-Binomial Distribution to Model the Effect of Policy Changes
on Appropriateness of Hospital Stays.” *Applied Statistics* 45 (3):
371–82.

Korosteleva, Olga. 2019. *Advanced Regression Models with SAS and R*.
CRC Press.

Mazucheli, Josmar, Bruna Alves, Mehmet Ç. Korkmaz, and Víctor Leiva.
2022. “Vasicek Quantile and Mean Regression Models for Bounded Data: New
Formulation, Mathematical Derivations, and Numerical Applications.”
*Mathematics* 10 (9): 1389. <https://doi.org/10.3390/math10091389>.

Mazucheli, Josmar, Víctor Leiva, Bruna Alves, and André F. B. Menezes.
2021. “A New Quantile Regression for Modeling Bounded Data Under a Unit
Birnbaum–Saunders Distribution with Applications in Medicine and
Politics.” *Symmetry* 13 (4): 682.
<https://doi.org/10.3390/sym13040682>.

Menezes, André F. B. 2026. *uwquantreg: Unit-Weibull Quantile
Regression*. <https://github.com/AndrMenezes/uwquantreg>.

Menezes, André F. B., Josmar Mazucheli, and Marcelo Bourguignon. 2021.
“A Parametric Quantile Regression Approach for Modeling Zero- or
One-Inflated Double Bounded Data.” *Biometrical Journal* 63 (4): 841–58.
<https://doi.org/10.1002/bimj.202000126>.

Moral, Rafael A., John Hinde, and Clarice G. B. Demétrio. 2017.
“Half-Normal Plots and Overdispersed Models in R: The hnp Package.”
*Journal of Statistical Software* 81 (10): 1–23.
<https://doi.org/10.18637/jss.v081.i10>.

Ospina, Raydonal, and Silvia L. P. Ferrari. 2010. “Inflated Beta
Distributions.” *Statistical Papers* 51 (1): 111–26.
<https://doi.org/10.1007/s00362-008-0125-4>.

Ospina, Raydonal, and Silvia L. P. Ferrari. 2012. “A General Class of
Zero-or-One Inflated Beta Regression Models.” *Computational Statistics
& Data Analysis* 56 (6): 1609–23.
<https://doi.org/10.1016/j.csda.2011.10.005>.

Petterle, Ricardo R., Wagner H. Bonat, Cassius T. Scarpin, Thaísa
Jonasson, and Victória Z. C. Borba. 2020. “Multivariate Quasi-Beta
Regression Models for Continuous Bounded Data.” *The International
Journal of Biostatistics* 17 (1): 39–53.
<https://doi.org/10.1515/ijb-2019-0163>.

Rigby, Robert A., and D. Mikis Stasinopoulos. 2005. “Generalized
Additive Models for Location, Scale and Shape.” *Journal of the Royal
Statistical Society: Series C (Applied Statistics)* 54 (3): 507–54.
<https://doi.org/10.1111/j.1467-9876.2005.00510.x>.

Stasinopoulos, D. Mikis, and Robert A. Rigby. 2007. “Generalized
Additive Models for Location Scale and Shape (GAMLSS) in R.” *Journal of
Statistical Software* 23 (7): 1–46.
<https://doi.org/10.18637/jss.v023.i07>.

Stasinopoulos, Mikis, and Robert Rigby. 2025. *gamlss.data: Data for
Generalized Additive Models for Location Scale and Shape*.
<https://CRAN.R-project.org/package=gamlss.data>.

Vasicek, Oldrich A. 2002. “The Distribution of Loan Portfolio Value.”
*Risk* 15 (12): 160–62.

Witzany, Jiří. 2013. “A Note on the Vasicek’s Model with the Logistic
Distribution.” *Ekonomický Časopis (Journal of Economics)* 61 (10):
1053–66. <https://ideas.repec.org/p/fau/wpaper/wp2013_01.html>.

Zhao, Yan, Andy H. Lee, Kelvin K. W. Yau, and Geoffrey J. McLachlan.
2011. “Assessing the Adequacy of Weibull Survival Models: A Simulated
Envelope Approach.” *Journal of Applied Statistics* 38 (10): 2089–97.
<https://doi.org/10.1080/02664763.2010.545115>.
