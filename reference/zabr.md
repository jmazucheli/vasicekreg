# Zero-Augmented Beta random-intercept model

Fits a two-component model for responses in \\\[0,1)\\ with a point mass
at zero (discrete component, parameterized by \\\gamma\\) and a Beta
distribution for positive values (continuous component, parameterized by
\\\beta\\).

## Usage

``` r
zabr(
  data,
  y,
  formula_bin = NULL,
  formula_cont = NULL,
  random = NULL,
  logistic_cov = NULL,
  beta_cov = NULL,
  subject_ind = NULL,
  time_ind,
  component_wise_test = TRUE,
  quad_n = 30,
  verbose = FALSE,
  joint_test = NULL,
  sd_lower = 1e-05,
  start = NULL,
  control = list(),
  hessian = TRUE
)
```

## Arguments

- data:

  A data.frame.

- y:

  Character name of the response column.

- formula_bin:

  One-sided formula for the presence component.

- formula_cont:

  One-sided formula for the Beta component.

- random:

  Formula `~ 1 | Subject` for the random intercept.

- logistic_cov:

  Deprecated: character vector of presence covariates.

- beta_cov:

  Deprecated: character vector of Beta covariates.

- subject_ind:

  Deprecated: subject column name.

- time_ind:

  Character name of the time column.

- component_wise_test:

  Logical; compute component-wise LRTs?

- quad_n:

  Integer; number of Gauss-Hermite quadrature points.

- verbose:

  Logical; print progress?

- joint_test:

  NULL, TRUE, or FALSE; compute joint LRTs?

- sd_lower:

  Lower bound for the random-effect SD.

- start:

  Optional list of starting values.

- control:

  List of control parameters.

- hessian:

  Logical; compute Hessian-based standard errors?

## Value

An object of class `"zabr"` with components including
`logistic_est_table`, `beta_est_table`, `precision_table`,
`random_effects`, `loglikelihood`, `joint_p`, `fit_statistics` (joint
model), `fit_statistics_bin` and `fit_statistics_cont` (per component;
see section *Fit statistics*), and `vcov`. Fixed-effect names in
[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) use the prefix `gamma_`
for the discrete component and `beta_` for the continuous component,
including their intercepts.

## Details

Formula and legacy arguments are mutually exclusive for each model
component: use either `formula_bin` or `logistic_cov`, either
`formula_cont` or `beta_cov`, and either `random` or `subject_ind`.
Supplying both arguments in any pair is an error.

For a response \\Y\_{it} \in \[0,1)\\ observed on subject \\i\\ (\\i =
1, \ldots, N\\) at time \\t\\ (\\t = 1, \ldots, T\\), the model places a
point mass at zero and a continuous component on \\(0,1)\\:

\$\$ Y\_{it} = 0 \quad \mbox{with probability } 1 - p\_{it}, \$\$ \$\$
Y\_{it} \sim \mathrm{Beta}\left(\mu\_{it}\phi,\\
(1-\mu\_{it})\phi\right) \quad \mbox{with probability } p\_{it}, \$\$

where \\0 \< p\_{it} \< 1\\, \\0 \< \mu\_{it} \< 1\\, and \\\phi \> 0\\.
Let \\X\_{it}\\ and \\Z\_{it}\\ be the covariate vectors entering the
discrete and continuous components, respectively; they may share columns
or be disjoint. Both components are modeled on the logit scale:

\$\$ \mathrm{logit}(p\_{it}) = \log\left(\frac{p\_{it}}{1 -
p\_{it}}\right) = a_i + \gamma_0 + X\_{it}^\top \gamma, \$\$ \$\$
\mathrm{logit}(\mu\_{it}) = \log\left(\frac{\mu\_{it}}{1 -
\mu\_{it}}\right) = b_i + \beta_0 + Z\_{it}^\top \beta, \$\$

where \\a_i\\ and \\b_i\\ are subject-specific random intercepts that
induce correlation across repeated measurements on the same subject,

\$\$ a_i \sim N(0, \sigma_1^2), \qquad b_i \sim N(0, \sigma_2^2). \$\$

The fixed-effect coefficients are denoted by \\\gamma\\ (discrete
component) and \\\beta\\ (continuous component); \\\phi\\ is the Beta
precision parameter. The two random intercepts are independent, and the
components have no shared parameters. Consequently, the marginal
likelihood factorizes into discrete and continuous components. The two
components are optimized separately; this is equivalent to maximizing
their joint likelihood. Random effects are integrated out using
non-adaptive Gauss–Hermite quadrature.

## Fit statistics

`fit_statistics` refers to the joint model and follows the "Fit
Statistics" table of SAS PROC NLMIXED: \\AIC = -2\ell + 2k\\ and \\BIC =
-2\ell + k\log(s)\\, where \\\ell\\ is the maximized log-likelihood,
\\k\\ the total number of parameters and \\s\\ the number of subjects.
BIC uses subjects, not observations.
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) return the same values as
`fit_statistics`.

`fit_statistics_bin` and `fit_statistics_cont` are informational. The
continuous component counts only observations with \\Y \> 0\\ and
subjects with at least one such observation, as if it were fitted alone
to the positive responses. Consequently, AIC is additive across
components, whereas BIC is additive only when every subject has at least
one positive response. To compare models, use `fit_statistics` or
[`BIC()`](https://rdrr.io/r/stats/AIC.html); do not sum component-wise
values.

[`logLik()`](https://rdrr.io/r/stats/logLik.html) sets the `nobs`
attribute to the number of subjects, so that
[`BIC()`](https://rdrr.io/r/stats/AIC.html) follows PROC NLMIXED.

The values coincide numerically with PROC NLMIXED only when the
log-likelihood coincides: the same model with independent random
intercepts and the same number of parameters, non-adaptive quadrature
with the same number of points (`NOAD` and `QPOINTS=` equal to
`quad_n`), and data sorted by subject.

## References

Chen, E. Z. and Li, H. (2016). A two-part mixed-effects model for
analyzing longitudinal microbiome compositional data. *Bioinformatics*,
**32**(17), 2611–2617.
[doi:10.1093/bioinformatics/btw308](https://doi.org/10.1093/bioinformatics/btw308)

## See also

[`zavr`](https://jmazucheli.github.io/vasicekreg/reference/zavr.md)

## Examples

``` r
# \donttest{
data(please_microbiome)

d <- subset(please_microbiome, Genus == "g__Bifidobacterium")
d$Subject <- factor(d$Subject)

fit <- zabr(
  data         = d,
  y            = "Y",
  formula_bin  = ~ Baseline + Time + Treat,
  formula_cont = ~ Baseline + Time + Treat,
  random       = ~ 1 | Subject,
  time_ind     = "Time"
)
fit
#> Zero-augmented Beta random-intercept model
#> Discrete component (gamma; Pvalue = LRT; Wald_Pvalue = normal reference):
#>               Estimate         SE      Pvalue Wald_Pvalue
#> (Intercept)  2.2183470  0.7988457 0.000930524  0.00548732
#> Baseline    27.9528419 17.7653473 0.011359847  0.11561523
#> Time         0.0569866  0.0932645 0.538603209  0.54118563
#> TreatEEN    -2.7560317  1.0486498 0.003454401  0.00858452
#> 
#> Continuous component (beta):
#>               Estimate        SE      Pvalue Wald_Pvalue
#> (Intercept) -2.8468204 0.2136013 3.99538e-27 1.59677e-40
#> Baseline     2.3064359 0.6424036 6.71193e-04 3.30269e-04
#> Time         0.0287498 0.0282817 3.11436e-01 3.09367e-01
#> TreatEEN    -0.7491682 0.3608348 3.35286e-02 3.78747e-02
#> 
#> Beta precision:
#>     Estimate      SE
#> phi  8.45022 1.50465
#> 
#> Random effects:
#>                   Estimate       SE
#> Presence_SD       2.011728 0.661405
#> Presence_variance 4.047049 2.661134
#> Positive_SD       0.612680 0.118277
#> Positive_variance 0.375377 0.144932
#> 
#> Joint likelihood-ratio p-values (2 df):
#>   Baseline       Time   TreatEEN 
#> 0.00012493 0.49591229 0.00145244 
#> 
#> Fit statistics, joint model (BIC uses subjects, as in PROC NLMIXED):
#>  N_obs N_subjects  K Neg2LogLik      AIC      BIC
#>    177         59 11   -428.969 -406.969 -384.116
#> 
#> Component-wise fit statistics (informational; AIC is additive,
#> BIC is additive only if every subject has at least one Y > 0):
#>  Component   N N_subjects K Neg2LogLik      AIC      BIC
#>   presence 177         59 5    126.583  136.583  146.971
#>   positive 142         53 6   -555.552 -543.552 -531.730
#> 
#> Full-model diagnostics:
#>            Fit Convergence Boundary    LogLik Hessian                  Message
#>  presence_full           0    FALSE -63.29143      OK relative convergence (4)
#>      beta_full           0    FALSE 277.77609      OK relative convergence (4)
# }
```
