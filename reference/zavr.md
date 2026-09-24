# Zero-Augmented Vasicek random-intercept model

Fits a two-component model for responses in \\\[0,1)\\ with a point mass
at zero (discrete component, parameterized by \\\gamma\\) and a Vasicek
distribution for positive values (continuous component, parameterized by
\\\beta\\).

## Usage

``` r
zavr(
  data,
  y,
  formula_bin = NULL,
  formula_cont = NULL,
  random = NULL,
  logistic_cov = NULL,
  vasicek_cov = NULL,
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

  One-sided formula for the Vasicek component.

- random:

  Formula `~ 1 | Subject` for the random intercept.

- logistic_cov:

  Deprecated: character vector of presence covariates.

- vasicek_cov:

  Deprecated: character vector of Vasicek covariates.

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

An object of class `"zavr"` with components including
`logistic_est_table`, `vasicek_est_table`, `shape_table`,
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
`formula_cont` or `vasicek_cov`, and either `random` or `subject_ind`.
Supplying both arguments in any pair is an error.

For a response \\Y\_{it} \in \[0,1)\\ observed on subject \\i\\ (\\i =
1, \ldots, N\\) at time \\t\\ (\\t = 1, \ldots, T\\), the model places a
point mass at zero and a continuous component on \\(0,1)\\:

\$\$ Y\_{it} = 0 \quad \mbox{with probability } 1 - p\_{it}, \$\$ \$\$
Y\_{it} \sim \mathrm{NVASIM}\left(\mu\_{it}, \sigma\right) \quad
\mbox{with probability } p\_{it}, \$\$

where \\0 \< p\_{it} \< 1\\, \\0 \< \mu\_{it} \< 1\\, and \\\sigma \in
(0,1)\\. The continuous component is the normal-kernel Vasicek
(mean-parameterized) distribution `NVASIM`, for which \\\mu\_{it} =
E(Y\_{it} \mid Y\_{it} \> 0)\\. Let \\X\_{it}\\ and \\Z\_{it}\\ be the
covariate vectors entering the discrete and continuous components,
respectively; they may share columns or be disjoint. Both components are
modeled on the logit scale:

\$\$ \mathrm{logit}(p\_{it}) = \log\left(\frac{p\_{it}}{1 -
p\_{it}}\right) = a_i + \gamma_0 + X\_{it}^\top \gamma, \$\$ \$\$
\mathrm{logit}(\mu\_{it}) = \log\left(\frac{\mu\_{it}}{1 -
\mu\_{it}}\right) = b_i + \beta_0 + Z\_{it}^\top \beta, \$\$

where \\a_i\\ and \\b_i\\ are subject-specific random intercepts that
induce correlation across repeated measurements on the same subject,

\$\$ a_i \sim N(0, \sigma_1^2), \qquad b_i \sim N(0, \sigma_2^2). \$\$

The fixed-effect coefficients are denoted by \\\gamma\\ (discrete
component) and \\\beta\\ (continuous component); \\\sigma\\ is the shape
parameter of the Vasicek component. The two random intercepts are
independent, and the components have no shared parameters. Consequently,
the marginal likelihood factorizes into discrete and continuous
components. The two components are optimized separately; this is
equivalent to maximizing their joint likelihood. Random effects are
integrated out using non-adaptive Gauss–Hermite quadrature.

This function differs from
[`zabr`](https://jmazucheli.github.io/vasicekreg/reference/zabr.md) only
in the continuous component: `zabr` uses a Beta distribution with
precision \\\phi\\, whereas `zavr` uses the Vasicek `NVASIM`
distribution with shape \\\sigma\\. The two models share the same
discrete component and the same random-intercept structure.

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

Mazucheli, J., Alves, B., Korkmaz, M. Ç., and Leiva, V. (2022). Vasicek
quantile and mean regression models for bounded data: New formulation,
mathematical derivations, and numerical applications. *Mathematics*,
**10**, 1389.
[doi:10.3390/math10091389](https://doi.org/10.3390/math10091389)

## See also

[`zabr`](https://jmazucheli.github.io/vasicekreg/reference/zabr.md),
[`NVASIM`](https://jmazucheli.github.io/vasicekreg/reference/NVASIM.md)

## Examples

``` r
# \donttest{
data(please_microbiome)

d <- subset(please_microbiome, Genus == "g__Bifidobacterium")
d$Subject <- factor(d$Subject)

fit <- zavr(
  data         = d,
  y            = "Y",
  formula_bin  = ~ Baseline + Time + Treat,
  formula_cont = ~ Baseline + Time + Treat,
  random       = ~ 1 | Subject,
  time_ind     = "Time"
)
fit
#> Zero-augmented Vasicek random-intercept model
#> Discrete component (gamma; Pvalue = LRT; Wald_Pvalue = normal reference):
#>               Estimate         SE      Pvalue Wald_Pvalue
#> (Intercept)  2.2183470  0.7988457 0.000930524  0.00548732
#> Baseline    27.9528419 17.7653473 0.011359847  0.11561523
#> Time         0.0569866  0.0932645 0.538603209  0.54118563
#> TreatEEN    -2.7560317  1.0486498 0.003454401  0.00858452
#> 
#> Continuous component (beta):
#>               Estimate        SE      Pvalue Wald_Pvalue
#> (Intercept) -3.3079885 0.2926297 2.46941e-20 1.24875e-29
#> Baseline     3.2474703 1.0361603 2.64694e-03 1.72359e-03
#> Time         0.0424343 0.0354811 2.33831e-01 2.31708e-01
#> TreatEEN    -1.3035290 0.5684559 2.21549e-02 2.18420e-02
#> 
#> Vasicek shape:
#>       Estimate        SE
#> shape 0.310053 0.0326822
#> 
#> Random effects:
#>                   Estimate       SE
#> Presence_SD        2.01173 0.661405
#> Presence_variance  4.04705 2.661134
#> Positive_SD        1.07124 0.187361
#> Positive_variance  1.14756 0.401417
#> 
#> Joint likelihood-ratio p-values (2 df):
#>    Baseline        Time    TreatEEN 
#> 0.000442912 0.407477250 0.001015875 
#> 
#> Fit statistics, joint model (BIC uses subjects, as in PROC NLMIXED):
#>  N_obs N_subjects  K Neg2LogLik      AIC      BIC
#>    177         59 11   -451.842 -429.842 -406.989
#> 
#> Component-wise fit statistics (informational; AIC is additive,
#> BIC is additive only if every subject has at least one Y > 0):
#>  Component   N N_subjects K Neg2LogLik      AIC      BIC
#>   presence 177         59 5    126.583  136.583  146.971
#>   positive 142         53 6   -578.425 -566.425 -554.603
#> 
#> Full-model diagnostics:
#>            Fit Convergence Boundary    LogLik Hessian                  Message
#>  presence_full           0    FALSE -63.29143      OK relative convergence (4)
#>   vasicek_full           0    FALSE 289.21239      OK relative convergence (4)
# }
```
