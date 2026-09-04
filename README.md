# vasicekreg

[![CRAN status](https://www.r-pkg.org/badges/version/vasicekreg)](https://CRAN.R-project.org/package=vasicekreg)

`vasicekreg` provides distribution functions and GAMLSS regression families
for Vasicek-type distributions on the unit interval. The package supports mean
and quantile regression under the standard normal kernel and quantile
regression under the logistic and hyperbolic-secant kernels. Normal-kernel families augmented at
zero, at one, or at both boundaries are available for responses containing
exact boundary values.

## Available families

| Family | Kernel | Response range | Parameterization | Interpretation of `mu` |
|---|---|---|---|---|
| `NVASIM` | Standard normal | `(0, 1)` | Mean | `mu = E(Y)` |
| `NVASIQ` | Standard normal | `(0, 1)` | Quantile | `mu = Q_Y(quantile)` |
| `LVASIQ` | Logistic | `(0, 1)` | Quantile | `mu = Q_Y(quantile)` |
| `HVASIQ` | Hyperbolic secant | `(0, 1)` | Quantile | `mu = Q_Y(quantile)` |
| `ZANVASIM` | Standard normal | `[0, 1)` | Zero-augmented mean | `mu = E(Y | Y > 0)` |
| `OANVASIM` | Standard normal | `(0, 1]` | One-augmented mean | `mu = E(Y | Y < 1)` |
| `ZOANVASIM` | Standard normal | `[0, 1]` | Zero-and-one-augmented | `mu = E(Y | 0 < Y < 1)` |

The documentation uses *augmented* for boundary mixtures and *Vasicek-type*
for the kernel-based constructions. The established exported names are kept
for backward compatibility and retain the familiar GAMLSS abbreviations `ZA`
and `OA`, historically read as zero- and one-adjusted.

For all seven families, `sigma` is a shape parameter in `(0, 1)` that
controls dispersion.
For `ZANVASIM`, `nu = P(Y = 0)` and the marginal mean is
`E(Y) = (1 - nu) * mu`.
For `OANVASIM`, `nu = P(Y = 1)` and the marginal mean is
`E(Y) = nu + (1 - nu) * mu`. In `ZOANVASIM`,
`nu = P(Y = 0)` and `tau = P(Y = 1 | Y > 0)`.
Its marginal mean is
`E(Y) = (1 - nu) * (tau + (1 - tau) * mu)`.
The logistic- and hyperbolic-secant-kernel means do not have closed-form
expressions and generally differ from `mu`; consequently, the package does not
provide mean-regression families for these kernels.

Each parameterization includes density, cumulative distribution, quantile,
and random generation functions:

- `dNVASIM()`, `pNVASIM()`, `qNVASIM()`, and `rNVASIM()`;
- `dNVASIQ()`, `pNVASIQ()`, `qNVASIQ()`, and `rNVASIQ()`;
- `dLVASIQ()`, `pLVASIQ()`, `qLVASIQ()`, and `rLVASIQ()`;
- `dHVASIQ()`, `pHVASIQ()`, `qHVASIQ()`, and `rHVASIQ()`;
- `d0NVASIM()`, `p0NVASIM()`, `q0NVASIM()`, and `r0NVASIM()`;
- `d1NVASIM()`, `p1NVASIM()`, `q1NVASIM()`, and `r1NVASIM()`;
- `d01NVASIM()`, `p01NVASIM()`, `q01NVASIM()`, and `r01NVASIM()`.

## Installation

Install the released version from CRAN:

```r
install.packages("vasicekreg")
```

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("jmazucheli/vasicekreg")
```

Because the package contains compiled C++ code, installation from source
requires an appropriate compiler toolchain.

## Included datasets

| Object | Main bounded response | Boundary pattern | Observations |
|---|---|---|---:|
| `bodyfat` | `ARMS`, `LEGS`, `BODY`, `ANDROID`, and `GYNECOID` | All in `(0, 1)` | 298 |
| `aep` | `noinap / los` | Zero and one | 1,383 |
| `transport` | `propbiked` | Zero | 60 |
| `trees` | `prop` | One | 26 |

The help pages preserve the original data provenance and the citations
associated with each dataset. Use `?bodyfat`, `?aep`, `?transport`, or
`?trees` for variable definitions, sources, and references.

## Distribution functions

The following example uses the normal-kernel mean parameterization:

```r
library(vasicekreg)

set.seed(123)
y <- rNVASIM(n = 1000, mu = 0.50, sigma = 0.25)

dNVASIM(x = 0.50, mu = 0.50, sigma = 0.25)
pNVASIM(q = 0.50, mu = 0.50, sigma = 0.25)
qNVASIM(p = 0.50, mu = 0.50, sigma = 0.25)
```

For the quantile parameterizations, the fixed level is supplied through the
`quantile` argument:

```r
qNVASIQ(p = 0.25, mu = 0.60, sigma = 0.25, quantile = 0.25)
qLVASIQ(p = 0.25, mu = 0.60, sigma = 0.25, quantile = 0.25)
qHVASIQ(p = 0.25, mu = 0.60, sigma = 0.25, quantile = 0.25)
```

All three calls return `mu = 0.60` because `mu` represents the quantile at
the fixed level `quantile = 0.25`.

## GAMLSS mean regression

```r
library(gamlss)
library(vasicekreg)

set.seed(123)
dat_mean <- data.frame(
  y = rNVASIM(n = 300, mu = 0.60, sigma = 0.25)
)

fit_mean <- gamlss(
  y ~ 1,
  data = dat_mean,
  family = NVASIM(
    mu.link = "logit",
    sigma.link = "logit"
  ),
  control = gamlss.control(trace = FALSE)
)

fitted(fit_mean, what = "mu")[1]
```

## GAMLSS zero-augmented mean regression

```r
set.seed(123)
dat_zero <- data.frame(
  y = r0NVASIM(n = 500, mu = 0.60, sigma = 0.25, nu = 0.20)
)

fit_zero <- gamlss(
  y ~ 1,
  sigma.formula = ~ 1,
  nu.formula = ~ 1,
  data = dat_zero,
  family = ZANVASIM(
    mu.link = "logit",
    sigma.link = "logit",
    nu.link = "logit"
  ),
  control = gamlss.control(trace = FALSE)
)

fitted(fit_zero, what = "mu")[1]
fitted(fit_zero, what = "sigma")[1]
fitted(fit_zero, what = "nu")[1]
```

Here `mu` is the conditional mean of the positive component. The fitted
marginal mean is `(1 - fitted(fit_zero, what = "nu")) *
fitted(fit_zero, what = "mu")`.

## GAMLSS one-augmented and zero-and-one-augmented regression

```r
dat_one <- data.frame(
  y = r1NVASIM(500, mu = 0.60, sigma = 0.25, nu = 0.20)
)

fit_one <- gamlss(
  y ~ 1,
  sigma.formula = ~ 1,
  nu.formula = ~ 1,
  data = dat_one,
  family = OANVASIM(),
  control = gamlss.control(trace = FALSE)
)

dat_boundary <- data.frame(
  y = r01NVASIM(
    500, mu = 0.60, sigma = 0.25, nu = 0.20, tau = 0.25
  )
)

fit_boundary <- gamlss(
  y ~ 1,
  sigma.formula = ~ 1,
  nu.formula = ~ 1,
  tau.formula = ~ 1,
  data = dat_boundary,
  family = ZOANVASIM(),
  control = gamlss.control(trace = FALSE)
)
```

For `ZOANVASIM`, `p0 = nu`, `p1 = (1 - nu) * tau`, and
`pc = (1 - nu) * (1 - tau)`. This parameterization guarantees valid
probabilities while retaining logit links for both boundary parameters. Here
`tau` is an estimated model parameter and is unrelated to the fixed
`quantile` argument of `NVASIQ()`, `LVASIQ()`, and `HVASIQ()`.

The included `aep` data provide a real example with both boundary values:

```r
data("aep", package = "vasicekreg")
aep$inappropriate <- with(aep, noinap / los)

fit_aep <- gamlss(
  inappropriate ~ sex + ward + year + age + loglos,
  sigma.formula = ~ loglos,
  nu.formula = ~ loglos,
  tau.formula = ~ loglos,
  family = ZOANVASIM(),
  data = aep,
  control = gamlss.control(n.cyc = 200, trace = FALSE)
)
```

## GAMLSS quantile regression

For `NVASIQ()`, `LVASIQ()`, and `HVASIQ()`, the fixed quantile level is passed
directly to the family constructor through `quantile`. Each fitted object is
self-contained, so models at different levels can coexist without changing
global variables.

### Normal kernel

```r
library(gamlss)
library(vasicekreg)

set.seed(123)
quantile_level <- 0.50

dat_normal <- data.frame(
  y = rNVASIQ(
    n = 300,
    mu = 0.60,
    sigma = 0.25,
    quantile = quantile_level
  )
)

fit_normal <- gamlss(
  y ~ 1,
  data = dat_normal,
  family = NVASIQ(
    quantile = quantile_level,
    mu.link = "logit",
    sigma.link = "logit"
  ),
  control = gamlss.control(trace = FALSE)
)

fitted(fit_normal, what = "mu")[1]
```

### Logistic kernel

```r
set.seed(123)
quantile_level <- 0.25

dat_logistic <- data.frame(
  y = rLVASIQ(
    n = 300,
    mu = 0.60,
    sigma = 0.25,
    quantile = quantile_level
  )
)

fit_logistic <- gamlss(
  y ~ 1,
  data = dat_logistic,
  family = LVASIQ(
    quantile = quantile_level,
    mu.link = "logit",
    sigma.link = "logit"
  ),
  control = gamlss.control(trace = FALSE)
)

fitted(fit_logistic, what = "mu")[1]
```

### Hyperbolic-secant kernel

```r
set.seed(123)
quantile_level <- 0.25

dat_hs <- data.frame(
  y = rHVASIQ(
    n = 300,
    mu = 0.60,
    sigma = 0.25,
    quantile = quantile_level
  )
)

fit_hs <- gamlss(
  y ~ 1,
  sigma.formula = ~ 1,
  data = dat_hs,
  family = HVASIQ(
    quantile = quantile_level,
    mu.link = "logit",
    sigma.link = "logit"
  ),
  control = gamlss.control(trace = FALSE)
)

fitted(fit_hs, what = "mu")[1]
```

The fixed level is embedded in the family and residual definitions.
`vasicek_envelope()` also inserts that stored value into every bootstrap refit,
so it does not depend on a global variable or on a symbol used in the original
family call. The generic `update()` method retains the usual R call semantics;
when direct later use of `update()` is planned, use a literal level in the
original call (for example, `NVASIQ(quantile = 0.25)`) or keep the referenced
symbol available and unchanged.

## Simulated residual envelopes

The function `vasicek_envelope()` constructs pointwise simulated envelopes
for normalized randomized quantile residuals and generalized Cox--Snell
residuals. In every bootstrap replication, a new response is generated under
the fitted model, the same GAMLSS model is re-estimated, and the residuals are
recalculated and ordered. This accounts for parameter-estimation variability.

For augmented families, the probability integral transform is randomized over
the fitted CDF jump at zero and/or one. Quantile residuals are compared with
the standard normal distribution, whereas Cox--Snell residuals are compared
with the unit-mean exponential distribution. Cox--Snell residuals are
evaluated as `-pnorm(rq, lower.tail = FALSE, log.p = TRUE)`, which is
equivalent to `-log(1 - pnorm(rq))` but remains stable in the upper tail.

```r
envelope_fit <- vasicek_envelope(
  object = fit_mean,
  residual = c("quantile", "cox-snell"),
  nsim = 200,
  level = 0.95,
  envelope = "quantile",
  seed = 123
)

plot(envelope_fit, which = "quantile")
plot(envelope_fit, which = "cox-snell")
```

The plots are full quantile--quantile plots rather than half-normal plots.
The default envelope uses pointwise empirical percentiles and is not a
simultaneous confidence band. To use the minimum and maximum of the ordered
simulated residuals at each position, set `envelope = "minmax"`. The blue
mean curve is calibrated to the fitted model and sample size; for Cox--Snell
residuals it may be more informative than the theoretical identity line in
the finite-sample upper tail. Failed or nonconverged bootstrap fits are
discarded and replaced up to the limit specified by `max.attempts`.
Automatic refitting requires the supplied data to contain exactly the
observations used in the original fit.

## Implementation

The density, cumulative distribution, and quantile functions call compiled
C++ routines through `Rcpp`. Random generation for `NVASIM` and `NVASIQ`
uses inverse transformation with the corresponding compiled quantile
functions, whereas `rLVASIQ()` and `rHVASIQ()` call compiled
random-generation routines directly.
The boundary-augmented functions reuse the compiled `NVASIM` functions and
add the required point masses in R.

The GAMLSS family definitions and analytical log-likelihood derivatives are
implemented in R. Numerical differentiation is not used. The mean and
variance components of `LVASIQ()` and `HVASIQ()` are evaluated by numerical
quadrature because these moments have no closed-form expressions.

## Testing

From the package root directory:

```bash
Rscript -e 'testthat::test_local(path = ".", reporter = "summary")'

cd ..
R CMD build vasicekreg
R CMD check vasicekreg_1.2.0.tar.gz
```

## Citation

To cite the package, use:

```r
citation("vasicekreg")
```

## References

### Vasicek distribution and regression

- Alves, B., Mazucheli, J., Leiva, V., and Bustillos, F. (2026). Vasicek mean regression: Estimation and diagnostics. *Under review*.

- Fischer, M. J., Hui, A. and Hösle, S. (2017). wHS-type distributions
  with application to finance. *Journal of Statistics and Management
  Systems*, **20**(1), 67--89.
  [doi:10.1080/09720510.2016.1190575](https://doi.org/10.1080/09720510.2016.1190575)

- Mazucheli, J. and Alves, B. (2026a). A hyperbolic-secant-kernel Vasicek-type quantile regression for rates and proportions. *Under review*.

- Mazucheli, J. and Alves, B. (2026b). A logistic-kernel Vasicek-type quantile regression for rates and proportions. *Under review*.

- Mazucheli, J. and Alves, B. (2026c). Augmented Vasicek mean regression models for rates and proportions. *Under review*.

- Mazucheli, J., Alves, B., Korkmaz, M. Ç., and Leiva, V. (2022).
  Vasicek quantile and mean regression models for bounded data: New
  formulation, mathematical derivations, and numerical applications.
  *Mathematics*, **10**, 1389.
  [doi:10.3390/math10091389](https://doi.org/10.3390/math10091389)

- Vasicek, O. A. (1987). *Probability of Loss on Loan Portfolio*. KMV
  Corporation.

- Vasicek, O. A. (2002). The distribution of loan portfolio value.
  *Risk*, **15**(12), 160--162.

- Witzany, J. (2013). A note on the Vasicek's model with the logistic
  distribution. *Ekonomický časopis (Journal of Economics)*, **61**(10),
  1053--1066.
  [IES Working Paper 1/2013](https://ideas.repec.org/p/fau/wpaper/wp2013_01.html)

### GAMLSS framework

- Rigby, R. A., and Stasinopoulos, D. M. (2005). Generalized additive
  models for location, scale and shape. *Journal of the Royal Statistical
  Society: Series C (Applied Statistics)*, **54**(3), 507--554.
  [doi:10.1111/j.1467-9876.2005.00510.x](https://doi.org/10.1111/j.1467-9876.2005.00510.x)

- Stasinopoulos, D. M., and Rigby, R. A. (2007). Generalized additive
  models for location scale and shape (GAMLSS) in R. *Journal of
  Statistical Software*, **23**(7), 1--46.
  [doi:10.18637/jss.v023.i07](https://doi.org/10.18637/jss.v023.i07)

### Boundary-augmented distributions

- Ospina, R., and Ferrari, S. L. P. (2010). Inflated beta distributions.
  *Statistical Papers*, **51**(1), 111--126.
  [doi:10.1007/s00362-008-0125-4](https://doi.org/10.1007/s00362-008-0125-4)

- Ospina, R., and Ferrari, S. L. P. (2012). A general class of zero-or-one
  inflated beta regression models. *Computational Statistics & Data
  Analysis*, **56**(6), 1609--1623.
  [doi:10.1016/j.csda.2011.10.005](https://doi.org/10.1016/j.csda.2011.10.005)

### Residuals and simulated envelopes

- Cox, D. R., and Snell, E. J. (1968). A general definition of residuals.
  *Journal of the Royal Statistical Society: Series B (Methodological)*,
  **30**(2), 248--275.
  [doi:10.1111/j.2517-6161.1968.tb00724.x](https://doi.org/10.1111/j.2517-6161.1968.tb00724.x)

- Dunn, P. K., and Smyth, G. K. (1996). Randomized quantile residuals.
  *Journal of Computational and Graphical Statistics*, **5**(3), 236--244.
  [doi:10.1080/10618600.1996.10474708](https://doi.org/10.1080/10618600.1996.10474708)

- Moral, R. A., Hinde, J., and Demetrio, C. G. B. (2017). Half-normal plots
  and overdispersed models in R: The hnp package. *Journal of Statistical
  Software*, **81**(10), 1--23.
  [doi:10.18637/jss.v081.i10](https://doi.org/10.18637/jss.v081.i10)

- Atkinson, A. C. (1985). *Plots, Transformations, and Regression: An
  Introduction to Graphical Methods of Diagnostic Regression Analysis*.
  Clarendon Press.

- Zhao, Y., Lee, A. H., Yau, K. K. W., and McLachlan, G. J. (2011).
  Assessing the adequacy of Weibull survival models: A simulated envelope
  approach. *Journal of Applied Statistics*, **38**(10), 2089--2097.
  [doi:10.1080/02664763.2010.545115](https://doi.org/10.1080/02664763.2010.545115)

### Included data sources

- Gange, S. J., Munoz, A., Saez, M., and Alonso, J. (1996). Use of the
  beta-binomial distribution to model the effect of policy changes on
  appropriateness of hospital stays. *Applied Statistics*, **45**(3),
  371--382.

- Korosteleva, O. (2019). *Advanced Regression Models with SAS and R*.
  CRC Press.

- Mazucheli, J., Leiva, V., Alves, B., and Menezes, A. F. B. (2021). A new
  quantile regression for modeling bounded data under a unit
  Birnbaum--Saunders distribution with applications in medicine and politics.
  *Symmetry*, **13**(4), 1--21.

- Menezes, A. F. B., Mazucheli, J., and Bourguignon, M. (2021). A parametric
  quantile regression approach for modeling zero- or one-inflated double
  bounded data. *Biometrical Journal*, **63**(4), 841--858.
  [doi:10.1002/bimj.202000126](https://doi.org/10.1002/bimj.202000126)

- Petterle, R. R., Bonat, W. H., Scarpin, C. T., Jonasson, T., and
  Borba, V. Z. C. (2020). Multivariate quasi-beta regression models for
  continuous bounded data. *The International Journal of Biostatistics*,
  **17**(1), 39--53.

- Stasinopoulos, M., and Rigby, R. (2025). *gamlss.data: Data for Generalized
  Additive Models for Location Scale and Shape*. R package version 6.0-7.

- Menezes, A. F. B. (2026). *uwquantreg: Unit-Weibull Quantile Regression*.
  R package version 0.1.0.

### Computational implementation

- Eddelbuettel, D., and François, R. (2011). Rcpp: Seamless R and C++
  integration. *Journal of Statistical Software*, **40**(8), 1--18.
  [doi:10.18637/jss.v040.i08](https://doi.org/10.18637/jss.v040.i08)

## License

`vasicekreg` is distributed under the MIT License.
