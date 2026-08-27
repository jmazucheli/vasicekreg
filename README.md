# vasicekreg

[![CRAN status](https://www.r-pkg.org/badges/version/vasicekreg)](https://CRAN.R-project.org/package=vasicekreg)

`vasicekreg` provides distribution functions and GAMLSS regression families
for Vasicek-type distributions on the unit interval. The package supports mean
and quantile regression under the standard normal kernel and quantile
regression under the logistic kernel. Normal-kernel families augmented at
zero, at one, or at both boundaries are available for responses containing
exact boundary values.

## Available families

| Family | Kernel | Response range | Parameterization | Interpretation of `mu` |
|---|---|---|---|---|
| `NVASIM` | Standard normal | `(0, 1)` | Mean | `mu = E(Y)` |
| `NVASIQ` | Standard normal | `(0, 1)` | Quantile | `mu = Q_Y(tau)` |
| `LVASIQ` | Logistic | `(0, 1)` | Quantile | `mu = Q_Y(tau)` |
| `ZANVASIM` | Standard normal | `[0, 1)` | Zero-adjusted mean | `mu = E(Y | Y > 0)` |
| `OANVASIM` | Standard normal | `(0, 1]` | One-adjusted mean | `mu = E(Y | Y < 1)` |
| `ZOANVASIM` | Standard normal | `[0, 1]` | Zero-and-one-adjusted mean | `mu = E(Y | 0 < Y < 1)` |

For all six families, `sigma` lies in `(0, 1)` and controls dispersion.
For `ZANVASIM`, `nu = P(Y = 0)` and the marginal mean is
`E(Y) = (1 - nu) * mu`.
For `OANVASIM`, `nu = P(Y = 1)` and the marginal mean is
`E(Y) = nu + (1 - nu) * mu`. In `ZOANVASIM`,
`nu = P(Y = 0)` and `tau = P(Y = 1 | Y > 0)`.
Its marginal mean is
`E(Y) = (1 - nu) * (tau + (1 - tau) * mu)`.
The logistic-kernel mean does not have a closed-form expression and generally
differs from `mu`; consequently, the package does not provide a
logistic-kernel mean-regression family.

Each parameterization includes density, cumulative distribution, quantile,
and random generation functions:

- `dNVASIM()`, `pNVASIM()`, `qNVASIM()`, and `rNVASIM()`;
- `dNVASIQ()`, `pNVASIQ()`, `qNVASIQ()`, and `rNVASIQ()`;
- `dLVASIQ()`, `pLVASIQ()`, `qLVASIQ()`, and `rLVASIQ()`;
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

For the quantile parameterizations, `tau` is supplied directly to the
distribution functions:

```r
qNVASIQ(p = 0.25, mu = 0.60, sigma = 0.25, tau = 0.25)
qLVASIQ(p = 0.25, mu = 0.60, sigma = 0.25, tau = 0.25)
```

Both calls return `mu = 0.60` because `mu` represents the `tau`-th quantile.

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

## GAMLSS zero-adjusted mean regression

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

## GAMLSS one-adjusted and zero-and-one-adjusted regression

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
`tau` is a model parameter and is unrelated to the global quantile level used
by `NVASIQ()` and `LVASIQ()`.

## GAMLSS quantile regression

For `NVASIQ()` and `LVASIQ()`, the quantile level must be defined as a scalar
variable named `tau` in the global environment. It is not passed as an
argument to the GAMLSS family constructor.

### Normal kernel

```r
library(gamlss)
library(vasicekreg)

set.seed(123)
tau <- 0.50

dat_normal <- data.frame(
  y = rNVASIQ(
    n = 300,
    mu = 0.60,
    sigma = 0.25,
    tau = tau
  )
)

fit_normal <- gamlss(
  y ~ 1,
  data = dat_normal,
  family = NVASIQ(
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
tau <- 0.25

dat_logistic <- data.frame(
  y = rLVASIQ(
    n = 300,
    mu = 0.60,
    sigma = 0.25,
    tau = tau
  )
)

fit_logistic <- gamlss(
  y ~ 1,
  data = dat_logistic,
  family = LVASIQ(
    mu.link = "logit",
    sigma.link = "logit"
  ),
  control = gamlss.control(trace = FALSE)
)

fitted(fit_logistic, what = "mu")[1]
```

The global value of `tau` must remain equal to the quantile level associated
with the fitted model when residuals or other post-fit quantities are
computed. Reset `tau` before working with a model fitted at another quantile
level.

## Simulated residual envelopes

The function `vasicek_envelope()` constructs pointwise simulated envelopes
for normalized randomized quantile residuals and generalized Cox--Snell
residuals. In every bootstrap replication, a new response is generated under
the fitted model, the same GAMLSS model is re-estimated, and the residuals are
recalculated and ordered. This accounts for parameter-estimation variability.

For adjusted families, the probability integral transform is randomized over
the fitted CDF jump at zero and/or one. Quantile residuals are compared with
the standard normal distribution, whereas Cox--Snell residuals are compared
with the unit-mean exponential distribution.

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

The default envelope uses pointwise empirical percentiles. To use the minimum
and maximum of the ordered simulated residuals at each position, set
`envelope = "minmax"`. Failed or nonconverged bootstrap fits are discarded
and replaced up to the limit specified by `max.attempts`.

## Implementation

The density, cumulative distribution, and quantile functions call compiled
C++ routines through `Rcpp`. Random generation for `NVASIM` and `NVASIQ`
uses inverse transformation with the corresponding compiled quantile
functions, whereas `rLVASIQ()` calls a compiled random-generation routine
directly.
The boundary-adjusted functions reuse the compiled `NVASIM` functions and
add the required point masses in R.

The GAMLSS family definitions and analytical log-likelihood derivatives are
implemented in R. Numerical differentiation is not used. The mean and
variance components of `LVASIQ()` are evaluated by numerical quadrature
because these moments have no closed-form expressions.

## Testing

From the package root directory:

```bash
Rscript -e 'testthat::test_local(path = ".", reporter = "summary")'

cd ..
R CMD build vasicekreg
R CMD check vasicekreg_1.1.0.tar.gz
```

## Citation

To cite the package, use:

```r
citation("vasicekreg")
```

The main methodological reference is:

> Mazucheli, J., Alves, B., Korkmaz, M. Ç., and Leiva, V. (2022).
> Vasicek quantile and mean regression models for bounded data: New
> formulation, mathematical derivations, and numerical applications.
> *Mathematics*, **10**, 1389.
> https://doi.org/10.3390/math10091389

## License

`vasicekreg` is distributed under the MIT License.
