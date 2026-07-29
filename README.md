# vasicekreg

[![CRAN status](https://www.r-pkg.org/badges/version/vasicekreg)](https://CRAN.R-project.org/package=vasicekreg)

`vasicekreg` provides distribution functions and GAMLSS regression families
for Vasicek-type distributions on the unit interval. The package supports mean
and quantile regression under the standard normal kernel and quantile
regression under the logistic kernel.

## Available families

| Family | Kernel | Parameterization | Interpretation of `mu` |
|---|---|---|---|
| `NVASIM` | Standard normal | Mean | `mu = E(Y)` |
| `NVASIQ` | Standard normal | Quantile | `mu = Q_Y(tau)` |
| `LVASIQ` | Logistic | Quantile | `mu = Q_Y(tau)` |

For all three families, `sigma` lies in `(0, 1)` and controls dispersion.
The logistic-kernel mean does not have a closed-form expression and generally
differs from `mu`; consequently, the package does not provide a
logistic-kernel mean-regression family.

Each parameterization includes density, cumulative distribution, quantile,
and random generation functions:

- `dNVASIM()`, `pNVASIM()`, `qNVASIM()`, and `rNVASIM()`;
- `dNVASIQ()`, `pNVASIQ()`, `qNVASIQ()`, and `rNVASIQ()`;
- `dLVASIQ()`, `pLVASIQ()`, `qLVASIQ()`, and `rLVASIQ()`.

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

## Implementation

The density, cumulative distribution, and quantile functions call compiled
C++ routines through `Rcpp`. Random generation for `NVASIM` and `NVASIQ`
uses inverse transformation with the corresponding compiled quantile
functions, whereas `rLVASIQ()` calls a compiled random-generation routine
directly.

The GAMLSS family definitions and analytical log-likelihood derivatives are
implemented in R. Numerical differentiation is not used. The mean and
variance components of `LVASIQ()` are evaluated by numerical quadrature
because these moments have no closed-form expressions.

## Testing

From the package root directory:

```r
devtools::test()
devtools::check()
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
