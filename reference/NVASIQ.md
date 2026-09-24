# Normal-kernel Vasicek-type distribution with quantile parameterization

The function `NVASIQ()` defines the normal-kernel Vasicek-type
distribution as a `gamlss.family` object. In this parameterization,
\\\mu\\ corresponds to the fixed \\\tau\\-th quantile and \\\sigma\\ is
a shape parameter. The fixed level is supplied through the `quantile`
argument. The functions `dNVASIQ`, `pNVASIQ`, `qNVASIQ`, and `rNVASIQ`
define the density, distribution function, quantile function, and random
generation for the Vasicek distribution, respectively.

## Usage

``` r
dNVASIQ(x, mu, sigma, quantile = 0.5, log = FALSE)

pNVASIQ(q, mu, sigma, quantile = 0.5, lower.tail = TRUE, log.p = FALSE)

qNVASIQ(p, mu, sigma, quantile = 0.5, lower.tail = TRUE, log.p = FALSE)

rNVASIQ(n, mu, sigma, quantile = 0.5)

NVASIQ(quantile = 0.5, mu.link = "logit", sigma.link = "logit")
```

## Arguments

- x:

  Vector of quantiles in the interval \\(0,1)\\.

- mu:

  Vector of \\\tau\\-th quantile parameter values.

- sigma:

  Vector of shape parameter values.

- quantile:

  Fixed quantile level \\\tau\in(0,1)\\ used in the distribution
  functions and in the `NVASIQ()` GAMLSS family.

- log, log.p:

  Logical; if `TRUE`, probabilities are returned on the log scale.

- q:

  Vector of values in \\\[0,1\]\\ at which the cumulative distribution
  function is evaluated.

- lower.tail:

  Logical; if `TRUE` (default), \\P(X \le x)\\ is returned; otherwise,
  \\P(X \> x)\\.

- p:

  Vector of probabilities in \\\[0,1\]\\ on the probability scale.

- n:

  Number of observations. If `length(n) > 1`, the length is taken to be
  the number required.

- mu.link:

  Link function for the \\\mu\\ parameter.

- sigma.link:

  Link function for the \\\sigma\\ parameter.

## Value

`NVASIQ()` returns a `gamlss.family` object that can be used to fit a
Vasicek-type distribution using the
[`gamlss`](https://rdrr.io/pkg/gamlss/man/gamlss.html) function.

## Details

Probability density function: \$\$f\left(x \mid \mu, \sigma, \tau\right)
= \sqrt{\frac{1-\sigma}{\sigma}}
\exp\left\\\frac{1}{2}\left\[\Phi^{-1}(x)^2 -
\left(\frac{\sqrt{1-\sigma}\left(\Phi^{-1}(x)-\Phi^{-1}(\mu)\right) -
\sqrt{\sigma}\\\Phi^{-1}(\tau)}{\sqrt{\sigma}}\right)^2\right\]\right\\.\$\$

Cumulative distribution function: \$\$F\left(x \mid \mu, \sigma,
\tau\right) =
\Phi\left(\frac{\sqrt{1-\sigma}\left(\Phi^{-1}(x)-\Phi^{-1}(\mu)\right) -
\sqrt{\sigma}\\\Phi^{-1}(\tau)}{\sqrt{\sigma}}\right).\$\$

where \\0\<x\<1\\, \\0\<\mu\<1\\, \\0\<\sigma\<1\\, and \\0\<\tau\<1\\;
\\\mu\\ is the \\\tau\\-th quantile and \\\sigma\\ is the shape
parameter.

## Note

For `NVASIQ()`, \\\mu\\ corresponds to the \\\tau\\-th quantile and
\\\sigma\\ is a shape parameter. The level supplied through `quantile`
is stored in the family definition and embedded as a numeric literal in
the family components used by GAMLSS; no global variable is required.

## References

Hastie, T. J. and Tibshirani, R. J. (1990). *Generalized Additive
Models*. Chapman and Hall, London.

Mazucheli, J., Alves, B., Korkmaz, M. Ç., and Leiva, V. (2022). Vasicek
quantile and mean regression models for bounded data: New formulation,
mathematical derivations, and numerical applications. *Mathematics*,
**10**, 1389.

Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive
models for location, scale and shape (with discussion). *Applied
Statistics*, **54**(3), 507–554.

Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F.
(2019). *Distributions for Modeling Location, Scale, and Shape: Using
GAMLSS in R*. Chapman and Hall/CRC.

Stasinopoulos, D. M. and Rigby, R. A. (2007). Generalized additive
models for location, scale and shape (GAMLSS) in R. *Journal of
Statistical Software*, **23**(7), 1–46.

Stasinopoulos, D. M., Rigby, R. A., Heller, G., Voudouris, V., and De
Bastiani, F. (2017). *Flexible Regression and Smoothing: Using GAMLSS in
R*. Chapman and Hall/CRC.

Vasicek, O. A. (1987). Probability of loss on loan portfolio. *KMV
Corporation*.

Vasicek, O. A. (2002). The distribution of loan portfolio value. *Risk*,
**15**(12), 160–162.

## See also

[`NVASIM`](https://jmazucheli.github.io/vasicekreg/reference/NVASIM.md)

## Author

Josmar Mazucheli <jmazucheli@gmail.com>

Bruna Alves <pg402900@uem.br>

## Examples

``` r
set.seed(123)
x <- rNVASIQ(n = 1000, mu = 0.50, sigma = 0.69, quantile = 0.50)
R <- range(x)
S <- seq(from = R[1], to = R[2], length.out = 1000)

hist(x, prob = TRUE, main = "Vasicek")
lines(S, dNVASIQ(x = S, mu = 0.50, sigma = 0.69, quantile = 0.50), col = 2)


plot(ecdf(x))
lines(S, pNVASIQ(q = S, mu = 0.50, sigma = 0.69, quantile = 0.50), col = 2)


plot(quantile(x, probs = S), type = "l")
lines(qNVASIQ(p = S, mu = 0.50, sigma = 0.69, quantile = 0.50), col = 2)


library(gamlss)
set.seed(123)
data <- data.frame(
  y = rNVASIQ(n = 100, mu = 0.50, sigma = 0.69, quantile = 0.50)
)

fit <- gamlss(y ~ 1, data = data,
              family = NVASIQ(quantile = 0.50,
                             mu.link = "logit",
                             sigma.link = "logit"))
#> GAMLSS-RS iteration 1: Global Deviance = -35.9846 
#> GAMLSS-RS iteration 2: Global Deviance = -35.9846 
1 / (1 + exp(-fit$mu.coefficients))
#> (Intercept) 
#>   0.4967521 
1 / (1 + exp(-fit$sigma.coefficients))
#> (Intercept) 
#>   0.6777598 

set.seed(123)
n <- 100
x <- rbinom(n, size = 1, prob = 0.5)
eta <- 0.5 + 1 * x
mu <- 1 / (1 + exp(-eta))
sigma <- 0.5
y <- rNVASIQ(n, mu, sigma, quantile = 0.5)
data <- data.frame(y, x)

fit_median <- gamlss(
  y ~ x, data = data, family = NVASIQ(quantile = 0.50)
)
#> GAMLSS-RS iteration 1: Global Deviance = -55.4043 
#> GAMLSS-RS iteration 2: Global Deviance = -55.4042 

summary(fit_median)
#> ******************************************************************
#> Family:  c("NVASIQ", "Normal-kernel Vasicek-type quantile") 
#> 
#> Call:  gamlss(formula = y ~ x, family = NVASIQ(quantile = 0.5),  
#>     data = data) 
#> 
#> Fitting method: RS() 
#> 
#> ------------------------------------------------------------------
#> Mu link function:  logit
#> Mu Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)   
#> (Intercept)   0.6554     0.1936   3.385  0.00103 **
#> x             0.9284     0.2975   3.120  0.00238 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> ------------------------------------------------------------------
#> Sigma link function:  logit
#> Sigma Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)  
#> (Intercept)  -0.2929     0.1414  -2.071    0.041 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> ------------------------------------------------------------------
#> No. of observations in the fit:  100 
#> Degrees of Freedom for the fit:  3
#>       Residual Deg. of Freedom:  97 
#>                       at cycle:  2 
#>  
#> Global Deviance:     -55.40423 
#>             AIC:     -49.40423 
#>             SBC:     -41.58872 
#> ******************************************************************
```
