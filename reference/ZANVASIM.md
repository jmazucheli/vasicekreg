# Zero-augmented normal-kernel Vasicek-type distribution

Defines a zero-augmented normal-kernel Vasicek-type distribution for
responses in \\\[0,1)\\. The parameter \\\nu\\ is the probability of a
structural zero. Conditional on a positive response, the distribution is
`NVASIM` with mean \\\mu\\ and shape parameter \\\sigma\\.

## Usage

``` r
d0NVASIM(x, mu = 0.5, sigma = 0.5, nu = 0.1, log = FALSE)

p0NVASIM(q, mu = 0.5, sigma = 0.5, nu = 0.1, lower.tail = TRUE, log.p = FALSE)

q0NVASIM(p, mu = 0.5, sigma = 0.5, nu = 0.1, lower.tail = TRUE, log.p = FALSE)

r0NVASIM(n, mu = 0.5, sigma = 0.5, nu = 0.1)

dZANVASIM(x, mu = 0.5, sigma = 0.5, nu = 0.1, log = FALSE)

pZANVASIM(q, mu = 0.5, sigma = 0.5, nu = 0.1, lower.tail = TRUE, log.p = FALSE)

qZANVASIM(p, mu = 0.5, sigma = 0.5, nu = 0.1, lower.tail = TRUE, log.p = FALSE)

rZANVASIM(n, mu = 0.5, sigma = 0.5, nu = 0.1)

ZANVASIM(mu.link = "logit", sigma.link = "logit", nu.link = "logit")
```

## Arguments

- x:

  Vector of values in \\\[0,1\]\\ at which the density or probability
  mass is evaluated. The distribution has support \\\[0,1)\\, and the
  returned value is zero at \\x=1\\.

- mu:

  Mean of the positive Vasicek component, in \\(0,1)\\.

- sigma:

  Shape parameter of the positive Vasicek component, in \\(0,1)\\.

- nu:

  Probability of a structural zero, in \\(0,1)\\.

- log:

  Logical; if `TRUE`, log probabilities or log densities are returned.

- q:

  Vector of values in \\\[0,1\]\\ at which the cumulative distribution
  function is evaluated.

- lower.tail:

  Logical; if `TRUE`, probabilities are \\P(Y\leq y)\\; otherwise, they
  are \\P(Y\>y)\\.

- log.p:

  Logical; if `TRUE`, probabilities are supplied or returned on the log
  scale.

- p:

  Vector of probabilities.

- n:

  Number of observations. If `length(n) > 1`, its length is taken to be
  the number required.

- mu.link:

  Link function for \\\mu\\.

- sigma.link:

  Link function for \\\sigma\\.

- nu.link:

  Link function for \\\nu\\.

## Value

`ZANVASIM()` returns a `gamlss.family` object. The functions
`d0NVASIM()`, `p0NVASIM()`, `q0NVASIM()`, and `r0NVASIM()` return
density or probability mass values, cumulative probabilities, quantiles,
and random observations, respectively. `dZANVASIM()`, `pZANVASIM()`,
`qZANVASIM()`, and `rZANVASIM()` are equivalent names following the
GAMLSS family-name convention.

## Details

Let \\Y\_+\sim\mathrm{NVASIM}(\mu,\sigma)\\ and let \\0\<\nu\<1\\. The
zero-augmented distribution is defined by \$\$P(Y=0)=\nu\$\$ and
\$\$f_Y(y)=(1-\nu)f\_{Y\_+}(y\mid\mu,\sigma),\quad 0\<y\<1.\$\$ Its
cumulative distribution function is
\$\$F_Y(y)=\nu+(1-\nu)F\_{Y\_+}(y\mid\mu,\sigma),\quad 0\<y\<1.\$\$
Consequently, \$\$E(Y)=(1-\nu)\mu\$\$ and
\$\$\mathrm{Var}(Y)=(1-\nu)\mathrm{Var}(Y\_+)+ \nu(1-\nu)\mu^2.\$\$

Thus, \\\mu\\ is the mean conditional on \\Y\>0\\; it is not the
marginal mean when \\\nu\>0\\. The marginal mean is \\(1-\nu)\mu\\.

## References

Mazucheli, J., Alves, B., Korkmaz, M. C., and Leiva, V. (2022). Vasicek
quantile and mean regression models for bounded data: New formulation,
mathematical derivations, and numerical applications. *Mathematics*,
**10**, 1389.
[doi:10.3390/math10091389](https://doi.org/10.3390/math10091389)

Ospina, R. and Ferrari, S. L. P. (2010). Inflated beta distributions.
*Statistical Papers*, **51**, 111–126.

Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive
models for location, scale and shape. *Applied Statistics*, **54**(3),
507–554.

## See also

[`NVASIM`](https://jmazucheli.github.io/vasicekreg/reference/NVASIM.md),
[`BEZI`](https://rdrr.io/pkg/gamlss.dist/man/BEZI.html)

## Examples

``` r
set.seed(123)
y <- r0NVASIM(1000, mu = 0.60, sigma = 0.30, nu = 0.20)
mean(y == 0)
#> [1] 0.198
mean(y)
#> [1] 0.478592
(1 - 0.20) * 0.60
#> [1] 0.48

library(gamlss)
fit <- gamlss(
  y ~ 1,
  sigma.formula = ~ 1,
  nu.formula = ~ 1,
  family = ZANVASIM(),
  control = gamlss.control(trace = FALSE)
)
fitted(fit, what = "mu")[1]
#> [1] 0.5972124
fitted(fit, what = "sigma")[1]
#> [1] 0.2947561
fitted(fit, what = "nu")[1]
#> [1] 0.198
```
