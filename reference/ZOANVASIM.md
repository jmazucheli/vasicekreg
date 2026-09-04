# Zero-and-one-augmented normal-kernel Vasicek-type distribution

Defines a normal-kernel Vasicek-type distribution augmented by point
masses at zero and one. Conditional on an observation in \\(0,1)\\, the
continuous component is `NVASIM` with mean \\\mu\\ and shape parameter
\\\sigma\\.

## Usage

``` r
d01NVASIM(x, mu = 0.5, sigma = 0.5, nu = 0.1, tau = 0.1, log = FALSE)

p01NVASIM(
  q,
  mu = 0.5,
  sigma = 0.5,
  nu = 0.1,
  tau = 0.1,
  lower.tail = TRUE,
  log.p = FALSE
)

q01NVASIM(
  p,
  mu = 0.5,
  sigma = 0.5,
  nu = 0.1,
  tau = 0.1,
  lower.tail = TRUE,
  log.p = FALSE
)

r01NVASIM(n, mu = 0.5, sigma = 0.5, nu = 0.1, tau = 0.1)

dZOANVASIM(x, mu = 0.5, sigma = 0.5, nu = 0.1, tau = 0.1, log = FALSE)

pZOANVASIM(
  q,
  mu = 0.5,
  sigma = 0.5,
  nu = 0.1,
  tau = 0.1,
  lower.tail = TRUE,
  log.p = FALSE
)

qZOANVASIM(
  p,
  mu = 0.5,
  sigma = 0.5,
  nu = 0.1,
  tau = 0.1,
  lower.tail = TRUE,
  log.p = FALSE
)

rZOANVASIM(n, mu = 0.5, sigma = 0.5, nu = 0.1, tau = 0.1)

ZOANVASIM(
  mu.link = "logit",
  sigma.link = "logit",
  nu.link = "logit",
  tau.link = "logit"
)
```

## Arguments

- x:

  Vector of values in \\\[0,1\]\\ at which the density or probability
  mass is evaluated. The distribution has support \\\[0,1\]\\.

- mu:

  Mean of the continuous Vasicek component, in \\(0,1)\\.

- sigma:

  Shape parameter of the continuous Vasicek component, in \\(0,1)\\.

- nu:

  Probability at zero, \\\nu=P(Y=0)\\, in \\(0,1)\\.

- tau:

  Conditional probability at one among nonzero observations,
  \\\tau=P(Y=1\mid Y\>0)\\, in \\(0,1)\\. This parameter is unrelated to
  the fixed quantile level used by
  [`NVASIQ()`](https://jmazucheli.github.io/vasicekreg/reference/NVASIQ.md),
  [`LVASIQ()`](https://jmazucheli.github.io/vasicekreg/reference/LVASIQ.md),
  and
  [`HVASIQ()`](https://jmazucheli.github.io/vasicekreg/reference/HVASIQ.md).

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

- tau.link:

  Link function for \\\tau\\.

## Value

`ZOANVASIM()` returns a four-parameter `gamlss.family` object. The
functions `d01NVASIM()`, `p01NVASIM()`, `q01NVASIM()`, and `r01NVASIM()`
return probability mass or density values, cumulative probabilities,
quantiles, and random observations, respectively. `dZOANVASIM()`,
`pZOANVASIM()`, `qZOANVASIM()`, and `rZOANVASIM()` are equivalent names
following the GAMLSS family-name convention.

## Details

Let \\Y_c\sim\mathrm{NVASIM}(\mu,\sigma)\\. Write \\p_0\\, \\p_1\\, and
\\p_c\\ for the probabilities of zero, one, and the continuous
component. The sequential BEOI-type parameterization is
\$\$\nu=P(Y=0),\qquad \tau=P(Y=1\mid Y\>0).\$\$ Hence,
\$\$p_0=\nu,\qquad p_1=(1-\nu)\tau,\qquad p_c=(1-\nu)(1-\tau).\$\$ The
distribution is \$\$P(Y=0)=p_0,\qquad P(Y=1)=p_1\$\$ and \$\$f_Y(y)=p_c
f\_{Y_c}(y\mid\mu,\sigma),\quad 0\<y\<1.\$\$ Its marginal mean and
variance are \$\$E(Y)=(1-\nu)\left\[\tau+(1-\tau)\mu\right\]\$\$ and
\$\$\mathrm{Var}(Y)=
(1-\nu)\left\[(1-\tau)\left\\\mathrm{Var}(Y_c)+\mu^2\right\\
+\tau\right\]
-\left\\(1-\nu)\left\[\tau+(1-\tau)\mu\right\]\right\\^2.\$\$ Logit
links for \\\nu\\ and \\\tau\\ guarantee valid probabilities. In the
closure of the parameter space, \\\nu=0\\ gives the BEOI-type
one-augmented model and \\\tau=0\\ gives the zero-augmented model. The
implemented distribution functions and GAMLSS family use the open
parameter space \\0\<\nu\<1\\ and \\0\<\tau\<1\\; these nested models
are therefore limiting cases rather than admissible interior parameter
values.

## References

Ospina, R. and Ferrari, S. L. P. (2010). Inflated beta distributions.
*Statistical Papers*, **51**, 111–126.

Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive
models for location, scale and shape. *Applied Statistics*, **54**(3),
507–554.

## See also

[`NVASIM`](https://jmazucheli.github.io/vasicekreg/reference/NVASIM.md),
[`ZANVASIM`](https://jmazucheli.github.io/vasicekreg/reference/ZANVASIM.md),
[`OANVASIM`](https://jmazucheli.github.io/vasicekreg/reference/OANVASIM.md),
[`BEOI`](https://rdrr.io/pkg/gamlss.dist/man/BEOI.html)

## Examples

``` r
set.seed(123)
y <- r01NVASIM(
  1000, mu = 0.60, sigma = 0.30, nu = 0.20, tau = 0.25
)
c(zero = mean(y == 0), one = mean(y == 1))
#>  zero   one 
#> 0.198 0.198 
mean(y)
#> [1] 0.5587152
(1 - 0.20) * (0.25 + (1 - 0.25) * 0.60)
#> [1] 0.56

if (FALSE) { # \dontrun{
library(gamlss)
fit <- gamlss(
  y ~ 1,
  sigma.formula = ~ 1,
  nu.formula = ~ 1,
  tau.formula = ~ 1,
  family = ZOANVASIM(),
  control = gamlss.control(trace = FALSE)
)

## Real-data example with observations at both boundaries
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
} # }
```
