# One-augmented normal-kernel Vasicek-type distribution

Defines a one-augmented normal-kernel Vasicek-type distribution for
responses in \\(0,1\]\\. The parameter \\\nu\\ is the probability at
one. Conditional on an observation in \\(0,1)\\, the distribution is
`NVASIM` with mean \\\mu\\ and shape parameter \\\sigma\\.

## Usage

``` r
d1NVASIM(x, mu = 0.5, sigma = 0.5, nu = 0.1, log = FALSE)

p1NVASIM(q, mu = 0.5, sigma = 0.5, nu = 0.1, lower.tail = TRUE, log.p = FALSE)

q1NVASIM(p, mu = 0.5, sigma = 0.5, nu = 0.1, lower.tail = TRUE, log.p = FALSE)

r1NVASIM(n, mu = 0.5, sigma = 0.5, nu = 0.1)

dOANVASIM(x, mu = 0.5, sigma = 0.5, nu = 0.1, log = FALSE)

pOANVASIM(q, mu = 0.5, sigma = 0.5, nu = 0.1, lower.tail = TRUE, log.p = FALSE)

qOANVASIM(p, mu = 0.5, sigma = 0.5, nu = 0.1, lower.tail = TRUE, log.p = FALSE)

rOANVASIM(n, mu = 0.5, sigma = 0.5, nu = 0.1)

OANVASIM(mu.link = "logit", sigma.link = "logit", nu.link = "logit")
```

## Arguments

- x:

  Vector of values in \\\[0,1\]\\ at which the density or probability
  mass is evaluated. The distribution has support \\(0,1\]\\, and the
  returned value is zero at \\x=0\\.

- mu:

  Mean of the continuous Vasicek component, in \\(0,1)\\.

- sigma:

  Shape parameter of the continuous Vasicek component, in \\(0,1)\\.

- nu:

  Probability at one, in \\(0,1)\\.

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

`OANVASIM()` returns a `gamlss.family` object. The functions
`d1NVASIM()`, `p1NVASIM()`, `q1NVASIM()`, and `r1NVASIM()` return
probability mass or density values, cumulative probabilities, quantiles,
and random observations, respectively. `dOANVASIM()`, `pOANVASIM()`,
`qOANVASIM()`, and `rOANVASIM()` are equivalent names following the
GAMLSS family-name convention.

## Details

Let \\Y_c\sim\mathrm{NVASIM}(\mu,\sigma)\\ and let \\0\<\nu\<1\\. The
BEOI-type one-augmented distribution is defined by \$\$P(Y=1)=\nu\$\$
and \$\$f_Y(y)=(1-\nu)f\_{Y_c}(y\mid\mu,\sigma),\quad 0\<y\<1.\$\$
Consequently, \$\$E(Y)=\nu+(1-\nu)\mu\$\$ and
\$\$\mathrm{Var}(Y)=(1-\nu)\mathrm{Var}(Y_c)+ \nu(1-\nu)(1-\mu)^2.\$\$
Thus, \\\mu=E(Y\mid 0\<Y\<1)\\ is the mean of the continuous component,
whereas \\\nu+(1-\nu)\mu\\ is the marginal mean.

## References

Ospina, R. and Ferrari, S. L. P. (2010). Inflated beta distributions.
*Statistical Papers*, **51**, 111–126.

Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive
models for location, scale and shape. *Applied Statistics*, **54**(3),
507–554.

## See also

[`NVASIM`](https://jmazucheli.github.io/vasicekreg/reference/NVASIM.md),
[`ZANVASIM`](https://jmazucheli.github.io/vasicekreg/reference/ZANVASIM.md),
[`BEOI`](https://rdrr.io/pkg/gamlss.dist/man/BEOI.html)

## Examples

``` r
set.seed(123)
y <- r1NVASIM(1000, mu = 0.60, sigma = 0.30, nu = 0.20)
mean(y == 1)
#> [1] 0.198
mean(y)
#> [1] 0.6782138
0.20 + (1 - 0.20) * 0.60
#> [1] 0.68

if (FALSE) { # \dontrun{
library(gamlss)
fit <- gamlss(
  y ~ 1,
  sigma.formula = ~ 1,
  nu.formula = ~ 1,
  family = OANVASIM(),
  control = gamlss.control(trace = FALSE)
)
} # }
```
