# Hyperbolic-secant-kernel Vasicek-type quantile distribution

`HVASIQ()` defines the hyperbolic-secant-kernel Vasicek-type
distribution as a `gamlss.family` object for conditional quantile
regression. The functions `dHVASIQ()`, `pHVASIQ()`, `qHVASIQ()`, and
`rHVASIQ()` provide the density, distribution function, quantile
function, and random generation.

## Usage

``` r
dHVASIQ(x, mu, sigma, quantile = 0.5, log = FALSE)

pHVASIQ(q, mu, sigma, quantile = 0.5, lower.tail = TRUE, log.p = FALSE)

qHVASIQ(p, mu, sigma, quantile = 0.5, lower.tail = TRUE, log.p = FALSE)

rHVASIQ(n, mu, sigma, quantile = 0.5)

HVASIQ(quantile = 0.5, mu.link = "logit", sigma.link = "logit")
```

## Arguments

- x:

  Vector of quantiles in \\(0,1)\\.

- mu:

  Vector of conditional \\\tau\\-th quantiles, \\0\<\mu\<1\\.

- sigma:

  Vector of shape parameter values, \\0\<\sigma\<1\\.

- quantile:

  Fixed quantile level \\\tau\in(0,1)\\ represented by \\\mu\\, used in
  the distribution functions and in `HVASIQ()`.

- log:

  Logical; if `TRUE`, returns the log-density.

- q:

  Vector of values in \\\[0,1\]\\ at which the cumulative distribution
  function is evaluated.

- lower.tail:

  Logical; if `TRUE`, probabilities are \\P(X\leq x)\\.

- log.p:

  Logical; if `TRUE`, probabilities are supplied or returned on the
  logarithmic scale.

- p:

  Vector of probabilities in \\\[0,1\]\\ on the probability scale.

- n:

  Number of observations.

- mu.link:

  Link function for \\\mu\\.

- sigma.link:

  Link function for \\\sigma\\.

## Value

`dHVASIQ()` returns the density, `pHVASIQ()` the distribution function,
`qHVASIQ()` the quantile function, and `rHVASIQ()` random deviates.
`HVASIQ()` returns a `gamlss.family` object.

## Details

Let \$\$H(w)=\frac{2}{\pi}\arctan\\\exp(w)\\,\qquad
h(w)=\frac{1}{\pi\cosh(w)},\$\$ and \$\$q(u)=H^{-1}(u)=
\log\left\\\tan\left(\frac{\pi u}{2}\right)\right\\.\$\$ For fixed
\\\tau\in(0,1)\\, define \$\$a=\sqrt{\frac{1-\sigma}{\sigma}},\qquad
z=a\\q(x)-q(\mu)\\+q(\tau).\$\$ The cumulative distribution function and
density are \$\$F(x\mid\mu,\sigma,\tau)=H(z)\$\$ and
\$\$f(x\mid\mu,\sigma,\tau)=a\frac{h(z)}{h\\q(x)\\}.\$\$ The quantile
function is \$\$Q(p\mid\mu,\sigma,\tau)=
H\left\[q(\mu)+\sqrt{\frac{\sigma}{1-\sigma}}
\\q(p)-q(\tau)\\\right\].\$\$ Consequently, \\Q(\tau)=\mu\\; \\\mu\\ is
exactly the conditional \\\tau\\-th quantile and is not the conditional
mean in general.

The GAMLSS family uses analytical derivatives. For one observation, let
\$\$d=q(y)-q(\mu),\quad z=ad+q(\tau),\quad T=\tanh(z),\quad
S=\operatorname{sech}^2(z),\$\$ \$\$b=\frac{1}{2\sigma(1-\sigma)},\qquad
g=\frac{\pi}{\sin(\pi\mu)}.\$\$ If \\\ell\\ is the individual
log-likelihood contribution, then
\$\$\frac{\partial\ell}{\partial\mu}=agT,\$\$
\$\$\frac{\partial\ell}{\partial\sigma}=b(adT-1),\$\$
\$\$\frac{\partial^2\ell}{\partial\mu^2}
=g^2\\-a\cos(\pi\mu)T-a^2S\\,\$\$
\$\$\frac{\partial^2\ell}{\partial\mu\\\partial\sigma} =-abg(T+adS),\$\$
and \$\$\frac{\partial^2\ell}{\partial\sigma^2}
=b^2\\2(1-2\sigma)-(3-4\sigma)adT-a^2d^2S\\.\$\$ The conditional mean
and variance do not have elementary closed forms and are evaluated by
numerical integration of the quantile function.

The fixed level is supplied through the `quantile` argument. It is
stored in the family definition and embedded as a numeric literal in the
family components used by GAMLSS; no global variable is required.

## References

Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals.
*Journal of Computational and Graphical Statistics*, **5**(3), 236–244.

Fischer, M. J., Hui, A. and Hösle, S. (2017). wHS-type distributions
with application to finance. *Journal of Statistics and Management
Systems*, **20**(1), 67–89.
[doi:10.1080/09720510.2016.1190575](https://doi.org/10.1080/09720510.2016.1190575)

Mazucheli, J., Alves, B., Korkmaz, M. C. and Leiva, V. (2022). Vasicek
quantile and mean regression models for bounded data: New formulation,
mathematical derivations, and numerical applications. *Mathematics*,
**10**, 1389.
[doi:10.3390/math10091389](https://doi.org/10.3390/math10091389)

## Author

Josmar Mazucheli <jmazucheli@gmail.com>

## Examples

``` r
set.seed(123)
y <- rHVASIQ(500, mu = 0.60, sigma = 0.30, quantile = 0.25)

fit <- gamlss::gamlss(
    y ~ 1,
    sigma.formula = ~ 1,
    family = HVASIQ(quantile = 0.25),
    control = gamlss::gamlss.control(trace = FALSE)
)
fitted(fit, what = "mu")[1]
#> [1] 0.5968944
```
