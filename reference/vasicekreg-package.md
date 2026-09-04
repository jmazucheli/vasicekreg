# Overview of the vasicekreg package

The vasicekreg package provides distribution functions and GAMLSS
families for Vasicek-type distributions on the unit interval. Four base
families are available:

- `NVASIM`: normal kernel with mean parameterization, where
  \\\mu=E(Y)\\.

- `NVASIQ`: normal kernel with quantile parameterization, where
  \\\mu=Q_Y(\tau)\\ for a fixed \\\tau\in(0,1)\\.

- `LVASIQ`: logistic kernel with quantile parameterization, where
  \\\mu=Q_Y(\tau)\\ for a fixed \\\tau\in(0,1)\\.

- `HVASIQ`: hyperbolic-secant kernel with quantile parameterization,
  where \\\mu=Q_Y(\tau)\\ for a fixed \\\tau\in(0,1)\\.

For responses observed at the boundaries, the normal-kernel mean model
is also available as:

- `ZANVASIM`: point mass at zero and a continuous component on
  \\(0,1)\\.

- `OANVASIM`: point mass at one and a continuous component on \\(0,1)\\.

- `ZOANVASIM`: point masses at zero and one and a continuous component
  on \\(0,1)\\.

The documentation uses *augmented* for these boundary mixtures and
*Vasicek-type* for kernel-based constructions. The established family
names are retained for backward compatibility and follow familiar GAMLSS
abbreviations in which `ZA` and `OA` historically denote zero- and
one-adjusted families. The shape parameter \\\sigma\in(0,1)\\ controls
dispersion in the continuous Vasicek component. The corresponding `d`,
`p`, `q`, and `r` functions provide density or probability mass values,
cumulative probabilities, quantiles, and random observations,
respectively.

## Details

Included datasets:

- [`bodyfat`](https://jmazucheli.github.io/vasicekreg/reference/bodyfat.md):
  body-fat proportions in \\(0,1)\\ and demographic covariates for 298
  individuals.

- [`aep`](https://jmazucheli.github.io/vasicekreg/reference/aep.md):
  hospital-stay data from 1,383 patients; `noinap / los` contains
  observations at zero and one.

- [`transport`](https://jmazucheli.github.io/vasicekreg/reference/transport.md):
  bicycle-trip proportions for 60 respondents, including observations at
  zero.

- [`trees`](https://jmazucheli.github.io/vasicekreg/reference/trees.md):
  two-year tree-survival proportions for 26 parks, including
  observations at one.

[`NVASIM`](https://jmazucheli.github.io/vasicekreg/reference/NVASIM.md):
Normal-kernel mean parameterization and GAMLSS family. In regression
models, covariates describe the conditional mean through \\\mu\\.

[`NVASIQ`](https://jmazucheli.github.io/vasicekreg/reference/NVASIQ.md):
Normal-kernel quantile parameterization and GAMLSS family. For a fixed
quantile level \\\tau\\, covariates describe the conditional \\\tau\\-th
quantile through \\\mu\\.

[`LVASIQ`](https://jmazucheli.github.io/vasicekreg/reference/LVASIQ.md):
Logistic-kernel quantile parameterization and GAMLSS family. For a fixed
quantile level \\\tau\\, covariates describe the conditional \\\tau\\-th
quantile through \\\mu\\. A logistic-kernel mean-regression family is
not provided because the mean has no closed-form expression and does not
equal \\\mu\\ under this parameterization.

[`HVASIQ`](https://jmazucheli.github.io/vasicekreg/reference/HVASIQ.md):
Hyperbolic-secant-kernel quantile parameterization and GAMLSS family.
For a fixed quantile level \\\tau\\, covariates describe the conditional
\\\tau\\-th quantile through \\\mu\\. Its conditional mean and variance
are obtained by numerical quadrature and \\\mu\\ must not be interpreted
as the mean.

[`ZANVASIM`](https://jmazucheli.github.io/vasicekreg/reference/ZANVASIM.md):
Zero-augmented normal-kernel mean family. Here \\\nu=P(Y=0)\\,
\\\mu=E(Y\mid Y\>0)\\, and the marginal mean is \\E(Y)=(1-\nu)\mu\\.

[`OANVASIM`](https://jmazucheli.github.io/vasicekreg/reference/OANVASIM.md):
One-augmented normal-kernel mean family. Here \\\nu=P(Y=1)\\,
\\\mu=E(Y\mid Y\<1)\\, and the marginal mean is \\E(Y)=\nu+(1-\nu)\mu\\.
The parameters \\\mu\\ and \\\nu\\ therefore have the same
interpretations as their counterparts in the one-inflated beta family
[`BEOI`](https://rdrr.io/pkg/gamlss.dist/man/BEOI.html). The shape
parameter \\\sigma\\ is distribution-specific and should not be compared
directly between these families.

[`ZOANVASIM`](https://jmazucheli.github.io/vasicekreg/reference/ZOANVASIM.md):
Zero-and-one-augmented normal-kernel mean family. Here \\\nu=P(Y=0)\\,
\\\tau=P(Y=1\mid Y\>0)\\, and \\\mu=E(Y\mid 0\<Y\<1)\\. Consequently,
\\P(Y=1)=(1-\nu)\tau\\ and \\E(Y)=(1-\nu)\[\tau+(1-\tau)\mu\]\\.

The distribution functions `dNVASIM`, `pNVASIM`, `qNVASIM`, `dNVASIQ`,
`pNVASIQ`, `qNVASIQ`, `dLVASIQ`, `pLVASIQ`, `qLVASIQ`, `dHVASIQ`,
`pHVASIQ`, and `qHVASIQ` call compiled C++ routines through Rcpp. The
boundary-augmented distribution functions are implemented in R and reuse
the compiled `NVASIM` functions for their continuous component.
Parameter validation, the GAMLSS family definitions, and all
log-likelihood derivatives are implemented in R. The mean and variance
components of the
[`LVASIQ()`](https://jmazucheli.github.io/vasicekreg/reference/LVASIQ.md)
and
[`HVASIQ()`](https://jmazucheli.github.io/vasicekreg/reference/HVASIQ.md)
family objects are obtained by numerical quadrature because these
moments have no closed-form expressions.

For the distribution functions and GAMLSS constructors associated with
`NVASIQ`, `LVASIQ`, and `HVASIQ`, the fixed quantile level is supplied
through the `quantile` argument. It is stored in the family definition
and embedded in the residual expression, so no global variable is
required. This fixed quantile level is distinct from the parameter `tau`
in
[`ZOANVASIM()`](https://jmazucheli.github.io/vasicekreg/reference/ZOANVASIM.md),
which represents the conditional probability at one among nonzero
observations.

## References

Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals.
*Journal of Computational and Graphical Statistics*, **5**(3), 236–244.
[doi:10.1080/10618600.1996.10474708](https://doi.org/10.1080/10618600.1996.10474708)

Fischer, M. J., Hui, A. and Hösle, S. (2017). wHS-type distributions
with application to finance. *Journal of Statistics and Management
Systems*, **20**(1), 67–89.
[doi:10.1080/09720510.2016.1190575](https://doi.org/10.1080/09720510.2016.1190575)

Mazucheli, J., Alves, B., Korkmaz, M. Ç., and Leiva, V. (2022). Vasicek
quantile and mean regression models for bounded data: New formulation,
mathematical derivations, and numerical applications. *Mathematics*,
**10**, 1389.
[doi:10.3390/math10091389](https://doi.org/10.3390/math10091389)

Moral, R. A., Hinde, J. and Demetrio, C. G. B. (2017). Half-normal plots
and overdispersed models in R: The hnp package. *Journal of Statistical
Software*, **81**(10), 1–23.
[doi:10.18637/jss.v081.i10](https://doi.org/10.18637/jss.v081.i10)

Ospina, R. and Ferrari, S. L. P. (2010). Inflated beta distributions.
*Statistical Papers*, **51**, 111–126.
[doi:10.1007/s00362-008-0125-4](https://doi.org/10.1007/s00362-008-0125-4)

Ospina, R. and Ferrari, S. L. P. (2012). A general class of zero-or-one
inflated beta regression models. *Computational Statistics & Data
Analysis*, **56**(6), 1609–1623.
[doi:10.1016/j.csda.2011.10.005](https://doi.org/10.1016/j.csda.2011.10.005)

Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive
models for location, scale and shape. *Applied Statistics*, **54**(3),
507–554.
[doi:10.1111/j.1467-9876.2005.00510.x](https://doi.org/10.1111/j.1467-9876.2005.00510.x)

Vasicek, O. A. (2002). The distribution of loan portfolio value. *Risk*,
**15**(12), 160–162.

Witzany, J. (2013). A note on the Vasicek's model with the logistic
distribution. *Ekonomický časopis (Journal of Economics)*, **61**(10),
1053–1066.

Zhao, Y., Lee, A. H., Yau, K. K. W. and McLachlan, G. J. (2011).
Assessing the adequacy of Weibull survival models: A simulated envelope
approach. *Journal of Applied Statistics*, **38**(10), 2089–2097.
[doi:10.1080/02664763.2010.545115](https://doi.org/10.1080/02664763.2010.545115)

## Author

Josmar Mazucheli <jmazucheli@gmail.com>

Bruna Alves <pg402900@uem.br>

## Examples

``` r
# \donttest{
if (requireNamespace("betareg", quietly = TRUE)) {
    data("ReadingSkills", package = "betareg")
    ReadingSkills$dyslexia <- stats::relevel(
        factor(ReadingSkills$dyslexia), ref = "no"
    )

    control <- gamlss::gamlss.control(n.cyc = 200, trace = FALSE)

    ## In both models, mu = E(Y | Y < 1) and nu = P(Y = 1).
    fit_oanvasim <- gamlss::gamlss(
        accuracy1 ~ dyslexia * iq,
        sigma.formula = ~ dyslexia + iq,
        nu.formula = ~ 1,
        family = OANVASIM(),
        data = ReadingSkills,
        control = control
    )

    fit_beoi <- gamlss::gamlss(
        accuracy1 ~ dyslexia * iq,
        sigma.formula = ~ dyslexia + iq,
        nu.formula = ~ 1,
        family = gamlss.dist::BEOI(),
        data = ReadingSkills,
        control = control
    )

    n <- nrow(ReadingSkills)
    comparison <- data.frame(
        family = c("OANVASIM", "BEOI"),
        logLik = -c(
            fit_oanvasim$G.deviance,
            fit_beoi$G.deviance
        ) / 2,
        AIC = c(
            gamlss::GAIC(fit_oanvasim, k = 2),
            gamlss::GAIC(fit_beoi, k = 2)
        ),
        BIC = c(
            gamlss::GAIC(fit_oanvasim, k = log(n)),
            gamlss::GAIC(fit_beoi, k = log(n))
        )
    )
    comparison
}
#>     family   logLik      AIC      BIC
#> 1 OANVASIM 10.46635 -4.93271 9.340807
#> 2     BEOI 10.95973 -5.91945 8.354067
# }
```
