# Simulated envelopes for residuals from Vasicek GAMLSS models

Constructs pointwise simulated envelopes for normalized randomized
quantile residuals and generalized Cox–Snell residuals from a Vasicek
model fitted by
[`gamlss()`](https://rdrr.io/pkg/gamlss/man/gamlss.html). Each simulated
response is generated under the fitted model, the model is re-estimated,
and the residuals from the re-estimated model are ordered before the
envelope is computed.

## Usage

``` r
vasicek_envelope(
  object,
  residual = c("quantile", "cox-snell"),
  nsim = 200L,
  level = 0.95,
  envelope = c("quantile", "minmax"),
  seed = NULL,
  data = NULL,
  max.attempts = 5L * nsim,
  refit = NULL,
  simulate = NULL,
  verbose = interactive()
)

# S3 method for class 'vasicek_envelope'
print(x, ...)

# S3 method for class 'vasicek_envelope'
plot(
  x,
  which = x$residual[1L],
  show.mean = TRUE,
  xlab = NULL,
  ylab = NULL,
  main = "",
  envelope.col = "grey85",
  mean.col = "blue",
  reference.col = "red",
  point.col = "black",
  ...
)
```

## Arguments

- object:

  A fitted object of class `"gamlss"` using one of the Vasicek families
  supplied by vasicekreg.

- residual:

  Character vector selecting `"quantile"`, `"cox-snell"`, or both.

- nsim:

  Number of successfully re-fitted simulated samples.

- level:

  Pointwise coverage probability used when `envelope = "quantile"`.

- envelope:

  Either `"quantile"`, for percentile limits, or `"minmax"`, for the
  minimum and maximum at each order position.

- seed:

  Optional integer seed. The previous random-number state is restored on
  exit.

- data:

  Optional data frame used in the original fit. It is normally recovered
  from `object$call$data`.

- max.attempts:

  Maximum number of simulation and re-fitting attempts. The default is
  five times `nsim`.

- refit:

  Optional function with arguments `object`, `y`, and `data`. It must
  return a re-fitted `"gamlss"` object. This is useful for transformed
  responses or nonstandard fitting calls.

- simulate:

  Optional function with the single argument `object` that returns one
  simulated response vector. By default, the random generator associated
  with the fitted Vasicek family is used.

- verbose:

  Logical; if `TRUE`, reports progress and failed fits.

- x:

  An object returned by `vasicek_envelope()`.

- ...:

  Further arguments. For the plot method, they are passed to
  [`points()`](https://rdrr.io/r/graphics/points.html).

- which:

  Residual type to be plotted.

- show.mean:

  Logical; if `TRUE`, draws the pointwise mean of the ordered simulated
  residuals.

- xlab, ylab, main:

  Graphical labels.

- envelope.col, mean.col, reference.col, point.col:

  Graphical colors.

## Value

An object of class `"vasicek_envelope"`. Its `results` component
contains, for each requested residual, the theoretical order statistics,
ordered observed residuals, pointwise lower, mean and upper curves, and
the matrix of ordered simulated residuals. The object also records the
number of successful simulations, attempts and failures.

## Details

Let \\U_i\\ denote the probability integral transform used by the
normalized randomized quantile residual. The two residuals are
\$\$r_i^q=\Phi^{-1}(U_i)\$\$ and \$\$r_i^{CS}=-\log(1-U_i).\$\$ For
augmented families, \\U_i\\ is randomized over the relevant jump of the
fitted distribution at zero or one, so the reported residuals inherit
whatever randomization
[`residuals()`](https://rdrr.io/r/stats/residuals.html) applies for
those families; the Cox–Snell residual is then also randomized at
boundary observations and has an \\\operatorname{Exp}(1)\\ reference
distribution under a correctly specified model. The Cox–Snell residual
is obtained from the quantile residual on the log survival scale,
\\-\log\\\Pr(Z\>r_i^q)\\\\, which is numerically stable in the upper
tail and preserves the ordering of the quantile residuals.

Each residual is displayed as a full quantile–quantile plot of the
ordered residuals against the corresponding theoretical quantiles
(normal or exponential), not as a half-normal plot. The envelope is
pointwise, not simultaneous: even under a correctly specified model a
fraction of points is expected to fall outside the band. With
`envelope = "quantile"`, the two tail probabilities are equal and sum to
\\1-\mathrm{level}\\. With `envelope = "minmax"`, the limits are the
minimum and maximum at each order position, following the construction
used by Zhao et al.

The pointwise mean of the ordered simulated residuals is the reference
calibrated to the fitted model and to the sample size. The identity line
drawn by the plot method is a theoretical idealization and, for
Cox–Snell residuals, may separate from the mean curve in the upper tail
at finite \\n\\; in that region the mean curve is the more reliable
reference.

Failed or nonconverged re-fits are discarded and replaced until `nsim`
successful samples are obtained or `max.attempts` is reached. The
simulation uses fitted values for every distribution parameter, so
covariate-dependent parameters are retained. Automatic simulation and
refitting of `NVASIQ`, `LVASIQ`, and `HVASIQ` recover the fixed level
embedded in the fitted object; no global quantile-level variable is
consulted. Consequently, fitted models at different quantile levels can
be used in the same R session safely.

Automatic refitting assumes that `data` contains exactly the
observations used by the original fit; if rows were dropped for missing
values or via `subset`, supply a matching `data` or a custom `refit`.

## References

Moral, R. A., Hinde, J., and Demetrio, C. G. B. (2017). Half-normal
plots and overdispersed models in R: The hnp package. *Journal of
Statistical Software*, **81**(10), 1–23.
[doi:10.18637/jss.v081.i10](https://doi.org/10.18637/jss.v081.i10)

Zhao, Y., Lee, A. H., Yau, K. K. W., and McLachlan, G. J. (2011).
Assessing the adequacy of Weibull survival models: A simulated envelope
approach. *Journal of Applied Statistics*, **38**, 2089–2097.
[doi:10.1080/02664763.2010.545115](https://doi.org/10.1080/02664763.2010.545115)

## Examples

``` r
if (FALSE) { # \dontrun{
library(gamlss)
set.seed(123)
dat <- data.frame(y = rNVASIM(100, mu = 0.55, sigma = 0.30))
fit <- gamlss(y ~ 1, sigma.formula = ~ 1, family = NVASIM(),
              data = dat, trace = FALSE)

env <- vasicek_envelope(fit, nsim = 200, seed = 123)
plot(env, which = "quantile")
plot(env, which = "cox-snell")
} # }
```
