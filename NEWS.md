# vasicekreg 1.2.0 (2026-09-03)

## Interface and documentation review

- Updated the package title and description to cover distribution functions,
  Vasicek-type models, and the hyperbolic-secant kernel and `HSVASIQ` family.
- Standardized boundary terminology as "augmented" and kernel terminology as
  "Vasicek-type" throughout the current documentation.
- Corrected the `bodyfat` documentation: it contains 10 variables, and its five
  response variables are proportions in `(0, 1)`, not values on a 0--100 scale.
- Replaced generic quantile notation `tau` by `p` in the `NVASIM` quantile
  formula and corrected the description of the bivariate-normal CDF.
- Added a real-data `ZOANVASIM` example based on the boundary-valued
  inappropriate-stay proportion in `aep`.
- Extended the numerical derivative checks for all three quantile families
  to a grid of response, quantile, shape, and fixed-quantile values.

## Hyperbolic-secant-kernel quantile regression

- Added the `HSVASIQ()` GAMLSS family for conditional quantile regression
  with a hyperbolic-secant kernel and fixed global quantile level `tau`.
- Added compiled `dHSVASIQ()`, `pHSVASIQ()`, `qHSVASIQ()`, and
  `rHSVASIQ()` functions, including stable log-probability and tail
  calculations.
- Implemented closed-form first, second, and cross log-likelihood derivatives
  with respect to `mu` and `sigma`.
- Added randomized quantile residual support and numerical conditional-moment
  functions for the GAMLSS family.
- Extended `vasicek_envelope()` to fitted `HSVASIQ` models.
- Corrected parametric-envelope simulation for quantile families so the
  global fixed quantile level `tau` is passed to their random generators.
- Added documentation and tests for distribution identities, normalization,
  fixed-quantile preservation, vectorization, analytical derivatives,
  conditional moments, random generation, and GAMLSS fitting.

## Boundary-augmented mean regression

- Added the GAMLSS-conventional aliases `d/p/q/rZANVASIM()`,
  `d/p/q/rOANVASIM()`, and `d/p/q/rZOANVASIM()`. The established numbered
  names remain available. This correspondence allows likelihood-based
  methods such as `vcov.gamlss()` to locate each family's density function.
- Added `ZANVASIM()`, a zero-augmented normal-kernel Vasicek mean family,
  with `nu` modeling the probability at zero.
- Added `OANVASIM()`, a one-augmented normal-kernel Vasicek mean family,
  with `nu` modeling the probability at one.
- Added `ZOANVASIM()`, a zero-and-one-augmented normal-kernel Vasicek mean
  family. In this family, `nu` is the probability at zero and `tau` is
  the conditional probability at one among nonzero observations.
- Added the corresponding `d0NVASIM()`, `p0NVASIM()`, `q0NVASIM()`,
  `r0NVASIM()`, `d1NVASIM()`, `p1NVASIM()`, `q1NVASIM()`,
  `r1NVASIM()`, `d01NVASIM()`, `p01NVASIM()`, `q01NVASIM()`, and
  `r01NVASIM()` functions.
- Added analytical derivatives, randomized quantile residuals, marginal
  moments, and GAMLSS fitting support for all three boundary-augmented
  families.

## Simulated residual envelopes

- Added `vasicek_envelope()` for parametric-bootstrap simulated envelopes
  based on fitted `gamlss` objects.
- Added support for normalized randomized quantile residuals and generalized
  Cox--Snell residuals, with the latter defined from the same fitted
  probability integral transform and evaluated on the log-survival scale for
  numerical stability in the upper tail.
- Added pointwise percentile and pointwise minimum--maximum envelopes.
- Each accepted bootstrap sample is simulated from the fitted model, refitted,
  and used to recalculate and order the requested residuals.
- Added print and plot methods for `vasicek_envelope` objects. The
  documentation now distinguishes these full quantile--quantile plots from
  half-normal plots and explains the finite-sample simulated mean curve.
- Added a row-count check to automatic bootstrap refitting so data altered by
  missing-value removal or `subset` must be supplied explicitly or handled
  by a custom refit function.

## Documentation and tests

- Standardized the R source filenames as `dpqr-0NvasicekM.R`,
  `dpqr-1NvasicekM.R`, and `dpqr-01NvasicekM.R`.
- Added numerical verification of all first, second, and cross derivatives
  for `ZANVASIM`, `OANVASIM`, and `ZOANVASIM`.
- Added tests of the randomized quantile residual expressions for the
  zero-augmented, one-augmented, and zero-and-one-augmented families, including
  the CDF jump intervals, interior observations, reproducibility, and
  simulation from the fitted family.
- Added tests for simulated residual envelopes, including argument
  validation, reproducibility, returned components, plotting, stable
  Cox--Snell upper-tail calculations, and incompatible refit data.
- Clarified the roles of `x` and `q` and distinguished the boundary
  behavior of the zero-augmented, one-augmented, and zero-and-one-augmented
  distributions.
- Updated the package overview to describe the boundary-augmented families
  and to distinguish the `tau` parameter in `ZOANVASIM` from the fixed
  quantile level used by `NVASIQ` and `LVASIQ`.
- Added a package-level example comparing `OANVASIM` and `BEOI` using the
  one-inflated `accuracy1` response from the `ReadingSkills` data.
- Expanded the README references for the Vasicek distributions, GAMLSS,
  boundary-augmented models, randomized quantile and Cox--Snell residuals,
  simulated envelopes, and Rcpp.

# vasicekreg 1.1.0 (2026-07-28)

## Interface reorganization
- Renamed the normal-kernel mean family and distribution functions from
  `VASIM` to `NVASIM`.
- Renamed the normal-kernel quantile family and distribution functions from
  `VASIQ` to `NVASIQ`.
- Retained `LVASIQ` for the logistic-kernel quantile distribution functions.
- Renamed the corresponding R and C++ source files to identify the kernel
  (`N` or `L`) and parameterization (`M` or `Q`) explicitly.
- Renamed the internal Rcpp routines and regenerated their registration
  interfaces.
- Updated the namespace, documentation, examples, package overview, and
  tests for the new interface.

## GAMLSS quantile level
- `NVASIQ()` does not accept `tau` as an argument. For GAMLSS fitting,
  `tau` must be defined as a scalar variable in the global environment.
- The distribution functions `dNVASIQ()`, `pNVASIQ()`, `qNVASIQ()`, and
  `rNVASIQ()` continue to accept `tau` explicitly.

## Logistic-kernel quantile regression
- Added the `LVASIQ()` GAMLSS family for conditional quantile regression
  with the logistic-kernel Vasicek distribution.
- Implemented closed-form first, second, and cross derivatives with respect
  to `mu` and `sigma`; model fitting does not use numerical differentiation.
- Added quantile residual support and numerical moment functions for the
  family object.

---

# vasicekreg 1.0.3 (2026-07-28)

## Bug fixes
- Restored the missing implementation of `qVASIM()`.
- Corrected `log.p` handling in all quantile functions.
- Corrected the bivariate-normal calculation used by the `variance`
  components of `VASIM()` and `VASIQ()`.
- Made the fixed quantile level an explicit argument of `VASIQ()` and
  removed dependence on a global `tau` object.
- Corrected the second derivative with respect to `mu` in `VASIQ()`.
- Added stable upper-tail calculations and consistent parameter-length
  checks to the C++ routines.

## Documentation and tests
- Clarified that the logistic-kernel implementation currently provides
  distribution functions but not a GAMLSS family.
- Added executable tests for distribution identities, moments, derivatives,
  vectorization, and fixed-quantile handling.

---

# vasicekreg 1.0.2 (2026-01-12)

## Improvements
- Performance optimizations in `dpqr-vasicekmean.R` and
  `dpqr-vasicekquant.R`.

---

# vasicekreg 1.0.1 (2021-12-05)

## Bug fixes
- Corrected the term `variance = function(mu, sigma)` in the VASIM family.
- Corrected the term `variance = function(mu, sigma)` in the VASIQ family.

---

# vasicekreg 1.0.0 (2021-07-05)

## Initial release
- First release of the vasicekreg package on CRAN.
