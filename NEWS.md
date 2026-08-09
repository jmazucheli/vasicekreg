# vasicekreg development version (2026-08-09)

## Boundary-adjusted mean regression

- Added the GAMLSS-conventional aliases `d/p/q/rZANVASIM()`,
  `d/p/q/rOANVASIM()`, and `d/p/q/rZOANVASIM()`. The established numbered
  names remain available. This correspondence allows likelihood-based
  methods such as `vcov.gamlss()` to locate each family's density function.
- Added `ZANVASIM()`, a zero-adjusted normal-kernel Vasicek mean family,
  with `nu` modeling the probability at zero.
- Added `OANVASIM()`, a one-adjusted normal-kernel Vasicek mean family,
  with `nu` modeling the probability at one.
- Added `ZOANVASIM()`, a zero-and-one-adjusted normal-kernel Vasicek mean
  family. In this family, `nu` is the probability at zero and `tau` is
  the conditional probability at one among nonzero observations.
- Added the corresponding `d0NVASIM()`, `p0NVASIM()`, `q0NVASIM()`,
  `r0NVASIM()`, `d1NVASIM()`, `p1NVASIM()`, `q1NVASIM()`,
  `r1NVASIM()`, `d01NVASIM()`, `p01NVASIM()`, `q01NVASIM()`, and
  `r01NVASIM()` functions.
- Added analytical derivatives, randomized quantile residuals, marginal
  moments, and GAMLSS fitting support for all three boundary-adjusted
  families.

## Documentation and tests

- Standardized the R source filenames as `dpqr-0NvasicekM.R`,
  `dpqr-1NvasicekM.R`, and `dpqr-01NvasicekM.R`.
- Added numerical verification of all first, second, and cross derivatives
  for `ZANVASIM`, `OANVASIM`, and `ZOANVASIM`.
- Clarified the roles of `x` and `q` and distinguished the boundary
  behavior of the zero-adjusted, one-adjusted, and zero-and-one-adjusted
  distributions.
- Updated the package overview to describe the boundary-adjusted families
  and to distinguish the `tau` parameter in `ZOANVASIM` from the fixed
  quantile level used by `NVASIQ` and `LVASIQ`.
- Added a package-level example comparing `OANVASIM` and `BEOI` using the
  one-inflated `accuracy1` response from the `ReadingSkills` data.

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

