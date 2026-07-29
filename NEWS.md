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
