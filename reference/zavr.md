# Zero-Augmented Vasicek random-intercept model

Fits a two-component model for responses in \\\[0,1)\\ with a point mass
at zero (presence component) and a Vasicek distribution for positive
values.

## Usage

``` r
zavr(
  data,
  y,
  formula_bin = NULL,
  formula_cont = NULL,
  random = NULL,
  logistic_cov = NULL,
  vasicek_cov = NULL,
  subject_ind = NULL,
  time_ind,
  component_wise_test = TRUE,
  quad_n = 30,
  verbose = FALSE,
  joint_test = NULL,
  sd_lower = 1e-05,
  start = NULL,
  control = list(),
  hessian = TRUE
)
```

## Arguments

- data:

  A data.frame.

- y:

  Character name of the response column.

- formula_bin:

  One-sided formula for the presence component.

- formula_cont:

  One-sided formula for the Vasicek component.

- random:

  Formula `~ 1 | Subject` for the random intercept.

- logistic_cov:

  Deprecated: character vector of presence covariates.

- vasicek_cov:

  Deprecated: character vector of Vasicek covariates.

- subject_ind:

  Deprecated: subject column name.

- time_ind:

  Character name of the time column.

- component_wise_test:

  Logical; compute component-wise LRTs?

- quad_n:

  Integer; number of Gauss-Hermite quadrature points.

- verbose:

  Logical; print progress?

- joint_test:

  NULL, TRUE, or FALSE; compute joint LRTs?

- sd_lower:

  Lower bound for the random-effect SD.

- start:

  Optional list of starting values.

- control:

  List of control parameters.

- hessian:

  Logical; compute Hessian-based standard errors?

## Value

An object of class `"zavr"`.
