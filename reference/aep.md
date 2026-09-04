# Hospital-stay appropriateness data

Data on 1,383 patients admitted to Hospital del Mar, Barcelona, during
1988 and 1990. The data record the total length of stay and the number
of days classified as inappropriate.

## Format

A data frame with 1,383 observations and 8 variables:

- los:

  Total number of days spent in hospital.

- noinap:

  Number of hospital-stay days classified as inappropriate.

- loglos:

  Logarithm of length of stay divided by 10, `log(los / 10)`.

- sex:

  Patient's sex: factor with levels `1` (male) and `2` (female).

- ward:

  Hospital ward: factor with levels `1` (medical), `2` (surgical), and
  `3` (other).

- year:

  Admission year: factor with levels `88` and `90`.

- age:

  Patient's age minus 55 years.

- y:

  Two-column matrix response whose first column is `noinap` and whose
  second column is `los - noinap`.

## Source

Gange, S. J., Munoz, A., Saez, M. and Alonso, J. (1996). Use of the
beta-binomial distribution to model the effect of policy changes on
appropriateness of hospital stays. *Applied Statistics*, **45**(3),
371–382.

## Details

Gange et al. (1996) modeled the number of inappropriate days conditional
on total length of stay using binomial and beta-binomial models. For
augmented continuous-response models, the patient-level proportion can
be constructed as `noinap / los`. The object distributed here retains
the structure and values supplied by the gamlss.data package.

## References

Stasinopoulos, M. and Rigby, R. (2025). gamlss.data: Data for
Generalized Additive Models for Location Scale and Shape. R package
version 6.0-7. <https://CRAN.R-project.org/package=gamlss.data>

## Examples

``` r
data("aep", package = "vasicekreg")
str(aep)
#> 'data.frame':    1383 obs. of  8 variables:
#>  $ los   : num  15 42 8 9 7 10 8 8 21 8 ...
#>  $ noinap: num  0 20 6 6 0 2 6 0 7 0 ...
#>  $ loglos: num  0.405 1.435 -0.223 -0.105 -0.357 ...
#>  $ sex   : Factor w/ 2 levels "1","2": 2 2 1 1 1 2 1 2 1 2 ...
#>  $ ward  : Factor w/ 3 levels "1","2","3": 2 1 1 2 2 2 2 2 1 1 ...
#>  $ year  : Factor w/ 2 levels "88","90": 1 1 1 1 1 1 1 2 2 2 ...
#>  $ age   : num  0 18 19 23 2 -8 15 -15 15 40 ...
#>  $ y     : num [1:1383, 1:2] 0 20 6 6 0 2 6 0 7 0 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:2] "noinap" "failures"
inappropriate <- with(aep, noinap / los)
table(inappropriate == 0, inappropriate == 1)
#>        
#>         FALSE TRUE
#>   FALSE   552   68
#>   TRUE    763    0
```
