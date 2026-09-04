# Transportation to campus

Data from a stratified sample of 60 respondents concerning their mode of
transportation to campus. The sampling design oversampled respondents
who sometimes traveled to campus by bicycle.

## Format

A data frame with 60 observations and 7 variables:

- ntrips:

  Number of trips to campus during the preceding four weeks.

- nbiked:

  Number of those trips made by bicycle.

- status:

  Respondent's institutional status: faculty, staff, or student.

- gender:

  Respondent's gender, recorded as `"F"` or `"M"`.

- parking:

  Duration of the parking permit, in months.

- distance:

  Distance to campus. The source documentation does not specify the
  measurement unit.

- propbiked:

  Proportion of trips to campus made by bicycle, calculated as
  `nbiked / ntrips`.

## Source

Korosteleva, O. (2019). *Advanced Regression Models with SAS and R*.
Boca Raton, FL: CRC Press.

<https://github.com/AndrMenezes/uwquantreg>

## Details

The data originate from a consulting study reported by Korosteleva
(2019). Menezes, Mazucheli, and Bourguignon (2021) analyzed them using
an inflated unit-Weibull quantile regression model. The object
distributed here retains the variable names and values supplied by the
uwquantreg package.

## References

Menezes, A. F. B., Mazucheli, J. and Bourguignon, M. (2021). A
parametric quantile regression approach for modeling zero- or
one-inflated double bounded data. *Biometrical Journal*, **63**(4),
841–858.
[doi:10.1002/bimj.202000126](https://doi.org/10.1002/bimj.202000126)

Menezes, A. F. B. (2026). uwquantreg: unit-Weibull quantile regression.
R package version 0.1.0. <https://github.com/AndrMenezes/uwquantreg>

## Examples

``` r
data("transport", package = "vasicekreg")
str(transport)
#> 'data.frame':    60 obs. of  7 variables:
#>  $ ntrips   : int  26 13 17 15 17 20 17 14 26 8 ...
#>  $ nbiked   : int  6 0 0 0 0 13 0 0 7 5 ...
#>  $ status   : chr  "student" "faculty" "faculty" "faculty" ...
#>  $ gender   : chr  "F" "M" "M" "M" ...
#>  $ parking  : int  6 9 9 9 6 6 9 9 6 12 ...
#>  $ distance : int  7 15 31 9 34 3 17 10 3 3 ...
#>  $ propbiked: num  0.231 0 0 0 0 ...
with(transport, all.equal(propbiked, nbiked / ntrips))
#> [1] TRUE
```
