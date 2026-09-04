# Survival of young trees

Data from a study conducted by a parks and recreation department on the
two-year survival of young trees planted in 26 parks.

## Format

A data frame with 26 observations and 7 variables:

- planted:

  Number of trees planted.

- survived:

  Number of trees that survived for two years.

- pest:

  Frequency of pest-control treatment.

- fertilization:

  Frequency of soil fertilization.

- precip:

  Average annual precipitation, in inches.

- wind:

  Average annual wind speed, in miles per hour.

- prop:

  Proportion of planted trees that survived for two years, calculated as
  `survived / planted`.

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
data("trees", package = "vasicekreg")
str(trees)
#> 'data.frame':    26 obs. of  7 variables:
#>  $ planted      : int  125 115 250 95 140 75 185 20 110 80 ...
#>  $ survived     : int  125 68 101 85 48 75 163 9 83 80 ...
#>  $ pest         : int  3 0 1 2 3 3 3 3 3 0 ...
#>  $ fertilization: int  1 0 1 2 1 2 3 0 1 1 ...
#>  $ precip       : int  18 8 17 22 15 27 15 18 24 18 ...
#>  $ wind         : num  9.6 13.4 12.8 10 15.1 6.3 12.3 9.4 13.1 7.8 ...
#>  $ prop         : num  1 0.591 0.404 0.895 0.343 ...
with(trees, all.equal(prop, survived / planted))
#> [1] TRUE
```
