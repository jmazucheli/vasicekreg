# Body Fat Proportions Dataset

Body fat proportions from individuals assisted in a public hospital in
Curitiba, Paraná, Brazil.

## Usage

``` r
bodyfat
```

## Format

A data frame with 298 observations and 10 variables. The five body-fat
responses are proportions in \\(0,1)\\ (for example, 0.163 represents
16.3 percent):

- `ID`: individual identifier.

- `ARMS`: arms fat proportion.

- `LEGS`: legs fat proportion.

- `BODY`: body fat proportion.

- `ANDROID`: android fat proportion.

- `GYNECOID`: gynoid fat proportion.

- `AGE`: age of individuals.

- `BMI`: body mass index.

- `SEX`: 1 for female and 2 for male.

- `IPAQ`: physical activity level according to IPAQ (0 = sedentary, 1 =
  insufficiently active, 2 = active).

## References

Mazucheli, J., Alves, B., Korkmaz, M. Ç., and Leiva, V. (2022). Vasicek
quantile and mean regression models for bounded data: New formulation,
mathematical derivations, and numerical applications. *Mathematics*,
**10**, 1389.

Mazucheli, J., Leiva, V., Alves, B., and Menezes, A. F. B. (2021). A new
quantile regression for modeling bounded data under a unit
Birnbaum-Saunders distribution with applications in medicine and
politics. *Symmetry*, **13**(4), 1–21.

Petterle, R. R., Bonat, W. H., Scarpin, C. T., Jonasson, T., and Borba,
V. Z. C. (2020). Multivariate quasi-beta regression models for
continuous bounded data. *The International Journal of Biostatistics*,
**17**(1), 39–53.

## Author

Josmar Mazucheli <jmazucheli@gmail.com>

Bruna Alves <pg402900@uem.br>

## Examples

``` r
data(bodyfat, package = "vasicekreg")

bodyfat$AGE <- bodyfat$AGE - 46.00
bodyfat$BMI <- bodyfat$BMI - 24.72
bodyfat$SEX <- as.factor(bodyfat$SEX)
bodyfat$IPAQ<- as.factor(bodyfat$IPAQ)

library(gamlss)

## Mean regression model
fitmean <- gamlss(
  ARMS ~ AGE + BMI + SEX + IPAQ,
  data = bodyfat,
  family = NVASIM(mu.link = "logit", sigma.link = "logit")
)
#> GAMLSS-RS iteration 1: Global Deviance = -667.3348 
#> GAMLSS-RS iteration 2: Global Deviance = -908.7608 
#> GAMLSS-RS iteration 3: Global Deviance = -911.2209 
#> GAMLSS-RS iteration 4: Global Deviance = -911.2213 

if (FALSE) { # \dontrun{
quantile_levels <- c(0.25, 0.50, 0.75)

## Quantile regression models with the normal kernel
fit_normal <- lapply(quantile_levels, function(level) {
  gamlss(
    ARMS ~ AGE + BMI + SEX + IPAQ,
    data = bodyfat,
    family = NVASIQ(
      quantile = level,
      mu.link = "logit",
      sigma.link = "logit"
    )
  )
})

## Quantile regression models with the logistic kernel
fit_logistic <- lapply(quantile_levels, function(level) {
  gamlss(
    ARMS ~ AGE + BMI + SEX + IPAQ,
    data = bodyfat,
    family = LVASIQ(
      quantile = level,
      mu.link = "logit",
      sigma.link = "logit"
    )
  )
})

## Quantile regression models with the Hyperbolic-secant-kernel
fit_hsk <- lapply(quantile_levels, function(level) {
  gamlss(
    ARMS ~ AGE + BMI + SEX + IPAQ,
    data = bodyfat,
    family = HVASIQ(
      quantile = level,
      mu.link = "logit",
      sigma.link = "logit"
    )
  )
})

lapply(fit_normal, summary)
lapply(fit_logistic, summary)
lapply(fit_hsk, summary)
} # }
```
