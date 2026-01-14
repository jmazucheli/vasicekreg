---
title: 'vasicekreg package: Mean and quantile regression models for bounded data based
  on the Vasicek distribution in R'
tags:
- R
- bounded data
- Vasicek distribution
- beta regression
- quantile regression
date: "2026-01-13"
output: pdf_document
affiliations:
- name: Department of Statistics, Universidade Estadual de Maringá, Brazil
  index: 1
bibliography: paper.bib
authors:
- name: Josmar Mazucheli
  orcid: "0000-0002-5405-6644"
  affiliation: 1
- name: Bruna Alves
  affiliation: 1
---

## Summary

Continuous response variables restricted to the unit interval are frequently observed in 
areas such as biomedicine, epidemiology, environmental sciences, and finance. Regression models based on 
the beta distribution are commonly employed to analyze such data. However, beta regression is primarily used 
to assess the effect of one or more explanatory variables on the mean of the response variable. In contrast, 
the Vasicek distribution offers a flexible alternative that extends beyond mean-based analysis by admitting 
direct quantile-based parameterizations. Consequently, the \texttt{vasicekreg} R package implements likelihood-based 
inference for both mean and quantile regression models, providing a unified framework for modeling 
bounded data. The package is intended for applied statisticians and researchers working with rates, proportions, and
bounded responses in empirical studies.

## Statement of need

Regression models for bounded responses are most commonly based on the beta 
distribution [@FerrariCribari2004; @SmithsonVerkuilen2006; @CribariNetoZeileis2010]. 
Although beta regression has become the standard approach for modeling rates and proportions,
it is primarily designed for conditional mean inference and does not naturally support
fully parametric quantile regression.

The Vasicek distribution is a two-parameter model defined on the unit interval, 
originally introduced in the context of credit risk modeling 
[@VasicekDistribution; @VasicekQuantilesCI]. It admits a closed-form quantile 
function and allows for direct quantile-based parameterizations, making it a 
natural alternative to the beta distribution when interest lies in modeling 
conditional quantiles in addition to conditional means. Despite these characteristics, 
Vasicek-based regression models have received considerably less attention in the statistical 
literature than beta-based approaches [@BeyondBetaVasicek; @VasicekLogisticNote].

Compared to existing beta regression tools, such as the \texttt{betareg} 
package [@CribariNetoZeileis2010], the \texttt{vasicekreg} package provides a fully 
parametric likelihood-based framework for both mean and quantile regression 
of bounded responses. To the best of our knowledge, the first comprehensive 
treatment of Vasicek-based mean and quantile regression models was presented 
in @Mazucheli2022Mathematics. The \texttt{vasicekreg} package implements these methodological 
developments by offering two 
complete regression frameworks within the \texttt{gamlss} environment [@RigbyStasinopoulos2005]: one for 
Vasicek mean regression and another for Vasicek quantile regression 
at any percentile order $\tau \in (0,1)$.




## Software description

The \texttt{vasicekreg} package provides a cohesive interface for fitting Vasicek-based regression models to 
responses bounded on the unit interval. The software supports:

- likelihood-based estimation for mean and quantile regression parameterizations;
- flexible regression predictors and link functions;
- model-based inference, prediction, and simulation tools.

The package is distributed through CRAN and integrates seamlessly with standard R modeling workflows.

## Implementation details

Estimation in \texttt{vasicekreg} is implemented via two custom distribution families within the `gamlss` framework, enabling both mean-based and quantile-based Vasicek regression models. The quantile regression framework allows inference at any fixed percentile order $\tau \in (0,1)$. Core computational routines, including evaluation of the density, distribution, quantile function, and random number generation, are implemented in C++ using the `Rcpp` interface, ensuring computational efficiency and numerical robustness.

## Research impact

The Vasicek distribution has been extensively studied in the context of credit portfolio losses and risk measures [@VasicekDistribution; @VasicekQuantilesCI], and has been compared with beta and other unit-interval distributions in simulation studies [@BeyondBetaVasicek]. Related methodological developments include logistic-type formulations and extensions in risk modeling [@VasicekLogisticNote].

The \texttt{vasicekreg} package has supported peer-reviewed methodological and applied studies
by providing likelihood-based estimation for Vasicek mean and quantile regression models. It 
provided the computational basis for the first systematic investigation of Vasicek-based mean and quantile regression 
models for bounded data published in *Mathematics* [@Mazucheli2022Mathematics], and was subsequently
employed in a 
review and application study involving biomedical and COVID-19 data [@Mazucheli2022CMPB].

Since its first CRAN release on 2021-05-01, and with the latest submission on 2026-01-12, the package has accumulated 12,321 downloads, indicating sustained adoption by the statistical modeling community.

## Availability

The package is available from CRAN at  
<https://cran.r-project.org/package=vasicekreg>.

The source code is publicly available in a version-controlled repository and is distributed under an OSI-approved open-source license.

## Acknowledgements

The author thanks collaborators and users for valuable feedback and testing, and acknowledges the R open-source community for foundational tools enabling this work.

## AI usage disclosure

No generative artificial intelligence tools were used in the development of the software or in the writing of this manuscript.

## References

