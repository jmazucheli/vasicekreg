---
title: "vasicekreg: An R package for Vasicek-based mean and quantile regression models for bounded data"
tags:
  - R
  - bounded data
  - Vasicek distribution
  - beta regression
  - quantile regression
  - parametric regression
authors:
  - name: Josmar Mazucheli
    orcid: 0000-0002-5405-6644
    affiliation: 1
affiliations:
  - name: Department of Statistics, Universidade Estadual de Maringá, Brazil
    index: 1
date: "2026-01-13"
bibliography: paper.bib
---

## Summary

Bounded continuous outcomes on the unit interval arise across many fields, including biomedicine, epidemiology, environmental sciences, and finance. Beta regression has become the default modeling tool for such data, largely due to its flexibility. However, beta-based approaches are primarily mean-oriented and can be less convenient when the scientific question targets conditional quantiles (e.g., tail behavior), which often requires nontrivial reparameterizations and may lead to numerical instability.

The Vasicek distribution provides a competitive alternative for modeling unit-interval data. Beyond supporting regression on the mean, it naturally enables quantile-based parameterizations, making fully parametric quantile regression feasible in a coherent likelihood framework. Despite these advantages, Vasicek-based regression models have been comparatively underexplored in the statistical literature.

The `vasicekreg` R package implements likelihood-based inference for both mean and quantile regression models using the Vasicek distribution, providing a unified and reproducible workflow for modeling bounded data.

## Statement of need

Regression modeling for bounded responses is typically carried out using the beta distribution and its variants. While effective in many applications, beta regression does not naturally yield parametric quantile regression, and extensions often rely on complex parameterizations or ad hoc strategies that can complicate interpretation and computation.

The Vasicek distribution is a two-parameter family on the unit interval originally arising in credit risk modeling, where it describes default or loss rates in large homogeneous portfolios [@VasicekDistribution; @VasicekQuantilesCI]. Its structure admits an analytically tractable quantile function and supports quantile-based parameterizations in a particularly direct way, which is advantageous when researchers need conditional quantiles and not only conditional means.

Even though the Vasicek distribution has been studied and compared with other unit-interval distributions [@BeyondBetaVasicek], the regression literature has used it far less than beta-based models. To the best of our knowledge, the first comprehensive development of *both* mean regression and quantile regression models based on the Vasicek distribution was presented in a peer-reviewed article published in *Mathematics* in 2022 [@Mazucheli2022Mathematics]. The `vasicekreg` package operationalizes these methodological developments and makes them accessible to applied researchers through a stable and documented implementation in R.

## Software description

The `vasicekreg` package provides a cohesive interface to fit Vasicek-based regression models for responses in (0, 1). The software includes:

- Maximum likelihood estimation for mean and quantile regression parameterizations;
- Flexible regression predictors and link functions for model components;
- Model-based inference and comparison via likelihood tools;
- Predictive routines and utilities supporting reproducible analysis workflows.

The package is distributed through CRAN and is designed to integrate seamlessly with standard R workflows for statistical modeling.

## Implementation details

Estimation in `vasicekreg` is performed by maximizing the Vasicek log-likelihood under regression structures. The package supports both mean-based and quantile-based parameterizations, allowing the analyst to choose the modeling target most aligned with the scientific question (central tendency versus tails). Particular attention is given to numerical robustness and stable optimization behavior, with implementation choices guided by practical experience in simulation and real-data analyses.

## Research impact

The Vasicek distribution has a long tradition in credit risk and portfolio loss modeling, including inferential work for quantiles and risk measures [@VasicekDistribution; @VasicekQuantilesCI]. More recently, comparative studies have evaluated Vasicek alongside beta and other unit-interval distributions, highlighting its flexibility but also the relative scarcity of accessible implementations for modern regression tasks [@BeyondBetaVasicek]. Related model developments also include logistic-type formulations and extensions in risk modeling contexts [@VasicekLogisticNote].

The `vasicekreg` package has been used in peer-reviewed research in both methodological and applied settings. It supported the first systematic investigation of Vasicek-based mean and quantile regression models for bounded data published in *Mathematics* [@Mazucheli2022Mathematics], and it has been employed in a broader review and applied study published in *Computer Methods and Programs in Biomedicine* with biomedical and COVID-19 data applications [@Mazucheli2022CMPB].

From its first CRAN release on 2021-05-01 to the latest CRAN submission on 2026-01-12, `vasicekreg` has accumulated 12,321 downloads (as measured via `dlstats::cran_stats("vasicekreg")`), providing quantitative evidence of sustained community adoption.

## Availability

The package is available from CRAN at:  
<https://cran.r-project.org/package=vasicekreg>

The source code is publicly available in a version-controlled repository and distributed under an OSI-approved open-source license.

## Acknowledgements

The author acknowledges collaborators and users for feedback and testing, and the open-source R community for foundational tools enabling reproducible statistical software.
