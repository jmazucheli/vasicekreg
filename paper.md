---
title: "vasicekreg: An R package for Vasicek-based mean and quantile regression models
  for bounded data"
tags:
- R
- bounded data
- Vasicek distribution
- quantile regression
- parametric regression
date: "2026-01-12"
output: pdf_document
authors:
- name: Josmar Mazucheli
  orcid: "0000-0002-5405-6644"
  affiliation: 1
- name: Bruna Alves
  affiliation: 1
- name: André F. B. Menezes
  affiliation: 1
- name: Víctor Leiva
  affiliation: 2
bibliography: paper.bib
affiliations:
- name: Department of Statistics, Universidade Estadual de Maringá, Brazil
  index: 1
- name: School of Industrial Engineering, Pontificia Universidad Católica de Valparaíso,
    Chile
  index: 2
---

## Summary

Modeling bounded continuous responses on the unit interval is a common problem in applied statistics, particularly in areas such as economics, political science, environmental sciences, and biomedical research. While beta regression and related approaches are widely used, they may present limitations when modeling conditional quantiles or when robust alternatives to mean-based inference are required.

The `vasicekreg` package provides a comprehensive implementation of mean and quantile regression models based on the Vasicek distribution for bounded data. By parameterizing the Vasicek distribution in terms of either its mean or a fixed quantile, the package enables flexible and interpretable regression modeling within a fully parametric likelihood-based framework.

## Statement of need

Existing regression models for bounded data often focus on conditional mean modeling and offer limited support for parametric quantile regression. Distribution-free quantile regression methods, although flexible, may suffer from issues such as quantile crossing and reduced interpretability in parametric settings.

The `vasicekreg` package addresses these limitations by providing a unified framework for modeling both conditional means and arbitrary conditional quantiles using the Vasicek distribution. The package is particularly suitable for applications where interest lies in different regions of the conditional distribution of bounded responses, offering an alternative to beta, Kumaraswamy, and related regression models.

## Software description

The `vasicekreg` package is implemented in the R programming language and provides functions for fitting mean and quantile regression models using maximum likelihood estimation. The software supports flexible link functions for regression parameters, likelihood-based inference, and model comparison.

The package is designed to integrate seamlessly with standard R modeling workflows and emphasizes reproducibility, numerical stability, and extensibility for methodological and applied research.

## Implementation details

Estimation procedures in `vasicekreg` are based on likelihood optimization routines tailored to the Vasicek regression framework. Both mean-based and quantile-based parameterizations are supported, allowing users to select the modeling strategy that best addresses their scientific objectives.

Special attention has been given to numerical robustness and interpretability of regression coefficients, ensuring reliable performance across a wide range of applications involving bounded data.

## Research impact

The `vasicekreg` package has been used in peer-reviewed methodological and applied studies focusing on parametric mean and quantile regression for bounded responses. These applications include analyses in biomedical research and COVID-19 data, demonstrating the practical relevance and versatility of the software.

Since its first release on CRAN in 2021, the package has been actively maintained and has accumulated over 12,000 downloads, indicating sustained adoption by the statistical modeling community.

## Availability

The `vasicekreg` package is openly available from the Comprehensive R Archive Network (CRAN) at\
<https://cran.r-project.org/package=vasicekreg>.

The source code is publicly hosted in a version-controlled repository and distributed under an OSI-approved open-source license.

## Acknowledgements

The authors acknowledge the support of their respective institutions and the open-source R community for fostering reproducible and transparent research software development.
