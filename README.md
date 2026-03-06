# Spatiotemporal Analysis of COVID-19 Incidence in Italy

Spatiotemporal analysis of COVID-19 incidence in Italy using a **Bayesian hierarchical model** combining Gaussian spatial processes and AR(1) temporal dynamics.  
The model captures regional dependencies, population effects, and epidemic phases to analyze and predict infection trends.

## Project Overview

This project investigates the regional evolution of COVID-19 incidence in Italy through a **Bayesian hierarchical spatiotemporal framework**.

The analysis combines:
- a **Gaussian spatial process** to model dependence between geographically close regions,
- an **AR(1) temporal process** to capture persistence over time,
- a **population effect** to account for regional demographic differences,
- a **phase-based decomposition** of the epidemic to compare different periods of the pandemic.

The model was implemented in **Stan** and fitted separately across epidemic phases to assess predictive performance, posterior behavior, and residual temporal structure.

## Dataset

The analysis uses official regional COVID-19 data released by the **Italian Civil Protection Department**, containing daily observations for Italian regions.

Main variables used:
- date
- region
- daily new positive cases
- latitude / longitude
- regional population

Sicily and Sardinia were excluded from the spatial analysis, and the Autonomous Provinces of Trento and Bolzano were aggregated into **Trentino-Alto Adige** when needed for map-based analyses.

## Methodology

The response variable is the **log-transformed number of daily new positive cases**.

The hierarchical model includes:
- a global intercept,
- day-of-week fixed effects,
- a population covariate,
- a latent **Gaussian spatial effect**,
- a latent **AR(1) temporal effect** for each region,
- a mixing parameter balancing spatial and temporal variability.

The workflow includes:
1. exploratory data analysis,
2. phase-specific preprocessing,
3. Bayesian inference in Stan,
4. posterior diagnostics,
5. predictive evaluation,
6. residual analysis of the temporal latent process.

## Main Results

The model captures major regional and temporal patterns of COVID-19 spread across Italy.

Key findings include:
- strong **temporal persistence** in regional contagion dynamics,
- a more limited but still relevant **spatial dependence** structure,
- heterogeneous predictive performance across epidemic phases,
- evidence that some phases require more flexible temporal structures than a simple AR(1) specification.

## Repository Structure

```text
.
├── README.md
├── covid19_spatiotemporal_model.stan
├── run_spatiotemporal_covid_model.R
├── plot_observed_vs_predicted_maps.R
├── exploratory_spatiotemporal_analysis.R
├── covid19_spatiotemporal_bayesian_analysis.pdf
└── figures/
