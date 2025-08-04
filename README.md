# Model-based comparison of universal and risk-based COVID-19 vaccine schedules in the United States, 2023-2025

This repository contains analytic code for comparing the risk of severe COVID-19 under different vaccination scenarios (including eligibility based on age/risk group, frequency, uptake), accounting for differential waning of protection against non-severe and severe disease and heterogeneous risk by age group and immunocompromised status.

Data sources used for this analysis is publicly available. All data used for this analysis can be found in the `data` folder.

This study is in-progress.

## Structure
* `data`: contains all data used in this analysis (both raw and processed forms)
* `analysis`:
  * `1-waning model`: contains code for constructing the waning protection curves by risk group, and waning curves used for sensitivity analyses
  * `2-model initialization`: contains code for constructing a hypothetical cohort of 50 million individuals (for main analysis) and 1 million individuals (for model calibration) that broadly resembles the demographic characteristics seen in the United States. Immune history of each individual (prior vaccination, prior infection) in the cohort is empirically estimated using historical data.
  * `3-model calibration`: contains code for calibrating the model to observed COVID-19 severe case data between July 1, 2023 - January 1, 2025. The calibration uses a MCMC approach to calibrate the model parameters (month-specific lambdas, baseline case-hospitalization fraction, time-since shift).
  * `4-main analysis`: contains code for running different vaccination scenarios. See the `README.md` inside this folder for more information (file not made yet).
  * `5-results processing`: contains all code needed for processing of simulation results.
  * `6-data cleaning`: contains all code needed for initial cleaning of raw data


## Software
Analysis was conducted in R (version 4.2.1). The installation time for R and library packages necessary for this study is only a few minutes. The model calibration and main analyses were run on a high-performance computing cluster (Sherlock). 

## Contact 
Please direct any questions to the study authors:

Hailey Park, Stanford University, contact: Hailey.Park@stanford.edu

Nathan Lo, Stanford University, contact: Nathan.Lo@stanford.edu