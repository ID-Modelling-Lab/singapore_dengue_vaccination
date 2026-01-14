# About
This repository contains all the scripts to generate the results presented in this modelling [study](https://www.medrxiv.org/content/10.1101/2025.10.01.25337064v1) on dengue vaccination strategy with Qdenga in Singapore.
The model is an age and serotype-stratified deterministic compartmental-based transmission model.
All the codes are written in `R (Version: R-4.3.1)`.

# Inference
Run `.../fitting/run_model_fitting_age.R` to perform model fitting to age-stratified annual dengue incidence data in a Bayesian framework using Metropolis-Hastings sampler.
The posterior distribution will be saved in `.../fitting/fitting_output`.
# Vaccine impact projection
