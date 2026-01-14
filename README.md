# Adult dengue vaccination in a low transmission setting: A modelling study in Singapore
This repository contains all the scripts to generate the results presented in this modelling [study](https://www.medrxiv.org/content/10.1101/2025.10.01.25337064v1) on dengue vaccination strategy with Qdenga in Singapore.
The model is an age and serotype-stratified deterministic compartmental-based transmission model.
All the codes are written in `R (Version: R-4.3.1)`. The libraries required to run the codes are including : [odin](https://mrc-ide.github.io/odin/), [odin.dust](https://mrc-ide.github.io/odin.dust/index.html), [dust](https://mrc-ide.github.io/dust/index.html), [mcstate](https://mrc-ide.github.io/mcstate/index.html). 

# Inference
Run `.../fitting/run_model_fitting_age.R` to perform model fitting to age-stratified annual dengue incidence data in a Bayesian framework using Metropolis-Hastings sampler.
The posterior distribution will be saved in `.../fitting/fitting_output`.
# Vaccine impact estimation
Run `.../projection/run_vac_projection.R` to generate the impact estimates presented in the study.
The output will be saved as a `.rds` file in `.../model_output`.
# Plotting
To generate the figures run `.../projection/plot_for_paper.rmd`.
Each chunk will generate a figure presented in the study. The descriptions with figure labels are given in each chunk.
