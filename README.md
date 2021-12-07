# wem
Wastewater Epidemic Model



[![Lifecycle:
development](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental-1)


## Introduction 

This repo provides the code for the R package `wem`, which stands for "Wastewater Epidemic Model". 

It implements an epidemic model that can include pathogen concentration in wastewater as a source of data in addition to the more traditional case incidence and hospitalization data. This model was designed with SARS-CoV-2/COVID-19 in mind, hence default parameters and model structure reflect this. The model and its implementation in `wem` is developed at the Public Health Agency of Canada / National Microbiology Laboratory. Please note this software is provivded "as is", without warranty of any kind; see the [license](LICENSE).

The epidemic model that represents pathogen transmission in the population is a deterministic [SEIR-type compartmental model](https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology). 

A simple advection-dispersion-decay model is used to simulate the fate of SARS-CoV-2
along its journey in wastewater from the shedding points to the sampling site. 
This model is a combination of an exponential viral decay and a dispersed plug-flow function representing all possible hydrodynamic processes (e.g., dilution, sedimentation and
resuspension) that leads to RNA degradation as well as decrease and delay of signal at the time of sampling.

For more details on the mathematical model see the file `model.pdf` in the package Help, "other documentation" section. See also [Nourbakhsh et al. (medRxiv 2021)](https://www.medrxiv.org/content/10.1101/2021.07.19.21260773v1). Note that this R package implements additional features compared to the ones presented in that publication (for example, vaccination is a new addition).

For practical examples, please check the vignettes.

## Installation

The R package `wem` is not on CRAN (yet), hence to install it simply type in a R console:

``` r
devtools::install_github(repo = 'phac-nml-phrsd/wem', build_vignettes = FALSE)
```

(we recommend to avoid building vignettes, as they take several minutes to compile. The HTML files for the vignettes are available in the `vignettes` folder in this repo )
