msocc is an R package for fitting and analyzing computationally efficient Bayesian multi-scale occupancy models. Its development was motivated by the use of environmental DNA (eDNA) for the monitoring of pathogens, parasites, and invasive species in the Yellowstone River. Consequently, the language in this package assumes this context. For more information on eDNA, see [Improved detection of rare, endangered and invasive trout in using a new large-volume sampling method for eDNA capture](https://onlinelibrary.wiley.com/doi/epdf/10.1002/edn3.23) or [Adding invasive species biosurveillance to the U.S. GeologicalSurvey streamgage network](https://esajournals.onlinelibrary.wiley.com/doi/epdf/10.1002/ecs2.2843). 

# Installation instructions
This package is still under development, but can be installed through GitHub using the following code:

```{r, eval = F}
install.packages('devtools') # only needed if devtools is not currently installed
devtools::install_github('StrattonCh/msocc')
```

# Contact information
Christian Stratton (christianstratton@montana.edu) developed this R package.

# Usage

The heavy lifting for this package is done with `msocc_mod`, which fits the models. The results of this function are then passed to `posterior_summary` to numerically summarize the posterior distribution and `cred_plot` to visually summarize it. Additionally, this package provides tools to simulate data from multi-scale occupancy models.

## Simulated examples
To showcase the utility of this package, we begin with a simple example where there are 10 sites of interest, from which we collect 5 samples each and analyze 5 PCR replicates for the presence of the target DNA. To simulated data consistent with this structure, we use the following code. 

```{r, include = F}
library(msocc)
```

```{r}
sim <- msocc_sim(M = 10, J = 5, K = 5)
```
