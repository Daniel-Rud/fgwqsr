---
title: 'fgwqsr: An R package for Frequentist Grouped Weighted Quantile Sum Regression '
tags:
  - R
  - Weighted Quantile Sum Regression 
  - chemical/pollutant mixture modeling 
  - correlated data
  - nonregular likelihood asymptotics 
authors:
  - name: Daniel Rud 
    orcid: 0000-0002-0508-4552
    equal-contrib: true
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Juan Pablo Lewinger 
    orcid: 0000-0002-0692-2570
    equal-contrib: true
    affiliation: "1" # (Multiple affiliations must be quoted)
    
affiliations:
 - name: Department of Population and Public Health Sciences, University of Southern California, USA
   index: 1
date: 3 September 2024
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary
**fix FG citation once finished in bib file **

Researchers in the area of environmental epidemiology are often interested in analyzing the effects of chemical/pollutant exposures on diverse health related outcomes.  Aside from fitting single constituent models, many of the popularized epidemiological methods for handling chemical/pollutant *mixtures* are designed to accommodate the joint adjustment of multiple exposures [@wqsr; @gwqsr; @bgwqsr; @Environmental_exposure_mixtures; @an_overview_of_methods; @sandy_cit].  Models for handling mixtures are often designed with consideration to the exposure data, which often possess high levels of correlation as a result of mixture constituents deriving from common sources.  Out of the existing methods, the Weighted Quantile Sum Regression (WQSR) has been heavily popularized and utilized to examine associations between exposure mixtures and health outcomes [@wqsr; @gwqsr; @bgwqsr].  WQSR estimates group effects that measure the magnitude of effect of a mixture group and sets of group weights that measure the relative contribution of each constituent from a mixture group. In the binary outcome scenario, the WQSR model can be mathematically formulated as 

$$
\begin{aligned}
y_i &\sim \text{Bernoulli}(\pi_i) \\ 
\text{logit}(\pi_i) &= c_0+ \sum\limits_{g = 1}^G \gamma_g \bigg( \sum\limits_{k = 1} ^{c_g} w_{g,k} \cdot q_{g,k,i}\bigg) + \sum_{r = 1}^R \phi_rz_{r,i} 
\end{aligned}
$$
where, for subject  $i$, $y_i$  represents the observed disease outcome,  $\pi_i$ the probability of disease,   $q_{g,k,i}$ the exposure to chemical $k$ in mixture group $g$, and $z_{r,i}$ represents the the $r^{th}$ confounder adjustment for individual $i$.  The weights for mixture group $g$ satisfy $\sum_{k=1}^{c_g} w_{g,k} = 1$ and $w_{g,k} \ge 0$.  WQSR models are constrained such that all constituents from a particular mixture group have effects in the same direction, which functions as a form of regularization to stabilize the effect estimates of the highly correlated exposures.  


In the `fgwqsr` package, we contribute software to fit the Frequentist Grouped Weighted Quantile Sum Regression (FGWQSR) model as described in @fgwqsr.  The main function of the package is the function `fgwqsr`.  FGWQSR can accommodate binary, continuous, and count outcome types.  To fit a FGWQSR model, aside from having created a dataframe containing all the relevant data, one needs to create a special model formula.  The model formula must contain vertical bars `|` in order to denote the separation between mixture group elements and a forward slash `/` to denote the separation between the mixture groups and unconstrained adjusting covariates.  In addition, the prefix `i.` is used in the model formula after the forward slash to indicate categorical adjusting covariates.  For example, consider an analysis on the outcome $Y$ with the following two mixture groups $\{A_1, A_2\}$ and $\{B_1, B_2\}$ and confounders $\{W_1, W_2, W_3\}$ where $W_1, W_2$ are numeric and $W_3$ is a three level categorical variable.  Then, with the data in a dataframe named `data`, the model formula will be 
```r
model_formula = Y ~ A1 + A2 | B1 + B2 / W_1 + W_2 + i.W_3
```

If no adjusting covariates are desired in the analysis, no forward slash is required.  In this case, the model formula would be 

```r
model_formula = Y ~ A1 + A2 | B1 + B2 
```

In the case of only fitting a FGWQSR model with a single mixture group, one does not need to use vertical bars.  For example, if we exclude $\{ B_1, B_2\}$, then the model formula would be 

```r
model_formula = Y ~ A1 + A2 / W_1 + W_2 + i.W_3
```

Given `model_formula`, the outcome family type `family` being one of `("binomial", "gaussian", "poisson")`, the number of quantiles `q` desired for the quantization of the mixture constituents, the number of multivariate normal simulations `n_mvn_rep` performed for each hypothesis test, and the number of cores `cores` one is willing to parallelize over, an FGWQSR model can be fitted with the call 

```r
fgwqsr_model = fgwqsr(formula = model_formula, 
                      data = data, 
                      quantiles = q, 
                      family = family, 
                      n_mvn_sims = n_mvn_rep, 
                      verbose = T, 
                      cores = cores)
```
The resulting FGWQSR model can be viewed using the summary function `summary(fgwqsr_model)`, where parameter estimates for group effects and sets of group weights are presented, along with statistical tests for group effects and single constituent effects. 

An optional tuning parameter is the `zero_threshold_cutoff`, which is used in the nonregular statistical testing procedure.  The `zero_threshold_cutoff` parameter defines how often parameters estimated close to the boundary of the parameter space are assigned a boundary cone in the constrained multivariate normal monte carlo inference procedure.  Reasonable values of the parameter may be within [0.05, 0.5]; however the default value of $0.5$ has been shown to work well in a variety of situations, so user tuning is not necessary.  For more information on this parameter, please see @fgwqsr. 

In addition to an implementation of the FGWQSR model, we also include in this package an implementation of the Bayesian Grouped Weighted Quantile Sum Regression (BGWQSR) model for binary outcomes. In contrast to the implementation of BGWQSR in the `BayesGWQS` package, our implementation leverages the `runjags` Markov Chain Monte Carlo (MCMC) package to run separate MCMC chain processes in parallel.  BGWQSR models can be fit in the `fgwqsr` package using the `bgwqsr` function, where several parameters to control the MCMC settings can be optionally supplied.  In addition, the `plot_results`, `plot_betas`, and `plot_weights` functions can be used to view the resulting group effect, weight, and confounder estimates and posterior credible intervals. 

 For more information on fitting `fgwqsr`  and `bgwqsr` models, we refer the reader to the vignette [here](https://github.com/Daniel-Rud/fgwqsr).


# Statement of need

FGWQSR was motivated by several limitations of existing WQSR approaches.  First off, many of the current methods [@wqsr; @gwqsr] require data splitting; that is, they require independent data sets to first estimate sets of group weights and then to estimate group effect parameters. FGWQSR outshines many of the existing methods as group effects and group weights are jointly fitted using a constrained optimization procedure detailed in @fgwqsr that does not require any data splitting.  Secondly, many of the existing implementations of WQSR do not easily accommodate large datasets.  FGWQSR was created with the intention to apply WQSR to a dataset containing $317,000$ observations, which existing implementations of WQSR could not handle.  Finally, FGWQSR further develops the WQSR statistical model by producing statistical tests for both group effects and single constituent effects, where WQSR predecessors only focused on statistical inference of group effects.  Thus, FGWQSR is a state of the art method for performing WQSR type analyses that does not require data splitting, can handle large datasets, and conducts statistical tests for group and single constituent effects.  




# Installation 

The most current version of `fgwqsr` package can be downloaded from github using the following instructions: 
```r
install.packages("devtools")
devtools::install_github("Daniel-Rud/fgwqsr")
```

# References







