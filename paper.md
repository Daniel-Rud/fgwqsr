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
date: 24 March 2025
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Environmental epidemiologists frequently study the effects of chemical and pollutant
exposures on health outcomes. Beyond single-constituent models, recent
epidemiological methods focus on modeling exposure mixtures jointly, accounting for the correlation
between exposures arising from common sources [@wqsr; @gwqsr; @bgwqsr; @Environmental_exposure_mixtures; @an_overview_of_methods; @sandy_cit].


Among existing approaches, Weighted Quantile Sum Regression (WQSR) has gained traction for evaluating associations between exposure mixtures and health outcomes [@wqsr; @gwqsr; @bgwqsr]. WQSR estimates both (1) group effects, which quantify the impact of a mixture group, and (2) sets of group weights, which represent the relative contributions of individual constituents within a mixture. In the binary outcome setting, the WQSR model is formulated as:

$$
\begin{aligned}
y_i &\sim \text{Bernoulli}(\pi_i) \\ 
\text{logit}(\pi_i) &= c_0+ \sum\limits_{g = 1}^G \gamma_g \bigg( \sum\limits_{k = 1} ^{c_g} w_{g,k} \cdot q_{g,k,i}\bigg) + \sum_{r = 1}^R \phi_rz_{r,i} 
\end{aligned}
$$
where, for subject  $i$, $y_i$  represents the observed disease outcome,  $\pi_i$ the probability of disease,   $q_{g,k,i}$ the exposure to chemical $k$ in mixture group $g$, and $z_{r,i}$ the $r^{th}$ confounder adjustment.  The weights for mixture group $g$ satisfy $\sum_{k=1}^{c_g} w_{g,k} = 1$ and $w_{g,k} \ge 0$. The parameter $\gamma_g$ represents the group effect for a given mixture group, capturing the impact of a one-quantile increase in all chemical constituents within the group.  WQSR models are constrained such that all constituents from a particular mixture group have effects in the same direction, which functions as a form of regularization to stabilize the effect estimates of the highly correlated exposures.  


The `fgwqsr` package implements the Frequentist Grouped Weighted Quantile Sum Regression (FGWQSR) model introduced in @fgwqsr. Its main function, `fgwqsr`, accommodates binary, continuous, and count outcome types. To fit a FGWQSR model, users must specify a special model formula using vertical bars (`|`) to separate mixture group elements and a forward slash (`/`) to separate mixture groups from unconstrained covariates. Categorical covariates must be prefixed with `i.`. For instance, in an analysis with outcome $Y$, mixture groups $\{A_1, A_2\}$ and $\{B_1, B_2\}$, and confounders $\{W_1, W_2, W_3\}$ (where $W_1, W_2$ are numeric and $W_3$ is categorical), the model formula is:
```r
model_formula = Y ~ A1 + A2 | B1 + B2 / W_1 + W_2 + i.W_3
```

If no adjusting covariates are included, no forward slash is required:

```r
model_formula = Y ~ A1 + A2 | B1 + B2 
```

For a single mixture group, vertical bars are not necessary:

```r
model_formula = Y ~ A1 + A2 / W_1 + W_2 + i.W_3
```

Given `model_formula`, the outcome family type `family` being one of (`"binomial"`, `"gaussian"`, `"poisson"`), the number of quantiles `q` desired for the quantization of the mixture constituents, the number of multivariate normal simulations `n_mvn_rep` performed for each hypothesis test, and the number of cores `cores` one is willing to parallelize over, an FGWQSR model can be fitted with the call:

```r
fgwqsr_model = fgwqsr(formula = model_formula, 
                      data = data, 
                      quantiles = q, 
                      family = family, 
                      n_mvn_sims = n_mvn_rep, 
                      verbose = T, 
                      cores = cores)
```
Results can be examined using `summary(fgwqsr_model)`, which provides parameter estimates for group effects, group weights, and statistical tests for both group and single-constituent effects.

An optional tuning parameter, `zero_threshold_cutoff`, is used in the non-regular statistical testing procedure. This parameter determines how often near-boundary estimates are assigned a boundary cone in the constrained multivariate normal Monte Carlo inference procedure. A default value of 0.5 has been shown to perform well across various scenarios, though reasonable values range from [0.05, 0.5]. More details are provided in @fgwqsr.

In addition to FGWQSR, the package includes an implementation of Bayesian Grouped Weighted Quantile Sum Regression (BGWQSR) for binary outcomes. Unlike the BayesGWQS package, our implementation leverages the `runjags` package for parallelized Markov Chain Monte Carlo (MCMC) sampling. BGWQSR models can be fitted using the `bgwqsr` function, with additional MCMC control parameters available. Visualization tools such as `plot_results`, `plot_betas`, and `plot_weights` provide graphical summaries of group effects, weights, and confounder estimates with posterior credible intervals.

For further guidance, see the package vignette [here](https://github.com/Daniel-Rud/fgwqsr).


# Statement of need

FGWQSR was developed to address several limitations of existing WQSR methods. Many existing approaches [@wqsr; @gwqsr] rely on data splitting, requiring separate datasets to first estimate group weights and then assess group effect parameters. In contrast, FGWQSR jointly estimates group effects and group weights using a constrained optimization procedure [@fgwqsr], eliminating the need for data splitting.

Moreover, existing WQSR implementations struggle with large datasets. FGWQSR was designed to efficiently handle large datasets and was successfully applied to a dataset with 317,767 observations, which previous implementations struggled to accommodate. Additionally, FGWQSR extends the statistical framework by introducing statistical tests for both group and single-constituent effects, whereas previous WQSR models focused solely on group effects.

Thus, FGWQSR represents a significant advancement in WQSR methodology, providing a scalable, statistically rigorous approach that does not require data splitting, handles large datasets, and enables statistical inference for both group and individual constituent effects.




# Installation 

The most current version of `fgwqsr` package can be downloaded from github using the following instructions: 
```r
install.packages("devtools")
devtools::install_github("Daniel-Rud/fgwqsr")
```

# References







