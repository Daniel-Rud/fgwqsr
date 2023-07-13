
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fgwqsr – Frequentist Grouped Weighted Quantile Sum Regression

<!-- badges: start -->
<!-- badges: end -->

# Overview

<br> FGWQSR estimates parameters to the following model:

$$y_i \sim \text{Bernoulli}(\pi_i) $$
$$\text{logit}(\pi_i) = c_0+ \sum\limits_{g = 1}^G \gamma_g \bigg( \sum\limits_{k = 1} ^{c_g} w_{g,k} \cdot q_{g,k,i}\bigg) + \sum_{r = 1}^R \phi_rz_{r,i} $$
where

- $y_i$ - outcome variable (coded 0 1)

- $c_0$ - offset variable (intercept)

- $G$ - total number of chemical mixture groups

- $\gamma_g$ - group index effect corresponding to mixture group $g$

- $c_g$ - total number of chemicals in group $g$

- $w_{g,k}$ - weight for chemical $k$ in group $g$

- $q_{g,k,i}$ - quantized chemical exposure $k$ in group $g$ for
  individual $i$

- $R$ - total number of adjusting covariates

- $\phi_r$ - effect of $r^{\text{th}}$ adjusting covariate

- $z_{r,i}$ - measure of adjusting covariate $r$ for individual $i$

subject to the constraints
$$\sum_{k = 1}^{c_g} w_{g,k} =1, \ w_{g,k} \in (0,1)$$

<br>

## Installation

You can install the development version of fgwqsr from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Daniel-Rud/fgwqsr")
```

``` r
library(fgwqsr)
```

# Data Generation

<br> To begin, let us generate some data to use. We will have 14
underlying chemical exposure variables, with mixture groups of sizes
(5,4,5) and weight distribution

- $\mathbf{w}_1 = (1/3,1/3,1/3,0,0)$

- $\mathbf{w}_2 = (1/2,1/2, 0, 0)$

- $\mathbf{w}_3 = (1/3,1/3,1/3,0,0)$.

We will use a fixed correlation structure such that the correlation of
chemicals within a group is .7 and the cross correlations between 2
chemicals in different groups is .5. <br> We first define ESS and CCS.
ESS is a list that contains the the sample size and the true underlying
weight distribution of the chemicals in each group. CCS is a two element
vector where the first index defines the within group correlation and
the second index controls the between group correlation. The sample size
for our dataset is fixed to $n = 10,000$ and we set the distribution of
the true underlying weights as follows:

<br> We set the group index effects for the three groups as follows (in
OR scale): $\exp(\gamma_1) = .5$, $\exp(\gamma_2) = 1$,
$\exp(\gamma_3) = 1.5$

``` r
n = 10000

gamma = log(c(.5,1,1.5)) # group index sizes in log odds scale

ess = list(n = n, 
           weights = list(w1 = c(1/3,1/3,1/3,0,0), 
                          w2 = c(1/2,1/2, 0, 0), 
                          w3 = c(1/3,1/3,1/3,0,0)
                          )
           )

ccs = c(.5,.1)
```

<br>

We create a function to create the desired correlation matrix.

``` r
# function to create correlation matrix 
create_corr_mat = function(ESS, CCS)
{
  num_pollutants = length(unlist(ESS$weights))
  num_groups = length(ESS$weights)
  
  current_index = 1
  
  corr_mat = diag(num_pollutants)
  
  # populate within group correlation blocks 
  for(i in 1:num_groups) # for each of the pollutant groups
  {
    corr_mat[current_index:(current_index + length(ESS$weights[[i]])-1), 
             current_index:(current_index + length(ESS$weights[[i]])-1) ] = CCS[1]
    current_index = current_index + length(ESS$weights[[i]])
  }
  
  # populate between group pollutant correlations 
  # set diagonal to 1 
  for(i in 1:num_pollutants)
  {
    corr_mat[i, which(corr_mat[i,] == 0)] = CCS[2]
    corr_mat[i,i] = 1
  }
  
  return(corr_mat)
}
```

<br>

Now we create a correlation matrix

``` r
corr = create_corr_mat(ESS = ess, CCS = ccs)

corr %>% data.frame # for viewing  
#>     X1  X2  X3  X4  X5  X6  X7  X8  X9 X10 X11 X12 X13 X14
#> 1  1.0 0.5 0.5 0.5 0.5 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
#> 2  0.5 1.0 0.5 0.5 0.5 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
#> 3  0.5 0.5 1.0 0.5 0.5 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
#> 4  0.5 0.5 0.5 1.0 0.5 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
#> 5  0.5 0.5 0.5 0.5 1.0 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
#> 6  0.1 0.1 0.1 0.1 0.1 1.0 0.5 0.5 0.5 0.1 0.1 0.1 0.1 0.1
#> 7  0.1 0.1 0.1 0.1 0.1 0.5 1.0 0.5 0.5 0.1 0.1 0.1 0.1 0.1
#> 8  0.1 0.1 0.1 0.1 0.1 0.5 0.5 1.0 0.5 0.1 0.1 0.1 0.1 0.1
#> 9  0.1 0.1 0.1 0.1 0.1 0.5 0.5 0.5 1.0 0.1 0.1 0.1 0.1 0.1
#> 10 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 1.0 0.5 0.5 0.5 0.5
#> 11 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.5 1.0 0.5 0.5 0.5
#> 12 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.5 0.5 1.0 0.5 0.5
#> 13 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.5 0.5 0.5 1.0 0.5
#> 14 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.5 0.5 0.5 0.5 1.0
```

<br>

We generate our exposure data based on the derived correlation
structure. We will use 5 quantiles, as that is the default for FGWQSR.

``` r
set.seed(1) # for reproducibility
chem_data = mvtnorm::rmvnorm(n, mean = rep(0, 14), sigma = corr)
```

<br>

Next, we quantize the exposure variables in order to make comparisons
between the true underlying oracle model and FGWQSR estimates. Note that
we will subtract 1 from our quantized variables, so that the first
quintile corresponds to the measure of 0. By doing this, we can ensure
that we have an interpretable model intercept.

``` r
chem_data = apply(chem_data, MARGIN = 2, statar::xtile, n = 5) - 1

head(chem_data %>% data.frame)
#>   X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 X14
#> 1  1  2  1  4  3  1  3  3  3   1   4   2   1   0
#> 2  4  3  3  4  4  4  4  4  3   0   2   1   1   0
#> 3  2  3  4  2  3  1  0  1  1   2   4   3   2   2
#> 4  3  3  1  1  3  4  3  4  3   2   3   1   4   4
#> 5  2  1  3  2  4  2  3  2  1   3   1   4   3   4
#> 6  2  0  2  0  0  2  1  2  2   1   1   1   3   0
```

<br>

Finally, we create the logistic outcome variable

``` r

intercept = 0

# create logit(pi) outcome 
logit_pi =  gamma[1] * (chem_data[, 1:5] %*% ess$weights$w1) + 
  gamma[2] * (chem_data[, 6:9] %*% ess$weights$w2) + 
  gamma[3] * (chem_data[, 10:14] %*% ess$weights$w3) + intercept

# transform to pi 
pi = stats::plogis(logit_pi)

# transform to bernoulli outcome 
y = sapply(pi, FUN = function(p) rbinom(1,1,p))

# create dataset 
data = data.frame(y = y, chem_data)

# view dataset
head(data)
#>   y X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 X14
#> 1 0  1  2  1  4  3  1  3  3  3   1   4   2   1   0
#> 2 0  4  3  3  4  4  4  4  4  3   0   2   1   1   0
#> 3 0  2  3  4  2  3  1  0  1  1   2   4   3   2   2
#> 4 0  3  3  1  1  3  4  3  4  3   2   3   1   4   4
#> 5 0  2  1  3  2  4  2  3  2  1   3   1   4   3   4
#> 6 1  2  0  2  0  0  2  1  2  2   1   1   1   3   0
```

<br> Now, we are ready to fit a FGWQSR model

<br>

# Fitting using FGWQSR

<br> From the FGWQSR package, the main function we will be using to fit
models if the fgwqsr() function. The model formula that we specify will
be different than traditional formulas in lm and glm, as we will need to
denote our mixture groups. Three special characters are used in fgwqsr
formulas: <br>

- `|` denotes the boundary of a mixture group, used to seperate
  chemicals within a mixture group.

- `/` - denotes the end of the mixture group specification, adjusting
  covariates can be added to the formula after this character. If no
  adjusting covariates, do not need to specify.

- `i.` - precedes categorical variables to denote a categorical
  variable. For example, if we have a categorical variable cat_var, we
  would denote this in the model formula by i.cat_var. This is similar
  to the stata syntax to declare categorical variables.

<br> The fgwqsr() function has other options too: <br>

- `data` - should be a data frame object containing variable columnames
  referenced in model formula (cat vars do not need to be named with i.
  in columnames)

- `quantiles` - number of quantiles to quantize the exposure variables
  in the mixture portion of the model

- `n_mvn_sims` - defines resolution for simulated null distribution for
  group index and single chemical LRTs. Default is 10,000.

- `zero_threshold_cutoff` - defines a .

- `verbose` - Displays messages and progress bar while fitting FGWQSR
  model. Default is TRUE.

- `cores` - number of cores to parallelize on for fitting nested models
  and simulated null LRT distributions. Default is number of available
  cores on user device.

- `optim_control_list` - option to supply control options to optim.

<!-- <br> -->
<!-- First, we will focus on fitting models with only chemical mixture variables.  To do this, we will supply the model formula as follows:  -->
<!-- ```{r} -->
<!-- mod_formula = y ~ X1 + X2 + X3 + X4 + X5 | X6 + X7 + X8 + X9 | X10 + X11 + X12 + X13 + X14 -->
<!-- ``` -->
<!-- <br> -->
<!-- Notice that we did not use the "$/$" character since we did not include adjusting covariates.  Now, we can fit the model using the function fgwqsr().   -->
<!-- ```{r} -->
<!-- fgwqsr_fit = fgwqsr(formula = mod_formula,  -->
<!--                     data = data,  -->
<!--                     quantiles = 5,  -->
<!--                     n_mvn_sims = 100,  -->
<!--                     verbose = T) -->
<!-- ``` -->
<!-- <br> -->
<!-- We can see the model summary using the call summary() -->
<!-- ```{r} -->
<!-- summary(fgwqsr_fit) -->
<!-- ``` -->
<!-- <br> -->
<!-- We can compare the true underlying group indices and chemical weights with their estimates, along with if the correct inference is made.  Note that dashes are inserted for the true weights of group 2, since the true underlying group index is null (thus, the weights are unidentifiable).   -->
<!-- ```{r} -->
<!-- group_index_frame = data.frame("True Group Index" = exp(gamma),  -->
<!--                           "Estimated Group Index" = exp(fgwqsr_fit$inference_frames$group_index_frame[,1]) %>% round(3),  -->
<!--                           "Signficiant?" = ifelse(fgwqsr_fit$inference_frames$group_index_frame[,3] < .05, "Yes", "No"), check.names = F) -->
<!-- rownames(group_index_frame) = rownames(fgwqsr_fit$inference_frames$group_index_frame) -->
<!-- true_weights = ess$weights %>% unlist %>% round(3) -->
<!-- true_weights[6:9] = "-" -->
<!-- weight_frame = data.frame("True Weight" = true_weights,  -->
<!--                           "Estimated Weight" = fgwqsr_fit$inference_frames$weight_frame[,1] %>% round(3),  -->
<!--                           "Signficiant?" = ifelse(fgwqsr_fit$inference_frames$weight_frame[,3] < .05, "Yes", "No"), check.names = F) -->
<!-- group_index_frame; weight_frame -->
<!-- ``` -->
<!-- <br> -->
<!-- # Fitting with Adjusting Continuous and Categorical Covariates  -->
<!-- As mentioned earlier, the "$/$" character in the model formula sent to fgwqsr() indicates that adjusting covariates follow this character.  Continuous characters can be referenced with their columnname from the dataset, while categorical variables need to be referenced using the prefix "i.".  We will illustrate this through an example where we create a continuous $weight$ variable and a categorical $city$ variable.  We will give an effect (in log odds scale) of .5 to weight, .2 to the comparison between city 2 and reference city 1, and -1 to the comparison between city 3 and city 1.    -->
<!-- <br> -->
<!-- ```{r} -->
<!-- set.seed(1) -->
<!-- # create adjusting covariates  -->
<!-- weight = rnorm(n = n, mean = 68, sd = 2.5)  -->
<!-- city = sample(c("city_1", "city_2", "city_3"), size = n, replace = T) -->
<!-- # need to adjust intercept for case control ratio  -->
<!-- intercept = -3.3 -->
<!-- # create logit(pi) outcome WITH adjusting covariates  -->
<!-- logit_pi =  gamma[1] * (chem_data[, 1:5] %*% ess$weights$w1) +  -->
<!--   gamma[2] * (chem_data[, 6:9] %*% ess$weights$w2) +  -->
<!--   gamma[3] * (chem_data[, 10:14] %*% ess$weights$w3) +  -->
<!--   .05 *weight +  -->
<!--   .2*ifelse(city == "city_2", 1, 0) + -1*ifelse(city == "city_3", 1, 0) +  -->
<!--   intercept -->
<!-- # transform to pi  -->
<!-- pi = stats::plogis(logit_pi) -->
<!-- # transform to bernoulli outcome  -->
<!-- y = sapply(pi, FUN = function(p) rbinom(1,1,p)) -->
<!-- # create dataset  -->
<!-- data = data.frame(y = y, chem_data, weight = weight, city = city) -->
<!-- head(data) -->
<!-- ``` -->
<!-- <br> -->
<!-- Now, we specify the formula to include the continuous and categorical variables and call fgwqsr().  -->
<!-- ```{r} -->
<!-- mod_formula_adj = y ~ X1 + X2 + X3 + X4 + X5 | X6 + X7 + X8 + X9 | X10 + X11 + X12 + X13 + X14 / weight + i.city -->
<!-- fgwqsr_fit_adj = fgwqsr(formula = mod_formula_adj,  -->
<!--                     data = data,  -->
<!--                     quantiles = 5,  -->
<!--                     n_mvn_sims = 100,  -->
<!--                     verbose = T) -->
<!-- ``` -->
<!-- <br> -->
<!-- Again, we use the function summary() to view the results of the fgwqsr() call.   -->
<!-- ```{r} -->
<!-- summary(fgwqsr_fit_adj) -->
<!-- ``` -->
<!-- <br> -->
<!-- Finally, we can compare the true parameter value with our estimates from FGWQSR.   -->
<!-- ```{r} -->
<!-- group_index_frame = data.frame("True Group Index" = exp(gamma),  -->
<!--                           "Estimated Group Index" = exp(fgwqsr_fit_adj$inference_frames$group_index_frame$Estimate) %>% round(3),  -->
<!--                           "Signficiant?" = ifelse(fgwqsr_fit_adj$inference_frames$group_index_frame$`P-value` < .05, "Yes", "No"), check.names = F) -->
<!-- rownames(group_index_frame) = rownames(fgwqsr_fit_adj$inference_frames$group_index_frame) -->
<!-- true_weights = ess$weights %>% unlist %>% round(3) -->
<!-- true_weights[6:9] = "-" -->
<!-- weight_frame = data.frame("True Weight" = true_weights,  -->
<!--                           "Estimated Weight" = fgwqsr_fit_adj$inference_frames$weight_frame$`Weight Estimate` %>% round(3),  -->
<!--                           "Signficiant?" = ifelse(fgwqsr_fit_adj$inference_frames$weight_frame$`P-value` < .05, "Yes", "No"), check.names = F) -->
<!-- adj_cov_frame = data.frame("True Covariate Effect" = c(-3.3, .5,.2,-1),  -->
<!--                            "Estimated Covariate Effect" = fgwqsr_fit_adj$inference_frames$adj_param_frame$Estimate,  -->
<!--                            "Significant?" = ifelse(fgwqsr_fit_adj$inference_frames$adj_param_frame$`P(Z > |z|)` < .05, "Yes", "No"), check.names = F) -->
<!-- rownames(adj_cov_frame) = rownames(fgwqsr_fit_adj$inference_frames$adj_param_frame) -->
<!-- group_index_frame; weight_frame; adj_cov_frame -->
<!-- ``` -->
<!-- # Fitting Models with BGWQSR  -->
<!-- We also provide in this package functions to run BGWQSR models using the runjags package.  The bgwqsr() function takes a model formulation similar to that in fgwqsr(), also utilizing the special characters "$|$", "$/$", and "i.".  Fitting the BGWQSR model using runjags allows us to leverage the "parallel" method in runjags that allows us to fit independent mcmc chains on multiple cores, speeding up fitting time.   -->
<!-- ```{r} -->
<!-- bgwqsr_fit = bgwqsr(formula = mod_formula_adj,  -->
<!--                     data = data,  -->
<!--                     quantiles = 5,  -->
<!--                     n.iter = 1000, -->
<!--                     n.burnin = 40,  -->
<!--                     n.thin = 1, n.chains = 3,  -->
<!--                     n.adapt = 10,  -->
<!--                     method = "parallel") -->
<!-- ``` -->
<!-- <br> -->
<!-- To retrieve the summary measures from the MCMC sampling, we can extract the summaries object from the bgwqsr model.    -->
<!-- ```{r} -->
<!-- bgwqsr_fit$model$summaries -->
<!-- ``` -->
<!-- We can analyze the mixing of the markov chains and corresponding autocorrelation plots using functions from the coda package.  We will  -->
<!-- output only the first few plots. -->
<!-- ```{r, figure.width = 10, figure.height = 10} -->
<!-- coda::traceplot(bgwqsr_fit$model$mcmc) # traceplot -->
<!-- ``` -->
<!-- ```{r, figure.width = 10, figure.height = 10} -->
<!-- coda::autocorr.plot(bgwqsr_fit$model$mcmc, auto.layout = F, ask = F) # autocorrelation plot  -->
<!-- ``` -->
<!-- ```{r, figure.width = 10, figure.height = 10} -->
<!-- coda::densplot(bgwqsr_fit$model$mcmc) # posterior density plots -->
<!-- ``` -->
<!-- ```{r, figure.width = 10, figure.height = 10} -->
<!-- plot(bgwqsr_fit$model$mcmc, auto.layout = F) # combines traceplot() and densplot() -->
<!-- ``` -->
<!-- Finally, we can look at plots for the posterior credible intervals for group indices, single chemical weights, and corresponding posterior means.  To do this, we can use the functions plotWeights(), plotBetas(), and plotResult().  plotResults() combines both plots generates by plotWeights and plotBetas into a side by side figure.  Below are code examples for each of the three function calls.   -->
<!-- ```{r} -->
<!-- plot_betas(bgwqsr_fit) -->
<!-- ``` -->
<!-- ```{r} -->
<!-- plot_weights(bgwqsr_fit) -->
<!-- ``` -->
<!-- ```{r, fig.width = 12} -->
<!-- plot_results(bgwqsr_fit) -->
<!-- ``` -->
