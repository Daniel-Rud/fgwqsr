
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fgwqsr – Frequentist Grouped Weighted Quantile Sum Regression

Fully frequentist model fitting procedure for the weighted quantile sum
regression formulation. Included in this package is an implementation of
Bayesian Grouped Weighted Quantile Sum Regression (bgwqsr) introduced by
Wheeler, David C et al. that fits Markov Chain Monte Carlo (MCMC) chains
in parallel through the runjags package.

<!-- badges: start -->
<!-- badges: end -->

# Overview

<br> FGWQSR estimates parameters to the following model:

<br>

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
$$\sum_{k = 1}^{c_g} w_{g,k} =1, \ w_{g,k} \in (0,1)$$ without requiring
.

<br>

# Installation

You can install the development version of fgwqsr from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Daniel-Rud/fgwqsr")
```

``` r
library(fgwqsr)
```

In order to use fgwqsr, let us first generate some data

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
for our dataset is fixed to $n = 100,000$ and we set the distribution of
the true underlying weights as follows:

<br> We set the group index effects for the three groups as follows (in
OR scale): $\exp(\gamma_1) = .5$, $\exp(\gamma_2) = 1$,
$\exp(\gamma_3) = 1.5$

``` r
n = 100000

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
set.seed(2023) # for reproducibility
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
#> 1  1  0  0  0  0  4  1  3  2   1   3   2   3   3
#> 2  2  3  3  3  4  3  2  2  4   1   0   0   0   0
#> 3  4  4  3  4  3  2  1  2  3   4   4   3   4   4
#> 4  4  3  4  4  4  4  4  4  4   2   3   4   4   4
#> 5  0  1  2  3  3  1  0  0  3   3   4   2   4   4
#> 6  0  0  1  0  2  2  1  1  2   1   3   4   2   3
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
#> 1 1  1  0  0  0  0  4  1  3  2   1   3   2   3   3
#> 2 1  2  3  3  3  4  3  2  2  4   1   0   0   0   0
#> 3 0  4  4  3  4  3  2  1  2  3   4   4   3   4   4
#> 4 0  4  3  4  4  4  4  4  4  4   2   3   4   4   4
#> 5 1  0  1  2  3  3  1  0  0  3   3   4   2   4   4
#> 6 0  0  0  1  0  2  2  1  1  2   1   3   4   2   3
```

<br> Now, we are ready to fit a FGWQSR model

<br>

# Fitting using FGWQSR

<br> From the FGWQSR package, the main function we will be using to fit
models if the `fgwqsr()` function. The model formula that we specify
will be different than traditional formulas in lm and glm, as we will
need to denote our mixture groups. Three special characters are used in
fgwqsr formulas: <br>

- `|` - denotes the boundary of a mixture group, used to seperate
  chemicals within a mixture group.

- `/` - denotes the end of the mixture group specification, adjusting
  covariates can be added to the formula after this character. If no
  adjusting covariates, do not need to specify.

- `i.` - precedes categorical variables to denote a categorical
  variable. For example, if we have a categorical variable cat_var, we
  would denote this in the model formula by i.cat_var. This is similar
  to the stata syntax to declare categorical variables.

<br> The `fgwqsr()` function has other options too: <br>

- `data` - a data frame object containing variable columnames referenced
  in model formula (cat vars do not need to be named with i. in column
  names).

- `quantiles` - number of quantiles to quantize the exposure variables
  in the mixture portion of the model.

- `n_mvn_sims` - defines resolution for simulated null distribution for
  group index and single chemical LRTs. Default is 10,000.

- `zero_threshold_cutoff` - Value within (0,.5\] that defines how often
  parameters estimated close to the boundary of the parameter space are
  assigned a boundary cone in the constrained multivariate normal monte
  carlo inference procedure. A value of 0.5 will assign parameters with
  FGWQSR maximum likelihood estimates of precisely 0 the boundary cone
  while a value of 0 will assign all parameters a boundary cone.
  Reasonable values may be within \[0.05, 0.5\] – all choices of are
  asymptotically equivalent. The default is set to
  zero_tolerance_threshold = 0.5.

- `verbose` - Displays messages and progress bar while fitting FGWQSR
  model. Default is TRUE.

- `cores` - number of cores to parallelize on for fitting nested models
  and simulated null LRT distributions. Default is number of available
  cores on user device.

- `optim_control_list` - option to supply control options to optim.

<br>

First, we will focus on fitting models with only chemical mixture
variables. To do this, we will supply the model formula as follows:

``` r

mod_formula = y ~ X1 + X2 + X3 + X4 + X5 | X6 + X7 + X8 + X9 | X10 + X11 + X12 + X13 + X14
```

<br>

Notice that we did not use the `/` character since we did not include
adjusting covariates. Now, we can fit the model using the function
`fgwqsr()`.

``` r
fgwqsr_fit = fgwqsr(formula = mod_formula,
                    data = data,
                    quantiles = 5,
                    n_mvn_sims = 10000,
                    verbose = T)
#> 
#> Fitting full model and nested models...
#> 
#> Generating LRT distributions under H0...
```

<br>

We can see the model summary using the call `summary()`

``` r

summary(fgwqsr_fit)
#> 
#> Call: 
#> FGWQSR with formula 'y ~ X1 + X2 + X3 + X4 + X5 | X6 + X7 + X8 + X9 | X10 + X11 + X12 + X13 + X14' on n = 100000 observations.
#> 
#> 10000 samples used for simulated LRT distirbution.
#> 
#> Log Likelihood: -59431.33 | AIC: 118892.7 | BIC: 119035.3
#> 
#> Estimates and Inference for Group Index Effects
#> 
#>                    Estimate         LRT P-value    
#> Mixture Effect 1 -0.6917669 11734.25081  <2e-16 ***
#> Mixture Effect 2 -0.0053746     1.15193  0.7152    
#> Mixture Effect 3  0.3896534  3721.18808  <2e-16 ***
#> 
#> Estimates and Inference for Weights
#> 
#>    Weight Estimate        LRT P-value    
#> X1      0.33902356 1509.16542  <2e-16 ***
#> X2      0.31486078 1303.76353  <2e-16 ***
#> X3      0.34439133 1555.44392  <2e-16 ***
#> X4      0.00172433    0.04012  0.3824    
#> X5      0.00000000    0.95207  0.1651    
#> -------------------------------------------------
#>    Weight Estimate     LRT P-value
#> X6               0 0.00000  1.0000
#> X7               0 0.95207  0.1150
#> X8               0 0.00000  0.4596
#> X9               1 0.95207  0.1305
#> -------------------------------------------------
#>     Weight Estimate       LRT P-value    
#> X10      0.33195025 447.74998  <2e-16 ***
#> X11      0.33481743 457.78386  <2e-16 ***
#> X12      0.32787708 436.99374  <2e-16 ***
#> X13      0.00000000   0.00000  1.0000    
#> X14      0.00535525   0.12052  0.3497    
#> -------------------------------------------------
#> 
#> Estimates and Inference for Intercept and Adjusting Covariates
#> 
#>            Estimate      SE       Z P(Z > |z|)  
#> intercept 0.0403881 0.02048 1.97179   0.048633 *
#> 
#> Significance Codes: <0.001 '***' <0.01 '**' <0.05 '*' <0.10 '.' 
#> 
#> Total runtime for FGWQSR:  2.17 minutes on 10 cores.
```

<br> We can compare the true underlying group indices and chemical
weights with their estimates, along with if the correct inference is
made. Note that dashes are inserted for the true weights of group 2,
since the true underlying group index is null (thus, the weights are
unidentifiable).

``` r
group_index_frame = data.frame("True Group Index" = exp(gamma),
                          "Estimated Group Index" = exp(fgwqsr_fit$inference_frames$group_index_frame[,1]) %>% round(3),
                          "Signficiant?" = ifelse(fgwqsr_fit$inference_frames$group_index_frame[,3] < .05, "Yes", "No"), check.names = F)

rownames(group_index_frame) = rownames(fgwqsr_fit$inference_frames$group_index_frame)

true_weights = ess$weights %>% unlist %>% round(3)
true_weights[6:9] = "-"


weight_frame = data.frame("True Weight" = true_weights,
                          "Estimated Weight" = fgwqsr_fit$inference_frames$weight_frame[,1] %>% round(3),
                          "Signficiant?" = ifelse(fgwqsr_fit$inference_frames$weight_frame[,3] < .05, "Yes", "No"), check.names = F)

group_index_frame; weight_frame
#>                  True Group Index Estimated Group Index Signficiant?
#> Mixture Effect 1              0.5                 0.501          Yes
#> Mixture Effect 2              1.0                 0.995           No
#> Mixture Effect 3              1.5                 1.476          Yes
#>     True Weight Estimated Weight Signficiant?
#> w11       0.333            0.339          Yes
#> w12       0.333            0.315          Yes
#> w13       0.333            0.344          Yes
#> w14           0            0.002           No
#> w15           0            0.000           No
#> w21           -            0.000           No
#> w22           -            0.000           No
#> w23           -            0.000           No
#> w24           -            1.000           No
#> w31       0.333            0.332          Yes
#> w32       0.333            0.335          Yes
#> w33       0.333            0.328          Yes
#> w34           0            0.000           No
#> w35           0            0.005           No
```

<br>

# Fitting with Adjusting Continuous and Categorical Covariates

As mentioned earlier, the `/` character in the model formula sent to
`fgwqsr()` indicates that adjusting covariates follow this character.
Continuous characters can be referenced with their columnname from the
dataset, while categorical variables need to be referenced using the
prefix `i.`. We will illustrate this through an example where we create
a continuous `weight` variable and a categorical `city` variable. We
will give an effect (in log odds scale) of .5 to weight, .2 to the
comparison between city 2 and reference city 1, and -1 to the comparison
between city 3 and city 1.

<br>

``` r
set.seed(2023)
# create adjusting covariates
weight = rnorm(n = n, mean = 68, sd = 2.5)
city = sample(c("city_1", "city_2", "city_3"), size = n, replace = T)

# need to adjust intercept for case control ratio
intercept = -3.3

# create logit(pi) outcome WITH adjusting covariates
logit_pi =  gamma[1] * (chem_data[, 1:5] %*% ess$weights$w1) +
  gamma[2] * (chem_data[, 6:9] %*% ess$weights$w2) +
  gamma[3] * (chem_data[, 10:14] %*% ess$weights$w3) +
  .05 *weight +
  .2*ifelse(city == "city_2", 1, 0) + -1*ifelse(city == "city_3", 1, 0) +
  intercept

# transform to pi

pi = stats::plogis(logit_pi)

# transform to bernoulli outcome
y = sapply(pi, FUN = function(p) rbinom(1,1,p))

# create dataset

data = data.frame(y = y, chem_data, weight = weight, city = city)

head(data)
#>   y X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 X14   weight   city
#> 1 0  1  0  0  0  0  4  1  3  2   1   3   2   3   3 67.79054 city_3
#> 2 0  2  3  3  3  4  3  2  2  4   1   0   0   0   0 65.54264 city_3
#> 3 1  4  4  3  4  3  2  1  2  3   4   4   3   4   4 63.31233 city_1
#> 4 1  4  3  4  4  4  4  4  4  4   2   3   4   4   4 67.53464 city_1
#> 5 1  0  1  2  3  3  1  0  0  3   3   4   2   4   4 66.41629 city_2
#> 6 0  0  0  1  0  2  2  1  1  2   1   3   4   2   3 70.72699 city_3
```

<br> Now, we specify the formula to include the continuous and
categorical variables and call `fgwqsr()`.

``` r

mod_formula_adj = y ~ X1 + X2 + X3 + X4 + X5 | X6 + X7 + X8 + X9 | X10 + X11 + X12 + X13 + X14 / weight + i.city

fgwqsr_fit_adj = fgwqsr(formula = mod_formula_adj,
                    data = data,
                    quantiles = 5,
                    n_mvn_sims = 10000,
                    verbose = T)
#> 
#> Fitting full model and nested models...
#> 
#> Generating LRT distributions under H0...
```

<br> Again, we use the function `summary()` to view the results of the
`fgwqsr()` call.

``` r

summary(fgwqsr_fit_adj)
#> 
#> Call: 
#> FGWQSR with formula 'y ~ X1 + X2 + X3 + X4 + X5 | X6 + X7 + X8 + X9 | X10 + X11 + X12 + X13 + X14/weight + i.city' on n = 100000 observations.
#> 
#> 10000 samples used for simulated LRT distirbution.
#> 
#> Log Likelihood: -56252.41 | AIC: 112540.8 | BIC: 112712
#> 
#> Estimates and Inference for Group Index Effects
#> 
#>                     Estimate         LRT P-value    
#> Mixture Effect 1 -0.69003498 10742.65584  <2e-16 ***
#> Mixture Effect 2  0.00419649     0.65689  0.8652    
#> Mixture Effect 3  0.40065868  3700.09222  <2e-16 ***
#> 
#> Estimates and Inference for Weights
#> 
#>    Weight Estimate     LRT P-value    
#> X1      0.33127985 1318.55  <2e-16 ***
#> X2      0.32808914 1294.74  <2e-16 ***
#> X3      0.33071400 1312.28  <2e-16 ***
#> X4      0.00000000    0.00       1    
#> X5      0.00991701    0.00       1    
#> -------------------------------------------------
#>    Weight Estimate LRT P-value
#> X6               0   0       1
#> X7               0   0       1
#> X8               1   0       1
#> X9               0   0       1
#> -------------------------------------------------
#>     Weight Estimate     LRT P-value    
#> X10      0.31015452 388.187  <2e-16 ***
#> X11      0.32041925 415.890  <2e-16 ***
#> X12      0.36287563 530.574  <2e-16 ***
#> X13      0.00000000   0.000       1    
#> X14      0.00655059   0.000       1    
#> -------------------------------------------------
#> 
#> Estimates and Inference for Intercept and Adjusting Covariates
#> 
#>                Estimate      SE        Z P(Z > |z|)    
#> intercept   -3.09764720 0.19810 -15.6371 < 2.22e-16 ***
#> weight       0.04711361 0.00289  16.2875 < 2.22e-16 ***
#> city_city_2  0.20325289 0.01693  12.0037 < 2.22e-16 ***
#> city_city_3 -1.00098441 0.01857 -53.8973 < 2.22e-16 ***
#> 
#> Significance Codes: <0.001 '***' <0.01 '**' <0.05 '*' <0.10 '.' 
#> 
#> Total runtime for FGWQSR:  2.4 minutes on 10 cores.
```

<br>

Finally, we can compare the true parameter value with our estimates from
FGWQSR.

``` r

group_index_frame = data.frame("True Group Index" = exp(gamma),
                          "Estimated Group Index" = exp(fgwqsr_fit_adj$inference_frames$group_index_frame$Estimate) %>% round(3),
                          "Signficiant?" = ifelse(fgwqsr_fit_adj$inference_frames$group_index_frame$`P-value` < .05, "Yes", "No"), check.names = F)

rownames(group_index_frame) = rownames(fgwqsr_fit_adj$inference_frames$group_index_frame)

true_weights = ess$weights %>% unlist %>% round(3)
true_weights[6:9] = "-"


weight_frame = data.frame("True Weight" = true_weights,
                          "Estimated Weight" = fgwqsr_fit_adj$inference_frames$weight_frame$`Weight Estimate` %>% round(3),
                          "Signficiant?" = ifelse(fgwqsr_fit_adj$inference_frames$weight_frame$`P-value` < .05, "Yes", "No"), check.names = F)


adj_cov_frame = data.frame("True Covariate Effect" = c(-3.3, .5,.2,-1),
                           "Estimated Covariate Effect" = fgwqsr_fit_adj$inference_frames$adj_param_frame$Estimate,
                           "Significant?" = ifelse(fgwqsr_fit_adj$inference_frames$adj_param_frame$`P(Z > |z|)` < .05, "Yes", "No"), check.names = F)

rownames(adj_cov_frame) = rownames(fgwqsr_fit_adj$inference_frames$adj_param_frame)

group_index_frame; weight_frame; adj_cov_frame
#>                  True Group Index Estimated Group Index Signficiant?
#> Mixture Effect 1              0.5                 0.502          Yes
#> Mixture Effect 2              1.0                 1.004           No
#> Mixture Effect 3              1.5                 1.493          Yes
#>     True Weight Estimated Weight Signficiant?
#> w11       0.333            0.331          Yes
#> w12       0.333            0.328          Yes
#> w13       0.333            0.331          Yes
#> w14           0            0.000           No
#> w15           0            0.010           No
#> w21           -            0.000           No
#> w22           -            0.000           No
#> w23           -            1.000           No
#> w24           -            0.000           No
#> w31       0.333            0.310          Yes
#> w32       0.333            0.320          Yes
#> w33       0.333            0.363          Yes
#> w34           0            0.000           No
#> w35           0            0.007           No
#>             True Covariate Effect Estimated Covariate Effect Significant?
#> intercept                    -3.3                -3.09764720          Yes
#> weight                        0.5                 0.04711361          Yes
#> city_city_2                   0.2                 0.20325289          Yes
#> city_city_3                  -1.0                -1.00098441          Yes
```

# Fitting Models with BGWQSR

We also provide in this package functions to run BGWQSR models using the
runjags package. The `bgwqsr()` function takes a model formulation
similar to that in `fgwqsr()`, also utilizing the special characters
`|`, `/`, and `i.`. Fitting the BGWQSR model using runjags allows us to
leverage the “parallel” method in runjags that allows us to fit
independent mcmc chains on multiple cores, speeding up fitting time. The
`bgwqsr()` function takes several arguments:

- `formula` A formula for model fitting of BGWQSR. Please see
  description for formula construction

- `data` dataframe that contains all covariates and the outcome data.
  Column names of dataframe should match those referenced int he model
  formula.

- `quantiles` number of quantiles to quantize the exposure variables in
  the mixture portion of the model. Default value is 5.

- `n.iter` number of mcmc iterations after burnin and adapt iterations
  PER mcmc chain.

- `n.burnin` number of mcmc burnin samples PER mcmc chain

- `n.thin` thinning interval for mcmc samples PER mcmc chain

- `n.chains` number of separate independent mcmc chains to use.

- `n.adapt` number of mcmc samples to perform mcmc adaption PER mcmc
  chain.

- `inits` initial values to provide for prior distributions.

- `method` method for mcmc fitting, a passed argument to run.jags
  function. Can be one of: `rjags`, `simple`, `interruptible`,
  `parallel`, `rjparallel`, `background`, `bgparallel`, or `snow`.

``` r

bgwqsr_fit = bgwqsr(formula = mod_formula_adj,
                    data = data,
                    quantiles = 5,
                    n.iter = 5000,
                    n.adapt = 1000,
                    n.burnin = 4000,
                    n.thin = 1, 
                    n.chains = 3,
                    method = "parallel")
#> Calling 3 simulations using the parallel method...
#> Following the progress of chain 1 (the program will wait for all chains
#> to finish before continuing):
#> Welcome to JAGS 4.3.2 (official binary) on Thu Aug 17 23:37:44 2023
#> JAGS is free software and comes with ABSOLUTELY NO WARRANTY
#> Loading module: basemod: ok
#> Loading module: bugs: ok
#> . Loading module: glm: ok
#> . . Reading data file data.txt
#> . Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 100000
#>    Unobserved stochastic nodes: 17
#>    Total graph size: 2113800
#> . Reading parameter file inits1.txt
#> . Initializing model
#> . Adapting 1000
#> -------------------------------------------------| 1000
#> ++++++++++++++++++++++++++++++++++++++++++++++++++ 100%
#> Adaptation incomplete.
#> . Updating 4000
#> -------------------------------------------------| 4000
#> ************************************************** 100%
#> . . . . . . . . . . . . . . . . . . . . . . Updating 5000
#> -------------------------------------------------| 5000
#> ************************************************** 100%
#> . . . . Updating 0
#> . Deleting model
#> All chains have finished
#> Simulation complete.  Reading coda files...
#> Coda files loaded successfully
#> Calculating summary statistics...
#> Calculating the Gelman-Rubin statistic for 21 variables....
#> Finished running the simulation
#> Calculating summary statistics...
#> Calculating the Gelman-Rubin statistic for 21 variables....
```

<br>

To retrieve the summary measures from the MCMC sampling, we can extract
the summaries object from the bgwqsr model.

``` r
bgwqsr_fit$model$summaries
#>                      Lower95       Median     Upper95         Mean          SD
#> B0              -3.43820e+00 -3.052330000 -2.65882000 -3.051973143 0.199682261
#> B1              -7.00969e-01 -0.686971000 -0.67276500 -0.687058138 0.007191698
#> B2              -2.39051e-02 -0.010599850  0.00390657 -0.010115060 0.007132723
#> B3               3.89593e-01  0.403025500  0.41640300  0.402976373 0.006837825
#> phi_weight       4.08947e-02  0.046712700  0.05226110  0.046722540 0.002915755
#> phi_city_city_2  1.67072e-01  0.200891500  0.23296100  0.200893781 0.016912649
#> phi_city_city_3 -1.03775e+00 -1.001410000 -0.96558700 -1.001556783 0.018491809
#> w1[1]            3.16707e-01  0.333155000  0.34927300  0.333250088 0.008489713
#> w1[2]            3.12449e-01  0.330164000  0.34621400  0.330204740 0.008439097
#> w1[3]            3.16598e-01  0.332582000  0.34901300  0.332341947 0.008310167
#> w1[4]            1.83418e-05  0.000592886  0.01025460  0.002293306 0.003717927
#> w1[5]            3.47763e-06  0.000821073  0.00958184  0.001909891 0.003185331
#> w2[1]            8.88496e-14  0.028174700  0.55762700  0.118956058 0.193476807
#> w2[2]            1.23444e-02  0.734034500  0.99999900  0.625282603 0.325737390
#> w2[3]            6.74466e-11  0.017633900  0.57631000  0.100166899 0.191278378
#> w2[4]            1.73558e-09  0.041523100  0.70487400  0.155594444 0.229513899
#> w3[1]            2.79829e-01  0.306624000  0.33599900  0.306916503 0.014570448
#> w3[2]            2.90683e-01  0.318362500  0.34547800  0.318432012 0.014075568
#> w3[3]            3.33094e-01  0.360208000  0.38876400  0.360050526 0.014534736
#> w3[4]            1.94065e-05  0.003501060  0.01463940  0.004894788 0.004826488
#> w3[5]            3.00510e-04  0.006384550  0.02742980  0.009706163 0.008848876
#>                 Mode        MCerr MC%ofSD SSeff        AC.10     psrf
#> B0                NA 2.236840e-03     1.1  7969  0.014276549 1.000419
#> B1                NA 8.622673e-05     1.2  6956  0.013129134 1.000188
#> B2                NA 2.222310e-04     3.1  1030  0.212066097 1.003390
#> B3                NA 8.057297e-05     1.2  7202  0.005349031 1.000261
#> phi_weight        NA 3.263117e-05     1.1  7984  0.012684100 1.000341
#> phi_city_city_2   NA 1.788536e-04     1.1  8942  0.006856053 1.000016
#> phi_city_city_3   NA 2.119176e-04     1.1  7614 -0.006321003 1.000217
#> w1[1]             NA 2.093849e-04     2.5  1644  0.124518169 1.007972
#> w1[2]             NA 2.173361e-04     2.6  1508  0.133670261 1.006177
#> w1[3]             NA 2.095775e-04     2.5  1572  0.120411922 1.004850
#> w1[4]             NA 9.218974e-04    24.8    16  0.981539119 2.642258
#> w1[5]             NA 5.270906e-04    16.5    37  0.952410466 1.418734
#> w2[1]             NA 9.615982e-03     5.0   405  0.504547653 1.046622
#> w2[2]             NA 1.691102e-02     5.2   371  0.609044469 1.022007
#> w2[3]             NA 1.030457e-02     5.4   345  0.596020164 1.064840
#> w2[4]             NA 1.082843e-02     4.7   449  0.495779989 1.014941
#> w3[1]             NA 3.956060e-04     2.7  1357  0.179787503 1.001499
#> w3[2]             NA 3.462830e-04     2.5  1652  0.143734543 1.002254
#> w3[3]             NA 3.775704e-04     2.6  1482  0.149803981 1.000577
#> w3[4]             NA 7.976633e-04    16.5    37  0.944523065 1.217704
#> w3[5]             NA 1.423914e-03    16.1    39  0.949863628 1.022293
```

We can analyze the mixing of the markov chains and corresponding
autocorrelation plots using functions from the coda package. We will
output only the first few plots.

``` r
par(mfrow = c(3,3))
coda::traceplot(bgwqsr_fit$model$mcmc) # traceplot
```

<img src="man/figures/README-unnamed-chunk-18-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-18-2.png" width="100%" /><img src="man/figures/README-unnamed-chunk-18-3.png" width="100%" />

``` r
coda::autocorr.plot(bgwqsr_fit$model$mcmc, ask = F) # autocorrelation plot
```

<img src="man/figures/README-unnamed-chunk-19-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-19-2.png" width="100%" /><img src="man/figures/README-unnamed-chunk-19-3.png" width="100%" /><img src="man/figures/README-unnamed-chunk-19-4.png" width="100%" /><img src="man/figures/README-unnamed-chunk-19-5.png" width="100%" /><img src="man/figures/README-unnamed-chunk-19-6.png" width="100%" /><img src="man/figures/README-unnamed-chunk-19-7.png" width="100%" />

``` r
par(mfrow = c(3,3))
coda::densplot(bgwqsr_fit$model$mcmc) # posterior density plots
```

<img src="man/figures/README-unnamed-chunk-20-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-20-2.png" width="100%" /><img src="man/figures/README-unnamed-chunk-20-3.png" width="100%" />

``` r
plot(bgwqsr_fit$model$mcmc) # combines traceplot() and densplot()
```

<img src="man/figures/README-unnamed-chunk-21-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-21-2.png" width="100%" /><img src="man/figures/README-unnamed-chunk-21-3.png" width="100%" /><img src="man/figures/README-unnamed-chunk-21-4.png" width="100%" /><img src="man/figures/README-unnamed-chunk-21-5.png" width="100%" /><img src="man/figures/README-unnamed-chunk-21-6.png" width="100%" />

Finally, we can look at plots for the posterior credible intervals for
group indices, single chemical weights, and corresponding posterior
means. To do this, we can use the functions `plot_weights()`,
`plot_betas()`, and `plot_result()`. `plot_results()` combines both
plots generates by plot_weights and plot_betas into a side by side
figure. Below are code examples for each of the three function calls.

``` r
plot_betas(bgwqsr_fit)
```

<img src="man/figures/README-unnamed-chunk-22-1.png" width="100%" />

``` r
plot_weights(bgwqsr_fit)
```

<img src="man/figures/README-unnamed-chunk-23-1.png" width="100%" />

``` r
plot_results(bgwqsr_fit)
```

<img src="man/figures/README-unnamed-chunk-24-1.png" width="100%" />
