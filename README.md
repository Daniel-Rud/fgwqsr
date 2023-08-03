
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
$$\sum_{k = 1}^{c_g} w_{g,k} =1, \ w_{g,k} \in (0,1)$$

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

- `|` - denotes the boundary of a mixture group, used to seperate
  chemicals within a mixture group.

- `/` - denotes the end of the mixture group specification, adjusting
  covariates can be added to the formula after this character. If no
  adjusting covariates, do not need to specify.

- `i.` - precedes categorical variables to denote a categorical
  variable. For example, if we have a categorical variable cat_var, we
  would denote this in the model formula by i.cat_var. This is similar
  to the stata syntax to declare categorical variables.

<br> The fgwqsr() function has other options too: <br>

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

Notice that we did not use the “$/$” character since we did not include
adjusting covariates. Now, we can fit the model using the function
fgwqsr().

``` r
fgwqsr_fit = fgwqsr(formula = mod_formula,
                    data = data,
                    quantiles = 5,
                    n_mvn_sims = 100,
                    verbose = T)
#> 
#> Fitting full model and nested models...
#> 
#> Generating LRT distributions under H0...
```

<br>

We can see the model summary using the call summary()

``` r

summary(fgwqsr_fit)
#> 
#> Call: 
#> FGWQSR with formula 'y ~ X1 + X2 + X3 + X4 + X5 | X6 + X7 + X8 + X9 | X10 + X11 + X12 + X13 + X14' on n = 10000 observations.
#> 
#> 100 samples used for simulated LRT distirbution.
#> 
#> Log Likelihood: -5960.159 | AIC: 11950.32 | BIC: 12058.47
#> 
#> Estimates and Inference for Group Index Effects
#> 
#>                    Estimate        LRT P-value    
#> Mixture Effect 1 -0.7023870 1168.18784  <2e-16 ***
#> Mixture Effect 2  0.0246758    1.49567    0.57    
#> Mixture Effect 3  0.4113445  408.88836  <2e-16 ***
#> 
#> Estimates and Inference for Weights
#> 
#>    Weight Estimate       LRT P-value    
#> X1       0.3055009 128.10791  <2e-16 ***
#> X2       0.3364256 156.35067  <2e-16 ***
#> X3       0.3223592 142.16491  <2e-16 ***
#> X4       0.0236343   0.75019    0.21    
#> X5       0.0120799   0.20227    0.35    
#> -------------------------------------------------
#>    Weight Estimate     LRT P-value
#> X6        0.248205 0.10221    0.42
#> X7        0.207174 0.07197    0.29
#> X8        0.111892 0.02121    0.47
#> X9        0.432729 0.30953    0.26
#> -------------------------------------------------
#>     Weight Estimate      LRT P-value    
#> X10      0.24593511 26.93181  <2e-16 ***
#> X11      0.31020577 43.50811  <2e-16 ***
#> X12      0.37606638 62.98216  <2e-16 ***
#> X13      0.00126201  0.00072    0.45    
#> X14      0.06653073  2.05122    0.09 .  
#> -------------------------------------------------
#> 
#> Estimates and Inference for Intercept and Adjusting Covariates
#> 
#>             Estimate     SE         Z P(Z > |z|)
#> intercept -0.0129252 0.0646 -0.200085    0.84141
#> 
#> Significance Codes: <0.001 '***' <0.01 '**' <0.05 '*' <0.10 '.' 
#> 
#> Total runtime for FGWQSR:  1.13 minutes on 10 cores.
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
#> Mixture Effect 1              0.5                 0.495          Yes
#> Mixture Effect 2              1.0                 1.025           No
#> Mixture Effect 3              1.5                 1.509          Yes
#>     True Weight Estimated Weight Signficiant?
#> w11       0.333            0.306          Yes
#> w12       0.333            0.336          Yes
#> w13       0.333            0.322          Yes
#> w14           0            0.024           No
#> w15           0            0.012           No
#> w21           -            0.248           No
#> w22           -            0.207           No
#> w23           -            0.112           No
#> w24           -            0.433           No
#> w31       0.333            0.246          Yes
#> w32       0.333            0.310          Yes
#> w33       0.333            0.376          Yes
#> w34           0            0.001           No
#> w35           0            0.067           No
```

<br>

# Fitting with Adjusting Continuous and Categorical Covariates

As mentioned earlier, the “$/$” character in the model formula sent to
fgwqsr() indicates that adjusting covariates follow this character.
Continuous characters can be referenced with their columnname from the
dataset, while categorical variables need to be referenced using the
prefix “i.”. We will illustrate this through an example where we create
a continuous $weight$ variable and a categorical $city$ variable. We
will give an effect (in log odds scale) of .5 to weight, .2 to the
comparison between city 2 and reference city 1, and -1 to the comparison
between city 3 and city 1.

<br>

``` r
set.seed(1)
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
#> 1 1  1  2  1  4  3  1  3  3  3   1   4   2   1   0 66.43387 city_3
#> 2 0  4  3  3  4  4  4  4  4  3   0   2   1   1   0 68.45911 city_3
#> 3 0  2  3  4  2  3  1  0  1  1   2   4   3   2   2 65.91093 city_3
#> 4 0  3  3  1  1  3  4  3  4  3   2   3   1   4   4 71.98820 city_3
#> 5 1  2  1  3  2  4  2  3  2  1   3   1   4   3   4 68.82377 city_1
#> 6 0  2  0  2  0  0  2  1  2  2   1   1   1   3   0 65.94883 city_3
```

<br> Now, we specify the formula to include the continuous and
categorical variables and call fgwqsr().

``` r

mod_formula_adj = y ~ X1 + X2 + X3 + X4 + X5 | X6 + X7 + X8 + X9 | X10 + X11 + X12 + X13 + X14 / weight + i.city

fgwqsr_fit_adj = fgwqsr(formula = mod_formula_adj,
                    data = data,
                    quantiles = 5,
                    n_mvn_sims = 100,
                    verbose = T)
#> 
#> Fitting full model and nested models...
#> 
#> Generating LRT distributions under H0...
```

<br> Again, we use the function summary() to view the results of the
fgwqsr() call.

``` r

summary(fgwqsr_fit_adj)
#> 
#> Call: 
#> FGWQSR with formula 'y ~ X1 + X2 + X3 + X4 + X5 | X6 + X7 + X8 + X9 | X10 + X11 + X12 + X13 + X14/weight + i.city' on n = 10000 observations.
#> 
#> 100 samples used for simulated LRT distirbution.
#> 
#> Log Likelihood: -5556.169 | AIC: 11148.34 | BIC: 11278.12
#> 
#> Estimates and Inference for Group Index Effects
#> 
#>                    Estimate       LRT P-value    
#> Mixture Effect 1 -0.7219736 1153.2447  <2e-16 ***
#> Mixture Effect 2  0.0708484   13.1311    0.01 ** 
#> Mixture Effect 3  0.4037247  375.4276  <2e-16 ***
#> 
#> Estimates and Inference for Weights
#> 
#>    Weight Estimate       LRT P-value    
#> X1      0.36719961 178.32227  <2e-16 ***
#> X2      0.28786631 110.03380  <2e-16 ***
#> X3      0.33762099 150.47260  <2e-16 ***
#> X4      0.00317254   0.01309    0.54    
#> X5      0.00414055   0.02615    0.43    
#> -------------------------------------------------
#>    Weight Estimate     LRT P-value  
#> X6       0.0279733 0.00981    0.53  
#> X7       0.6527946 5.39774    0.02 *
#> X8       0.1479473 0.28094    0.33  
#> X9       0.1712847 0.36739    0.24  
#> -------------------------------------------------
#>     Weight Estimate     LRT P-value    
#> X10        0.318641 42.0794  <2e-16 ***
#> X11        0.330624 45.7642  <2e-16 ***
#> X12        0.350735 50.9094  <2e-16 ***
#> X13        0.000000  0.0000    0.74    
#> X14        0.000000  0.0000    0.78    
#> -------------------------------------------------
#> 
#> Estimates and Inference for Intercept and Adjusting Covariates
#> 
#>                Estimate      SE         Z P(Z > |z|)    
#> intercept   -2.66161146 0.62443  -4.26245 2.0220e-05 ***
#> weight       0.03962030 0.00913   4.34090 1.4190e-05 ***
#> city_city_2  0.21681253 0.05404   4.01191 6.0229e-05 ***
#> city_city_3 -1.03548465 0.05909 -17.52474 < 2.22e-16 ***
#> 
#> Significance Codes: <0.001 '***' <0.01 '**' <0.05 '*' <0.10 '.' 
#> 
#> Total runtime for FGWQSR:  1.11 minutes on 10 cores.
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
#> Mixture Effect 1              0.5                 0.486          Yes
#> Mixture Effect 2              1.0                 1.073          Yes
#> Mixture Effect 3              1.5                 1.497          Yes
#>     True Weight Estimated Weight Signficiant?
#> w11       0.333            0.367          Yes
#> w12       0.333            0.288          Yes
#> w13       0.333            0.338          Yes
#> w14           0            0.003           No
#> w15           0            0.004           No
#> w21           -            0.028           No
#> w22           -            0.653          Yes
#> w23           -            0.148           No
#> w24           -            0.171           No
#> w31       0.333            0.319          Yes
#> w32       0.333            0.331          Yes
#> w33       0.333            0.351          Yes
#> w34           0            0.000           No
#> w35           0            0.000           No
#>             True Covariate Effect Estimated Covariate Effect Significant?
#> intercept                    -3.3                 -2.6616115          Yes
#> weight                        0.5                  0.0396203          Yes
#> city_city_2                   0.2                  0.2168125          Yes
#> city_city_3                  -1.0                 -1.0354846          Yes
```

# Fitting Models with BGWQSR

We also provide in this package functions to run BGWQSR models using the
runjags package. The bgwqsr() function takes a model formulation similar
to that in fgwqsr(), also utilizing the special characters “$|$”, “$/$”,
and “i.”. Fitting the BGWQSR model using runjags allows us to leverage
the “parallel” method in runjags that allows us to fit independent mcmc
chains on multiple cores, speeding up fitting time. The bgwqsr()
function takes several arguments:

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
  function. Can be one of: `rjags`, `simple`, `interruptible`,’
  `parallel`, `rjparallel`, `background`, `bgparallel`, or `snow`.

``` r

bgwqsr_fit = bgwqsr(formula = mod_formula_adj,
                    data = data,
                    quantiles = 5,
                    n.iter = 1000,
                    n.burnin = 40,
                    n.thin = 1, n.chains = 3,
                    n.adapt = 10,
                    method = "parallel")
#> Calling 3 simulations using the parallel method...
#> Following the progress of chain 1 (the program will wait for all chains
#> to finish before continuing):
#> Welcome to JAGS 4.3.2 (official binary) on Thu Aug  3 15:33:11 2023
#> JAGS is free software and comes with ABSOLUTELY NO WARRANTY
#> Loading module: basemod: ok
#> Loading module: bugs: ok
#> . Loading module: glm: ok
#> . . Reading data file data.txt
#> . Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 10000
#>    Unobserved stochastic nodes: 17
#>    Total graph size: 220906
#> . Reading parameter file inits1.txt
#> . Initializing model
#> . Adapting 10
#> Adaptation incomplete
#> . Updating 40
#> . . . . . . . . . . . . . . . . . . . . . . Updating 1000
#> -------------------------------------------------| 1000
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
#>                      Lower95      Median    Upper95        Mean          SD
#> B0              -3.41273e+00 -2.12428000 -0.9028630 -2.13380252 0.661837747
#> B1              -7.71074e-01 -0.72503650 -0.6806060 -0.72537830 0.023587364
#> B2               1.96039e-02  0.06073335  0.0999250  0.06062689 0.020883301
#> B3               3.65088e-01  0.40809650  0.4499350  0.40775981 0.021391236
#> phi_weight       1.46445e-02  0.03222635  0.0519251  0.03237327 0.009732472
#> phi_city_city_2  6.25986e-02  0.18485350  0.2929000  0.18382036 0.059122404
#> phi_city_city_3 -1.16384e+00 -1.04874500 -0.9312250 -1.04846931 0.059948509
#> w1[1]            3.10031e-01  0.35901100  0.4100450  0.35728153 0.025613528
#> w1[2]            2.31327e-01  0.27804900  0.3348240  0.27797406 0.026387614
#> w1[3]            2.73779e-01  0.32544600  0.3748780  0.32631108 0.025981286
#> w1[4]            4.75491e-05  0.01013670  0.0441801  0.01387641 0.014103587
#> w1[5]            2.45905e-04  0.02451640  0.0492053  0.02455690 0.014961278
#> w2[1]            7.37600e-04  0.13951000  0.4945270  0.18025597 0.161376595
#> w2[2]            1.49478e-04  0.45208800  0.9470420  0.42178710 0.333197559
#> w2[3]            9.02283e-05  0.06863615  0.5864340  0.14638910 0.189965132
#> w2[4]            1.51517e-03  0.18338550  0.7731790  0.25156782 0.233090798
#> w3[1]            2.11589e-01  0.30103150  0.3898180  0.29923034 0.047799525
#> w3[2]            2.29993e-01  0.31922250  0.4072830  0.31770958 0.046368542
#> w3[3]            2.51638e-01  0.34019600  0.4347630  0.34043288 0.047018854
#> w3[4]            3.39058e-04  0.01415870  0.0555459  0.01990778 0.018339622
#> w3[5]            1.25381e-04  0.01798075  0.0612273  0.02271944 0.019595252
#>                 Mode        MCerr MC%ofSD SSeff      AC.10     psrf
#> B0                NA 0.0366659158     5.5   326 0.16542471 1.058841
#> B1                NA 0.0006744405     2.9  1223 0.02434819 1.002149
#> B2                NA 0.0006137575     2.9  1158 0.04445830 1.117741
#> B3                NA 0.0005643491     2.6  1437 0.03361305 1.003649
#> phi_weight        NA 0.0005587788     5.7   303 0.16757459 1.066254
#> phi_city_city_2   NA 0.0023645464     4.0   625 0.07424356 1.005846
#> phi_city_city_3   NA 0.0017177313     2.9  1218 0.01416054 1.001210
#> w1[1]             NA 0.0013703750     5.4   349 0.13302945 1.004929
#> w1[2]             NA 0.0014779093     5.6   319 0.14283951 1.018734
#> w1[3]             NA 0.0016191259     6.2   257 0.16131391 1.028073
#> w1[4]             NA 0.0025918484    18.4    30 0.78934898 1.979973
#> w1[5]             NA 0.0027030364    18.1    31 0.80974423 1.260754
#> w2[1]             NA 0.0264491639    16.4    37 0.80478535 1.225026
#> w2[2]             NA 0.0577072214    17.3    33 0.78715046 3.535312
#> w2[3]             NA 0.0428085624    22.5    20 0.85828129 1.710840
#> w2[4]             NA 0.0405026270    17.4    33 0.82571197 1.659761
#> w3[1]             NA 0.0029453542     6.2   263 0.21781170 1.005608
#> w3[2]             NA 0.0028011142     6.0   274 0.17817684 1.002106
#> w3[3]             NA 0.0030271245     6.4   241 0.25073098 1.002907
#> w3[4]             NA 0.0034971990    19.1    28 0.82861071 1.424545
#> w3[5]             NA 0.0032998754    16.8    35 0.79147539 1.423923
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
means. To do this, we can use the functions plot_weights(),
plot_betas(), and plot_result(). plot_results() combines both plots
generates by plot_weights and plot_betas into a side by side figure.
Below are code examples for each of the three function calls.

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
