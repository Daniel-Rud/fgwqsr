
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fgwqsr – Frequentist Grouped Weighted Quantile Sum Regression <img src="man/figures/FGWQSR_hex.png" width="100" align="right" />

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
any data splitting

<br>

# Installation

You can install the development version of fgwqsr from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Daniel-Rud/fgwqsr")
#> Using github PAT from envvar GITHUB_PAT. Use `gitcreds::gitcreds_set()` and unset GITHUB_PAT in .Renviron (or elsewhere) if you want to use the more secure git credential store instead.
#> Downloading GitHub repo Daniel-Rud/fgwqsr@HEAD
#> 
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>      checking for file ‘/private/var/folders/b4/9kg7p6cj729_pzc5dggk_9jm0000gn/T/RtmpNaCE4k/remotes17c3f61554d60/Daniel-Rud-fgwqsr-52f0b8b/DESCRIPTION’ ...  ✔  checking for file ‘/private/var/folders/b4/9kg7p6cj729_pzc5dggk_9jm0000gn/T/RtmpNaCE4k/remotes17c3f61554d60/Daniel-Rud-fgwqsr-52f0b8b/DESCRIPTION’
#>   ─  preparing ‘fgwqsr’:
#>    checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
#>   ─  checking for LF line-endings in source and make files and shell scripts
#>   ─  checking for empty or unneeded directories
#>      Omitted ‘LazyData’ from DESCRIPTION
#>   ─  building ‘fgwqsr_0.1.0.tar.gz’
#>      
#> 
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
set.seed(11) # for reproducibility
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
#> 1  0  1  0  0  3  1  4  3  2   0   0   0   0   0
#> 2  0  1  1  3  1  0  0  1  1   2   2   2   1   1
#> 3  1  0  0  0  2  0  0  2  1   3   1   1   0   1
#> 4  3  0  3  1  0  4  4  3  4   2   1   3   2   3
#> 5  2  3  4  3  1  4  2  0  3   4   1   1   1   3
#> 6  3  2  1  3  4  2  1  2  0   2   2   2   1   0
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
#> 1 0  0  1  0  0  3  1  4  3  2   0   0   0   0   0
#> 2 1  0  1  1  3  1  0  0  1  1   2   2   2   1   1
#> 3 0  1  0  0  0  2  0  0  2  1   3   1   1   0   1
#> 4 0  3  0  3  1  0  4  4  3  4   2   1   3   2   3
#> 5 0  2  3  4  3  1  4  2  0  3   4   1   1   1   3
#> 6 1  3  2  1  3  4  2  1  2  0   2   2   2   1   0
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

- `family` - can be one of ‘binomial’, ‘gaussian’ or ‘poisson’ for
  binary, continuous, and count outcomes respectivley.

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
                    family = 'binomial',
                    n_mvn_sims = 10000,
                    verbose = T)
#> Fitting full and nested FGWQSR models...
#> 
#> Generating LRT distributions under H0...
```

<br>

We can see the model summary using the call `summary()`

``` r

summary(fgwqsr_fit)
#> 
#> Call: 
#> FGWQSR with formula 'y ~ X1 + X2 + X3 + X4 + X5 | X6 + X7 + X8 + X9 | X10 + X11 + X12 + X13 + X14' on n = 10000 observations and family = 'binomial'.
#> 
#> 10000 samples used for simulated LRT distirbution.
#> 
#> Log Likelihood: -5997.012 | AIC: 12024.02 | BIC: 12132.18
#> 
#> Estimates and Inference for Group Index Effects
#> 
#>                    Estimate        LRT P-value    
#> Mixture Effect 1 -0.6624635 1057.49341  <2e-16 ***
#> Mixture Effect 2  0.0195867    1.55513  0.5966    
#> Mixture Effect 3  0.4131177  431.13990  <2e-16 ***
#> 
#> Estimates and Inference for Weights
#> 
#>    Weight Estimate       LRT P-value    
#> X1       0.2999211 112.43088  <2e-16 ***
#> X2       0.3417464 143.95821  <2e-16 ***
#> X3       0.3225850 128.82142  <2e-16 ***
#> X4       0.0357475   1.58066  0.0867 .  
#> X5       0.0000000   0.00000  1.0000    
#> -------------------------------------------------
#>    Weight Estimate     LRT P-value  
#> X6               0 0.00000  1.0000  
#> X7               0 0.00000  1.0000  
#> X8               0 0.00000  1.0000  
#> X9               1 1.45345  0.0828 .
#> -------------------------------------------------
#>     Weight Estimate     LRT P-value    
#> X10        0.339334 54.3167  <2e-16 ***
#> X11        0.347920 57.0600  <2e-16 ***
#> X12        0.312745 46.7864  <2e-16 ***
#> X13        0.000000  0.0000       1    
#> X14        0.000000  0.0000       1    
#> -------------------------------------------------
#> 
#> Estimates and Inference for Intercept and Adjusting Covariates
#> 
#>             Estimate      SE        Z P(Z > |z|)
#> intercept -0.0866543 0.06543 -1.32441    0.18537
#> 
#> Significance Codes: <0.001 '***' <0.01 '**' <0.05 '*' <0.10 '.' 
#> 
#> Total runtime for FGWQSR:  1.67 minutes on 10 cores.
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
#> Mixture Effect 1              0.5                 0.516          Yes
#> Mixture Effect 2              1.0                 1.020           No
#> Mixture Effect 3              1.5                 1.512          Yes
#>     True Weight Estimated Weight Signficiant?
#> w11       0.333            0.300          Yes
#> w12       0.333            0.342          Yes
#> w13       0.333            0.323          Yes
#> w14           0            0.036           No
#> w15           0            0.000           No
#> w21           -            0.000           No
#> w22           -            0.000           No
#> w23           -            0.000           No
#> w24           -            1.000           No
#> w31       0.333            0.339          Yes
#> w32       0.333            0.348          Yes
#> w33       0.333            0.313          Yes
#> w34           0            0.000           No
#> w35           0            0.000           No
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
set.seed(11)
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
#> 1 1  0  1  0  0  3  1  4  3  2   0   0   0   0   0 66.52242 city_1
#> 2 1  0  1  1  3  1  0  0  1  1   2   2   2   1   1 68.06649 city_1
#> 3 1  1  0  0  0  2  0  0  2  1   3   1   1   0   1 64.20862 city_3
#> 4 1  3  0  3  1  0  4  4  3  4   2   1   3   2   3 64.59337 city_2
#> 5 0  2  3  4  3  1  4  2  0  3   4   1   1   1   3 70.94622 city_3
#> 6 1  3  2  1  3  4  2  1  2  0   2   2   2   1   0 65.66462 city_1
```

<br> Now, we specify the formula to include the continuous and
categorical variables and call `fgwqsr()`.

``` r

mod_formula_adj = y ~ X1 + X2 + X3 + X4 + X5 | X6 + X7 + X8 + X9 | X10 + X11 + X12 + X13 + X14 / weight + i.city

fgwqsr_fit_adj = fgwqsr(formula = mod_formula_adj,
                    data = data,
                    quantiles = 5,
                    family = 'binomial',
                    n_mvn_sims = 10000,
                    verbose = T)
#> Fitting full and nested FGWQSR models...
#> 
#> Generating LRT distributions under H0...
```

<br> Again, we use the function `summary()` to view the results of the
`fgwqsr()` call.

``` r

summary(fgwqsr_fit_adj)
#> 
#> Call: 
#> FGWQSR with formula 'y ~ X1 + X2 + X3 + X4 + X5 | X6 + X7 + X8 + X9 | X10 + X11 + X12 + X13 + X14/weight + i.city' on n = 10000 observations and family = 'binomial'.
#> 
#> 10000 samples used for simulated LRT distirbution.
#> 
#> Log Likelihood: -5584.27 | AIC: 11204.54 | BIC: 11334.33
#> 
#> Estimates and Inference for Group Index Effects
#> 
#>                    Estimate       LRT P-value    
#> Mixture Effect 1 -0.7183223 1155.8522  <2e-16 ***
#> Mixture Effect 2  0.0259615    2.2956  0.4172    
#> Mixture Effect 3  0.3965223  349.5661  <2e-16 ***
#> 
#> Estimates and Inference for Weights
#> 
#>    Weight Estimate     LRT P-value    
#> X1        0.352231 172.232  <2e-16 ***
#> X2        0.336441 155.602  <2e-16 ***
#> X3        0.311328 133.957  <2e-16 ***
#> X4        0.000000   0.000       1    
#> X5        0.000000   0.000       1    
#> -------------------------------------------------
#>    Weight Estimate     LRT P-value  
#> X6       0.0817587 0.01341  0.3686  
#> X7       0.0000000 0.00000  1.0000  
#> X8       0.0000000 0.00000  1.0000  
#> X9       0.9182413 1.60019  0.0652 .
#> -------------------------------------------------
#>     Weight Estimate      LRT P-value    
#> X10       0.3542949 50.28806  <2e-16 ***
#> X11       0.2789109 31.16497  <2e-16 ***
#> X12       0.2988009 35.97137  <2e-16 ***
#> X13       0.0679932  1.88293  0.0739 .  
#> X14       0.0000000  0.00000  1.0000    
#> -------------------------------------------------
#> 
#> Estimates and Inference for Intercept and Adjusting Covariates
#> 
#>                Estimate      SE         Z P(Z > |z|)    
#> intercept   -3.42824833 0.63539  -5.39553 6.8321e-08 ***
#> weight       0.05235603 0.00927   5.64542 1.6478e-08 ***
#> city_city_2  0.20128467 0.05345   3.76578 0.00016603 ***
#> city_city_3 -1.03218054 0.05940 -17.37727 < 2.22e-16 ***
#> 
#> Significance Codes: <0.001 '***' <0.01 '**' <0.05 '*' <0.10 '.' 
#> 
#> Total runtime for FGWQSR:  1.77 minutes on 10 cores.
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
#> Mixture Effect 1              0.5                 0.488          Yes
#> Mixture Effect 2              1.0                 1.026           No
#> Mixture Effect 3              1.5                 1.487          Yes
#>     True Weight Estimated Weight Signficiant?
#> w11       0.333            0.352          Yes
#> w12       0.333            0.336          Yes
#> w13       0.333            0.311          Yes
#> w14           0            0.000           No
#> w15           0            0.000           No
#> w21           -            0.082           No
#> w22           -            0.000           No
#> w23           -            0.000           No
#> w24           -            0.918           No
#> w31       0.333            0.354          Yes
#> w32       0.333            0.279          Yes
#> w33       0.333            0.299          Yes
#> w34           0            0.068           No
#> w35           0            0.000           No
#>             True Covariate Effect Estimated Covariate Effect Significant?
#> intercept                    -3.3                -3.42824833          Yes
#> weight                        0.5                 0.05235603          Yes
#> city_city_2                   0.2                 0.20128467          Yes
#> city_city_3                  -1.0                -1.03218054          Yes
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
#> Warning in system(paste0("lipo -archs ", gsub("/bin/jags",
#> "/libexec/jags-terminal", : running command 'lipo -archs
#> /usr/local/libexec/jags-terminal' had status 1
#> Warning in runjags::testjags(silent = T): Reading the binary type of JAGS
#> failed - do you need to update the Xcode Command Line Tools?
#> Calling 3 simulations using the parallel method...
#> Following the progress of chain 1 (the program will wait for all chains
#> to finish before continuing):
#> Welcome to JAGS 4.3.2 (official binary) on Wed Apr 23 17:30:28 2025
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
#>    Total graph size: 220912
#> . Reading parameter file inits1.txt
#> . Initializing model
#> . Adapting 1000
#> -------------------------------------------------| 1000
#> ++++++++++++++++++++++++++++++++++++++++++++++++++ 100%
#> Adaptation successful
#> . Updating 4000
#> -------------------------------------------------| 4000
#> ************************************************** 100%
#> . . . . . . . . . . . . . . . . . . . . . . Updating 5000
#> -------------------------------------------------| 5000
#> ************************************************** 100%
#> . . . . Updating 0
#> . Deleting model
#> . 
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
#>                      Lower95        Median     Upper95          Mean
#> B0              -4.34782e+00 -3.0878650000 -1.79175000 -3.0847462334
#> B1              -7.60202e-01 -0.7156250000 -0.67162100 -0.7156927088
#> B2              -1.95092e-02  0.0161372500  0.05685000  0.0167256253
#> B3               3.49472e-01  0.3920905000  0.43472400  0.3922761100
#> phi_weight       2.84115e-02  0.0479448000  0.06594520  0.0479189677
#> phi_city_city_2  1.30874e-02  0.1619150000  0.27285700  0.1579288804
#> phi_city_city_3 -1.17209e+00 -1.0464550000 -0.92675000 -1.0466077578
#> w1[1]            3.04020e-01  0.3514860000  0.40182500  0.3517599415
#> w1[2]            2.88097e-01  0.3333030000  0.38387100  0.3338621327
#> w1[3]            2.59910e-01  0.3102180000  0.35836200  0.3096683633
#> w1[4]            6.16440e-06  0.0005610100  0.02765660  0.0039372013
#> w1[5]            7.14493e-10  0.0000548948  0.00367727  0.0007723635
#> w2[1]            5.33781e-11  0.0931336000  0.90453900  0.2391847626
#> w2[2]            3.01559e-11  0.0641686000  0.90519100  0.2243244062
#> w2[3]            4.13718e-12  0.0295937000  0.85736100  0.1791966123
#> w2[4]            1.46942e-15  0.2530120000  0.97442900  0.3572942244
#> w3[1]            2.66967e-01  0.3630540000  0.45435900  0.3634170657
#> w3[2]            1.82342e-01  0.2813570000  0.37528500  0.2807463649
#> w3[3]            2.06964e-01  0.2976065000  0.39850800  0.2988328809
#> w3[4]            1.64042e-05  0.0396070000  0.13203800  0.0486607099
#> w3[5]            1.00739e-10  0.0019389200  0.04176720  0.0083429784
#>                          SD Mode        MCerr MC%ofSD SSeff        AC.10
#> B0              0.660950751   NA 0.0081843238     1.2  6522 -0.005036043
#> B1              0.022611448   NA 0.0002668354     1.2  7181  0.005340841
#> B2              0.019753114   NA 0.0004696327     2.4  1769  0.116928700
#> B3              0.021771162   NA 0.0002634621     1.2  6829  0.019492507
#> phi_weight      0.009655255   NA 0.0001196751     1.2  6509 -0.007746521
#> phi_city_city_2 0.066072957   NA 0.0016066911     2.4  1691  0.143373168
#> phi_city_city_3 0.062526347   NA 0.0008860734     1.4  4980  0.034062612
#> w1[1]           0.025269015   NA 0.0006129529     2.4  1700  0.080146087
#> w1[2]           0.024590205   NA 0.0005995454     2.4  1682  0.131792578
#> w1[3]           0.025358605   NA 0.0006213296     2.5  1666  0.102504953
#> w1[4]           0.009197432   NA 0.0014137568    15.4    42  0.918815445
#> w1[5]           0.002155731   NA 0.0003151204    14.6    47  0.922797114
#> w2[1]           0.296331789   NA 0.0149296845     5.0   394  0.512573392
#> w2[2]           0.296851761   NA 0.0157227851     5.3   356  0.584325126
#> w2[3]           0.275878758   NA 0.0154643527     5.6   318  0.600834420
#> w2[4]           0.353545039   NA 0.0213802938     6.0   273  0.640013439
#> w3[1]           0.048455381   NA 0.0012350432     2.5  1539  0.129312309
#> w3[2]           0.049648639   NA 0.0013841516     2.8  1287  0.204552118
#> w3[3]           0.048614964   NA 0.0013054725     2.7  1387  0.171244393
#> w3[4]           0.042137794   NA 0.0033850988     8.0   155  0.771835837
#> w3[5]           0.014595877   NA 0.0018318112    12.6    63  0.850216402
#>                     psrf
#> B0              1.000263
#> B1              1.000481
#> B2              1.001490
#> B3              1.001598
#> phi_weight      1.000031
#> phi_city_city_2 1.001370
#> phi_city_city_3 1.000487
#> w1[1]           1.001298
#> w1[2]           1.011604
#> w1[3]           1.003032
#> w1[4]           1.235120
#> w1[5]           1.926119
#> w2[1]           1.102668
#> w2[2]           1.011382
#> w2[3]           1.054547
#> w2[4]           1.184856
#> w3[1]           1.003092
#> w3[2]           1.002040
#> w3[3]           1.003471
#> w3[4]           1.163290
#> w3[5]           1.346692
```

We can analyze the mixing of the markov chains and corresponding
autocorrelation plots using functions from the coda package. We will
output only the first few plots.

``` r
par(mfrow = c(3,3))
coda::traceplot(bgwqsr_fit$model$mcmc) # traceplot
```

<img src="man/figures/README-unnamed-chunk-19-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-19-2.png" width="100%" /><img src="man/figures/README-unnamed-chunk-19-3.png" width="100%" />

``` r
coda::autocorr.plot(bgwqsr_fit$model$mcmc, ask = F) # autocorrelation plot
```

<img src="man/figures/README-unnamed-chunk-20-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-20-2.png" width="100%" /><img src="man/figures/README-unnamed-chunk-20-3.png" width="100%" /><img src="man/figures/README-unnamed-chunk-20-4.png" width="100%" /><img src="man/figures/README-unnamed-chunk-20-5.png" width="100%" /><img src="man/figures/README-unnamed-chunk-20-6.png" width="100%" /><img src="man/figures/README-unnamed-chunk-20-7.png" width="100%" />

``` r
par(mfrow = c(3,3))
coda::densplot(bgwqsr_fit$model$mcmc) # posterior density plots
```

<img src="man/figures/README-unnamed-chunk-21-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-21-2.png" width="100%" /><img src="man/figures/README-unnamed-chunk-21-3.png" width="100%" />

``` r
plot(bgwqsr_fit$model$mcmc) # combines traceplot() and densplot()
```

<img src="man/figures/README-unnamed-chunk-22-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-22-2.png" width="100%" /><img src="man/figures/README-unnamed-chunk-22-3.png" width="100%" /><img src="man/figures/README-unnamed-chunk-22-4.png" width="100%" /><img src="man/figures/README-unnamed-chunk-22-5.png" width="100%" /><img src="man/figures/README-unnamed-chunk-22-6.png" width="100%" />

Finally, we can look at plots for the posterior credible intervals for
group indices, single chemical weights, and corresponding posterior
means. To do this, we can use the functions `plot_weights()`,
`plot_betas()`, and `plot_result()`. `plot_results()` combines both
plots generates by plot_weights and plot_betas into a side by side
figure. Below are code examples for each of the three function calls.

``` r
plot_betas(bgwqsr_fit)
```

<img src="man/figures/README-unnamed-chunk-23-1.png" width="100%" />

``` r
plot_weights(bgwqsr_fit)
```

<img src="man/figures/README-unnamed-chunk-24-1.png" width="100%" />

``` r
plot_results(bgwqsr_fit)
```

<img src="man/figures/README-unnamed-chunk-25-1.png" width="100%" />
