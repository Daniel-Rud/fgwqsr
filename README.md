
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fgwqsr – Frequentist Grouped Weighted Quantile Sum Regression <img src="man/figures/FGWQSR_hex.png" width="100" align="right" />

Fully frequentist estimation and inference for Weighted Quantile Sum
Regression, enabling inference on both individual and grouped pollutants
without the need for data splitting. [Rud et al
(2025)](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.70078).

Also included in this package is an implementation of Bayesian Grouped
Weighted Quantile Sum Regression (bgwqsr) introduced by Wheeler, David C
et al. that fits Markov Chain Monte Carlo (MCMC) chains in parallel
through the runjags package.

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
[GitHub](https://github.com/Daniel-Rud/fgwqsr) with:

``` r
# install.packages("devtools")
remotes::install_github("Daniel-Rud/fgwqsr")
#> Using github PAT from envvar GITHUB_PAT. Use `gitcreds::gitcreds_set()` and unset GITHUB_PAT in .Renviron (or elsewhere) if you want to use the more secure git credential store instead.
#> Downloading GitHub repo Daniel-Rud/fgwqsr@HEAD
#> 
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>      checking for file ‘/private/var/folders/b4/9kg7p6cj729_pzc5dggk_9jm0000gn/T/RtmpelvZTp/remotes121c86043a123/Daniel-Rud-fgwqsr-8814766/DESCRIPTION’ ...  ✔  checking for file ‘/private/var/folders/b4/9kg7p6cj729_pzc5dggk_9jm0000gn/T/RtmpelvZTp/remotes121c86043a123/Daniel-Rud-fgwqsr-8814766/DESCRIPTION’
#>   ─  preparing ‘fgwqsr’:
#>      checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
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

The `fgwqsr` package depends on JAGS installation through the `runjags`
package, which is used to fit BGWQSR models. Please see the `runjags`
installation instructions for your operating system
[here](https://cran.r-project.org/web/packages/runjags/vignettes/Installation.html).
Note that JAGS is not required to fit FGWQSR models.

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
that we have an interpretable model intercept. For ease of
understanding, we give the chemical variables arbitrary names of
chemical constituents.

``` r
chem_data = apply(chem_data, MARGIN = 2, statar::xtile, n = 5) - 1
colnames(chem_data) = chemicals <- c(
  "PM2.5", 
  "NO2", 
  "O3", 
  "SO2", 
  "CO", 
  "Lead", 
  "Mercury", 
  "Arsenic", 
  "Cadmium", 
  "Benzene", 
  "Toluene", 
  "Phthalates", 
  "Bisphenol_A", 
  "Polychlorinated_Biphenyls"
)

head(chem_data %>% data.frame)
#>   PM2.5 NO2 O3 SO2 CO Lead Mercury Arsenic Cadmium Benzene Toluene Phthalates
#> 1     0   1  0   0  3    1       4       3       2       0       0          0
#> 2     0   1  1   3  1    0       0       1       1       2       2          2
#> 3     1   0  0   0  2    0       0       2       1       3       1          1
#> 4     3   0  3   1  0    4       4       3       4       2       1          3
#> 5     2   3  4   3  1    4       2       0       3       4       1          1
#> 6     3   2  1   3  4    2       1       2       0       2       2          2
#>   Bisphenol_A Polychlorinated_Biphenyls
#> 1           0                         0
#> 2           1                         1
#> 3           0                         1
#> 4           2                         3
#> 5           1                         3
#> 6           1                         0
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
#>   y PM2.5 NO2 O3 SO2 CO Lead Mercury Arsenic Cadmium Benzene Toluene Phthalates
#> 1 0     0   1  0   0  3    1       4       3       2       0       0          0
#> 2 1     0   1  1   3  1    0       0       1       1       2       2          2
#> 3 0     1   0  0   0  2    0       0       2       1       3       1          1
#> 4 0     3   0  3   1  0    4       4       3       4       2       1          3
#> 5 0     2   3  4   3  1    4       2       0       3       4       1          1
#> 6 1     3   2  1   3  4    2       1       2       0       2       2          2
#>   Bisphenol_A Polychlorinated_Biphenyls
#> 1           0                         0
#> 2           1                         1
#> 3           0                         1
#> 4           2                         3
#> 5           1                         3
#> 6           1                         0
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

mod_formula = y ~ PM2.5 + NO2 + O3 + SO2 + CO | Lead + Mercury + Arsenic + Cadmium | Benzene + Toluene + Phthalates + Bisphenol_A + Polychlorinated_Biphenyls
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
#> Performing Likelihood Ratio Test Inference...
```

<br>

We can see the model summary using the call `summary()`

``` r

summary(fgwqsr_fit)
#> 
#> Call: 
#> FGWQSR with formula 'y ~ PM2.5 + NO2 + O3 + SO2 + CO | Lead + Mercury + Arsenic + Cadmium | Benzene + Toluene + Phthalates + Bisphenol_A + Polychlorinated_Biphenyls' on n = 10000 observations and family = 'binomial'.
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
#>       Weight Estimate       LRT P-value    
#> PM2.5       0.2999211 112.43088  <2e-16 ***
#> NO2         0.3417464 143.95821  <2e-16 ***
#> O3          0.3225850 128.82142  <2e-16 ***
#> SO2         0.0357475   1.58066  0.0867 .  
#> CO          0.0000000   0.00000  1.0000    
#> -------------------------------------------------
#>         Weight Estimate     LRT P-value  
#> Lead                  0 0.00000  1.0000  
#> Mercury               0 0.00000  1.0000  
#> Arsenic               0 0.00000  1.0000  
#> Cadmium               1 1.45345  0.0828 .
#> -------------------------------------------------
#>                           Weight Estimate     LRT P-value    
#> Benzene                          0.339334 54.3167  <2e-16 ***
#> Toluene                          0.347920 57.0600  <2e-16 ***
#> Phthalates                       0.312745 46.7864  <2e-16 ***
#> Bisphenol_A                      0.000000  0.0000       1    
#> Polychlorinated_Biphenyls        0.000000  0.0000       1    
#> -------------------------------------------------
#> 
#> Estimates and Inference for Intercept and Adjusting Covariates
#> 
#>             Estimate      SE        Z P(Z > |z|)
#> intercept -0.0866543 0.06543 -1.32441    0.18537
#> 
#> Significance Codes: <0.001 '***' <0.01 '**' <0.05 '*' <0.10 '.' 
#> 
#> Total runtime for FGWQSR:  1.42 minutes on 10 cores.
```

<br>

We can plot a visual summary of the model output using `plot()`

``` r

plot(fgwqsr_fit)
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />
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
comparison between San Bernardino and reference Los Angeles, and -1 to
the comparison between Santa Barbara and Los Angeles.

<br>

``` r
set.seed(11)
# create adjusting covariates
weight = rnorm(n = n, mean = 68, sd = 2.5)
city = sample(c("Los_Angeles", "San_Bernardino", "Santa_Barbara"), size = n, replace = T)

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
#>   y PM2.5 NO2 O3 SO2 CO Lead Mercury Arsenic Cadmium Benzene Toluene Phthalates
#> 1 1     0   1  0   0  3    1       4       3       2       0       0          0
#> 2 1     0   1  1   3  1    0       0       1       1       2       2          2
#> 3 0     1   0  0   0  2    0       0       2       1       3       1          1
#> 4 1     3   0  3   1  0    4       4       3       4       2       1          3
#> 5 1     2   3  4   3  1    4       2       0       3       4       1          1
#> 6 1     3   2  1   3  4    2       1       2       0       2       2          2
#>   Bisphenol_A Polychlorinated_Biphenyls   weight           city
#> 1           0                         0 66.52242    Los_Angeles
#> 2           1                         1 68.06649    Los_Angeles
#> 3           0                         1 64.20862  Santa_Barbara
#> 4           2                         3 64.59337 San_Bernardino
#> 5           1                         3 70.94622  Santa_Barbara
#> 6           1                         0 65.66462    Los_Angeles
```

<br> Now, we specify the formula to include the continuous and
categorical variables and call `fgwqsr()`.

``` r

mod_formula_adj = y ~ PM2.5 + NO2 + O3 + SO2 + CO | Lead + Mercury + Arsenic + Cadmium | Benzene + Toluene + Phthalates + Bisphenol_A + Polychlorinated_Biphenyls / weight + i.city

fgwqsr_fit_adj = fgwqsr(formula = mod_formula_adj,
                        data = data,
                        quantiles = 5,
                        family = 'binomial',
                        n_mvn_sims = 10000,
                        verbose = T)
#> Fitting full and nested FGWQSR models...
#> 
#> Performing Likelihood Ratio Test Inference...
```

<br> Again, we use the function `summary()` to view the results of the
`fgwqsr()` call.

``` r

summary(fgwqsr_fit_adj)
#> 
#> Call: 
#> FGWQSR with formula 'y ~ PM2.5 + NO2 + O3 + SO2 + CO | Lead + Mercury + Arsenic + Cadmium | Benzene + Toluene + Phthalates + Bisphenol_A + Polychlorinated_Biphenyls/weight + i.city' on n = 10000 observations and family = 'binomial'.
#> 
#> 10000 samples used for simulated LRT distirbution.
#> 
#> Log Likelihood: -5981.247 | AIC: 11998.49 | BIC: 12128.28
#> 
#> Estimates and Inference for Group Index Effects
#> 
#>                    Estimate        LRT P-value    
#> Mixture Effect 1 -0.7118585 1249.38144  <2e-16 ***
#> Mixture Effect 2 -0.0236949    1.74268  0.5447    
#> Mixture Effect 3  0.4216609  413.91554  <2e-16 ***
#> 
#> Estimates and Inference for Weights
#> 
#>       Weight Estimate     LRT P-value    
#> PM2.5        0.354744 187.432  <2e-16 ***
#> NO2          0.358079 188.728  <2e-16 ***
#> O3           0.287177 122.582  <2e-16 ***
#> SO2          0.000000   0.000       1    
#> CO           0.000000   0.000       1    
#> -------------------------------------------------
#>         Weight Estimate     LRT P-value
#> Lead           0.000000 0.00000  1.0000
#> Mercury        0.700016 0.89830  0.1247
#> Arsenic        0.299984 0.16411  0.2555
#> Cadmium        0.000000 0.00000  1.0000
#> -------------------------------------------------
#>                           Weight Estimate      LRT P-value    
#> Benzene                          0.319127 49.15113  <2e-16 ***
#> Toluene                          0.295642 42.30203  <2e-16 ***
#> Phthalates                       0.260109 33.11722  <2e-16 ***
#> Bisphenol_A                      0.125123  7.82728  0.0022 ** 
#> Polychlorinated_Biphenyls        0.000000  0.00000  1.0000    
#> -------------------------------------------------
#> 
#> Estimates and Inference for Intercept and Adjusting Covariates
#> 
#>                        Estimate      SE         Z P(Z > |z|)    
#> intercept           -3.48925597 0.60888 -5.730591 1.0008e-08 ***
#> weight               0.05377299 0.00889  6.048480 1.4622e-09 ***
#> city_San_Bernardino  0.01105269 0.05374  0.205688    0.83703    
#> city_Santa_Barbara  -0.01047868 0.05414 -0.193551    0.84653    
#> 
#> Significance Codes: <0.001 '***' <0.01 '**' <0.05 '*' <0.10 '.' 
#> 
#> Total runtime for FGWQSR:  1.57 minutes on 10 cores.
```

<br>

<br>

Plotting a visual summary of the model output using `plot()` …

``` r

plot(fgwqsr_fit)
```

<img src="man/figures/README-unnamed-chunk-17-1.png" width="100%" />

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
#> Mixture Effect 1              0.5                 0.491          Yes
#> Mixture Effect 2              1.0                 0.977           No
#> Mixture Effect 3              1.5                 1.524          Yes
#>     True Weight Estimated Weight Signficiant?
#> w11       0.333            0.355          Yes
#> w12       0.333            0.358          Yes
#> w13       0.333            0.287          Yes
#> w14           0            0.000           No
#> w15           0            0.000           No
#> w21           -            0.000           No
#> w22           -            0.700           No
#> w23           -            0.300           No
#> w24           -            0.000           No
#> w31       0.333            0.319          Yes
#> w32       0.333            0.296          Yes
#> w33       0.333            0.260          Yes
#> w34           0            0.125          Yes
#> w35           0            0.000           No
#>                     True Covariate Effect Estimated Covariate Effect
#> intercept                            -3.3                -3.48925597
#> weight                                0.5                 0.05377299
#> city_San_Bernardino                   0.2                 0.01105269
#> city_Santa_Barbara                   -1.0                -0.01047868
#>                     Significant?
#> intercept                    Yes
#> weight                       Yes
#> city_San_Bernardino           No
#> city_Santa_Barbara            No
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

A note of **caution**: JAGS does not allow for variables with spaces in
their names – please name variables appropriately!

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
#> Welcome to JAGS 4.3.2 (official binary) on Sun Aug 31 12:13:59 2025
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
#>                              Lower95       Median    Upper95         Mean
#> B0                      -4.51673e+00 -3.234565000 -2.0556200 -3.240077827
#> B1                      -7.55912e-01 -0.712505500 -0.6699580 -0.712707515
#> B2                      -4.60343e-02 -0.006361070  0.0308357 -0.006839838
#> B3                       3.71839e-01  0.414750500  0.4566030  0.414595520
#> phi_weight               3.18810e-02  0.049787700  0.0679265  0.049852758
#> phi_city_San_Bernardino -4.42564e-02  0.001908540  0.0482238  0.002743681
#> phi_city_Santa_Barbara  -4.99132e-02 -0.002199065  0.0414722 -0.003038943
#> w1[1]                    3.04301e-01  0.353853000  0.3984010  0.353346952
#> w1[2]                    3.07197e-01  0.353501000  0.4004670  0.353768675
#> w1[3]                    2.38390e-01  0.284111000  0.3318060  0.284609213
#> w1[4]                    8.79485e-07  0.002227950  0.0242490  0.005940134
#> w1[5]                    3.09730e-08  0.000229071  0.0115766  0.002335033
#> w2[1]                    1.80818e-09  0.073018950  0.8729500  0.219533553
#> w2[2]                    2.53242e-10  0.124019000  0.9152900  0.283659605
#> w2[3]                    8.74800e-12  0.142057000  0.9123200  0.276559296
#> w2[4]                    2.88558e-08  0.079288600  0.8521550  0.220247550
#> w3[1]                    2.33887e-01  0.325652000  0.4128360  0.325486426
#> w3[2]                    2.16828e-01  0.300339000  0.3932940  0.301239531
#> w3[3]                    1.74056e-01  0.260643000  0.3567910  0.261618549
#> w3[4]                    9.07646e-04  0.109226500  0.1878680  0.106964416
#> w3[5]                    5.19473e-09  0.000166607  0.0245559  0.004691055
#>                                  SD Mode        MCerr MC%ofSD SSeff
#> B0                      0.626959500   NA 0.0071901968     1.1  7603
#> B1                      0.021856942   NA 0.0002490331     1.1  7703
#> B2                      0.019647701   NA 0.0005323049     2.7  1362
#> B3                      0.021767741   NA 0.0003117522     1.4  4875
#> phi_weight              0.009157264   NA 0.0001067935     1.2  7353
#> phi_city_San_Bernardino 0.022078554   NA 0.0002027716     0.9 11856
#> phi_city_Santa_Barbara  0.022033935   NA 0.0001996468     0.9 12180
#> w1[1]                   0.023802653   NA 0.0005695187     2.4  1747
#> w1[2]                   0.023936764   NA 0.0005643569     2.4  1799
#> w1[3]                   0.024260656   NA 0.0005943187     2.4  1666
#> w1[4]                   0.008534823   NA 0.0012082987    14.2    50
#> w1[5]                   0.004925966   NA 0.0004853203     9.9   103
#> w2[1]                   0.283558992   NA 0.0130010415     4.6   476
#> w2[2]                   0.319709746   NA 0.0153566404     4.8   433
#> w2[3]                   0.307104507   NA 0.0144209047     4.7   454
#> w2[4]                   0.280572692   NA 0.0124012762     4.4   512
#> w3[1]                   0.046271693   NA 0.0012691338     2.7  1329
#> w3[2]                   0.045269872   NA 0.0012065383     2.7  1408
#> w3[3]                   0.046845087   NA 0.0013317882     2.8  1237
#> w3[4]                   0.050891020   NA 0.0028554960     5.6   318
#> w3[5]                   0.012826129   NA 0.0016683974    13.0    59
#>                                AC.10     psrf
#> B0                      -0.019209941 1.000450
#> B1                       0.011448491 1.005023
#> B2                       0.156969017 1.001026
#> B3                       0.032430960 1.003198
#> phi_weight              -0.019678517 1.000449
#> phi_city_San_Bernardino -0.007921566 1.000239
#> phi_city_Santa_Barbara  -0.019274547 1.000421
#> w1[1]                    0.106027568 1.010543
#> w1[2]                    0.083389597 1.005137
#> w1[3]                    0.117619789 1.006981
#> w1[4]                    0.911847123 1.492687
#> w1[5]                    0.927279147 3.717224
#> w2[1]                    0.515293919 1.031223
#> w2[2]                    0.557117046 1.063476
#> w2[3]                    0.522599017 1.060445
#> w2[4]                    0.499738604 1.009667
#> w3[1]                    0.199308931 1.004396
#> w3[2]                    0.195566134 1.004504
#> w3[3]                    0.250774293 1.004680
#> w3[4]                    0.639730871 1.032039
#> w3[5]                    0.853183850 1.234071
```

We can analyze the mixing of the markov chains and corresponding
autocorrelation plots using functions from the coda package. We will
output only the first few plots.

``` r
par(mfrow = c(3,3))
coda::traceplot(bgwqsr_fit$model$mcmc) # traceplot
```

<img src="man/figures/README-unnamed-chunk-21-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-21-2.png" width="100%" /><img src="man/figures/README-unnamed-chunk-21-3.png" width="100%" />

``` r
coda::autocorr.plot(bgwqsr_fit$model$mcmc, ask = F) # autocorrelation plot
```

<img src="man/figures/README-unnamed-chunk-22-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-22-2.png" width="100%" /><img src="man/figures/README-unnamed-chunk-22-3.png" width="100%" /><img src="man/figures/README-unnamed-chunk-22-4.png" width="100%" /><img src="man/figures/README-unnamed-chunk-22-5.png" width="100%" /><img src="man/figures/README-unnamed-chunk-22-6.png" width="100%" /><img src="man/figures/README-unnamed-chunk-22-7.png" width="100%" />

``` r
par(mfrow = c(3,3))
coda::densplot(bgwqsr_fit$model$mcmc) # posterior density plots
```

<img src="man/figures/README-unnamed-chunk-23-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-23-2.png" width="100%" /><img src="man/figures/README-unnamed-chunk-23-3.png" width="100%" />

``` r
plot(bgwqsr_fit$model$mcmc) # combines traceplot() and densplot()
```

<img src="man/figures/README-unnamed-chunk-24-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-24-2.png" width="100%" /><img src="man/figures/README-unnamed-chunk-24-3.png" width="100%" /><img src="man/figures/README-unnamed-chunk-24-4.png" width="100%" /><img src="man/figures/README-unnamed-chunk-24-5.png" width="100%" /><img src="man/figures/README-unnamed-chunk-24-6.png" width="100%" />

Finally, we can look at plots for the posterior credible intervals for
group indices, single chemical weights, and corresponding posterior
means. To do this, we can use the functions `plot_weights()`,
`plot_betas()`, and `plot_result()`. `plot_results()` combines both
plots generates by plot_weights and plot_betas into a side by side
figure. Below are code examples for each of the three function calls.

``` r
plot_betas(bgwqsr_fit)
```

<img src="man/figures/README-unnamed-chunk-25-1.png" width="100%" />

``` r
plot_weights(bgwqsr_fit)
```

<img src="man/figures/README-unnamed-chunk-26-1.png" width="100%" />

``` r
plot_results(bgwqsr_fit)
```

<img src="man/figures/README-unnamed-chunk-27-1.png" width="100%" />

# Community Guidelines

We welcome contributions and engagement from the community to help
improve and extend the functionality of `fgwqsr`. To maintain a
collaborative and supportive environment, please follow these
guidelines:

## Contributing

- Contributions of all kinds are welcome, including code improvements,
  documentation updates, bug fixes, and new feature proposals.  
- Before submitting a pull request, please open an issue to discuss your
  idea or bug report.  
- Ensure that code contributions include appropriate documentation and
  tests where applicable.  
- Follow the existing code style and structure to maintain consistency.

## Reporting Issues

- If you encounter a problem, please [open an
  issue](https://github.com/Daniel-Rud/fgwqsr/issues) with a clear
  description of the bug or unexpected behavior.  
- Include reproducible examples (minimal code and sample data) whenever
  possible.  
- Clearly state your R version, operating system, and any other relevant
  session information (`sessionInfo()`).

### Seeking Support

- For usage questions or clarification about methods, please first
  consult the package vignettes and documentation.  
- If further support is needed, feel free to open a discussion or issue.
