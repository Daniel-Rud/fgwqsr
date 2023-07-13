
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fgwqsr

<!-- badges: start -->
<!-- badges: end -->

# Introduction

<br> FGWQSR estimates parameters to the following model: where

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

and $\sum_{k = 1}^{c_g} w_{g,k} =1$, $w_{g,k} \in (0,1)$.

<br>

## Installation

You can install the development version of fgwqsr from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Daniel-Rud/fgwqsr_package")
```

``` r
library(fgwqsr)
## basic example code
```

# Load Packages

First, we load the FGWQSR package and other libraries we will be using
for the excerise.

``` r
library(fgwqsr)
```

# Data Generation

<br> To begin, let us generate some data to use. We will have 14
underlying chemical exposure variables, where we will have mixture
groups of sizes (5,4,5). We will use a fixed correlation structure such
that the correlation of chemicals within a group is .7 and the cross
correlations between 2 chemicals in different groups is .5. <br> We
first define ESS and CCS. ESS is a list that contains the the sample
size and the true underlying weight distribution of the chemicals in
each group. CCS is a two element vector where the first index defines
the within group correlation and the second index controls the between
group correlation. The sample size for our dataset is fixed to
$n = 10,000$ and we set the distribution of the true underlying weights
as follows: *$\mathbf{w}_1 = (1/3,1/3,1/3,0,0)$
*$\mathbf{w}_2 = (1/2,1/2, 0, 0)$
\*$\mathbf{w}_3 = (1/3,1/3,1/3,0,0)$.  
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
formulas: <br> \* $|$ - denotes the boundary of a mixture group, used to
seperate chemicals within a mixture group \* $/$ - denotes the end of
the mixture group specification, adjusting covariates can be added to
the formula after this character. If no adjusting covariates, do not
need to specify \* i. - precedes categorical variables to denote a
categorical variable. For example, if we have a categorical variable
cat_var, we would denote this in the model formula by i.cat_var. This is
similar to the stata syntax to declare categorical variables.  
<br> The fgwqsr() function has other options too: <br> \* data - should
be a data frame object containing variable columnames referenced in
model formula (cat vars do not need to be named with i. in columnames)
\* quantiles - number of quantiles to quantize the exposure variables in
the mixture portion of the model \* n_mvn_sims - defines resolution for
simulated null distribution for group index and single chemical LRTs.
Default is 10,000.  
\* zero_threshold_cutoff - defines a . \* verbose - Displays messages
and progress bar while fitting FGWQSR model. Default is TRUE.  
\* cores - number of cores to parallelize on for fitting nested models
and simulated null LRT distributions. Default is number of available
cores on user device.  
\* optim_control_list - option to supply control options to optim.

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
#> Fitting nested models for Weight Inference:
#> 
#> Now Performing Inference...
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
#>                  Estimate        LRT P-value    
#> Mixture Index 1 -0.702387 1168.18784  <1e-06 ***
#> Mixture Index 2  0.024676    1.49567    0.57    
#> Mixture Index 3  0.411345  408.88836  <1e-06 ***
#> 
#> Estimates and Inference for Weights
#>    Weight Estimate        LRT P-value    
#> X1        0.305501 128.107913  <1e-06 ***
#> X2        0.336426 156.350666  <1e-06 ***
#> X3        0.322359 142.164908  <1e-06 ***
#> X4        0.023634   0.750190    0.21    
#> X5        0.012080   0.202268    0.35    
#> -------------------------------------------------
#>    Weight Estimate      LRT P-value  
#> X6        0.248205 0.102211    0.42  
#> X7        0.207174 0.071973    0.29  
#> X8        0.111892 0.021212    0.47  
#> X9        0.432729 0.309532    0.26  
#> -------------------------------------------------
#>     Weight Estimate       LRT P-value    
#> X10        0.245935 26.931807  <1e-06 ***
#> X11        0.310206 43.508110  <1e-06 ***
#> X12        0.376066 62.982156  <1e-06 ***
#> X13        0.001262  0.000717    0.45    
#> X14        0.066531  2.051216    0.09   .
#> -------------------------------------------------
#> 
#> Estimates and Inference for Intercept and Adjusting Covariates
#>            Estimate       SE         Z P(Z > |z|)                95% CI  
#> intercept -0.012925 0.064598 -0.200085   0.841414 (-0.139536, 0.113685)  
#> 
#> Significance Codes: <0.001 '***' <0.01 '**' <0.05 '*' <0.10 '.' 
#> 
#> Total runtime for FGWQSR:  1.36 minutes on 10 cores.
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
#>                 True Group Index Estimated Group Index Signficiant?
#> Mixture Index 1              0.5                 0.495          Yes
#> Mixture Index 2              1.0                 1.025           No
#> Mixture Index 3              1.5                 1.509          Yes
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
#> Fitting nested models for Weight Inference:
#> 
#> Now Performing Inference...
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
#>                  Estimate       LRT P-value    
#> Mixture Index 1 -0.721974 1153.2447  <1e-06 ***
#> Mixture Index 2  0.070848   13.1311    0.01  **
#> Mixture Index 3  0.403725  375.4276  <1e-06 ***
#> 
#> Estimates and Inference for Weights
#>    Weight Estimate        LRT P-value    
#> X1        0.367200 178.322269  <1e-06 ***
#> X2        0.287866 110.033798  <1e-06 ***
#> X3        0.337621 150.472596  <1e-06 ***
#> X4        0.003173   0.013089    0.54    
#> X5        0.004141   0.026154    0.43    
#> -------------------------------------------------
#>    Weight Estimate      LRT P-value  
#> X6        0.027973 0.009814    0.53  
#> X7        0.652795 5.397737    0.02 *
#> X8        0.147947 0.280936    0.33  
#> X9        0.171285 0.367393    0.24  
#> -------------------------------------------------
#>     Weight Estimate     LRT P-value    
#> X10        0.318641 42.0794  <1e-06 ***
#> X11        0.330624 45.7642  <1e-06 ***
#> X12        0.350735 50.9094  <1e-06 ***
#> X13        0.000000  0.0000    0.74    
#> X14        0.000000  0.0000    0.78    
#> -------------------------------------------------
#> 
#> Estimates and Inference for Intercept and Adjusting Covariates
#>              Estimate       SE          Z P(Z > |z|)                 95% CI    
#> intercept   -2.661611 0.624433  -4.262446      2e-05 (-3.885477, -1.437746) ***
#> weight       0.039620 0.009127   4.340902    1.4e-05   (0.021731, 0.057509) ***
#> city_city_2  0.216813 0.054042   4.011911      6e-05   (0.110892, 0.322733) ***
#> city_city_3 -1.035485 0.059087 -17.524742     <1e-06 (-1.151293, -0.919676) ***
#> 
#> Significance Codes: <0.001 '***' <0.01 '**' <0.05 '*' <0.10 '.' 
#> 
#> Total runtime for FGWQSR:  1.3 minutes on 10 cores.
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
#>                 True Group Index Estimated Group Index Signficiant?
#> Mixture Index 1              0.5                 0.486          Yes
#> Mixture Index 2              1.0                 1.073          Yes
#> Mixture Index 3              1.5                 1.497          Yes
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

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.
