
# fgwqsr test

testthat::test_that("Expectation: fgwqsr succeeds in fitting and summary works", {
  n = 10000

  gamma = log(c(.5,1,1.5)) # group index sizes in log odds scale

  ess = list(n = n,
             weights = list(w1 = c(1/3,1/3,1/3,0,0),
                            w2 = c(1/2,1/2, 0, 0),
                            w3 = c(1/3,1/3,1/3,0,0)
             )
  )

  ccs = c(.5,.1)

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

  corr = create_corr_mat(ESS = ess, CCS = ccs)

  corr %>% data.frame

  set.seed(1) # for reproducibility
  chem_data = mvtnorm::rmvnorm(n, mean = rep(0, 14), sigma = corr)
  chem_data = apply(chem_data, MARGIN = 2, statar::xtile, n = 5) - 1
  set.seed(1)

  # create adjusting covariates
  weight = stats::rnorm(n = n, mean = 68, sd = 2.5)
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

  pi = plogis(logit_pi)

  # transform to bernoulli outcome
  y = sapply(pi, FUN = function(p) rbinom(1,1,p))

  # create dataset

  data = data.frame(y = y, chem_data, weight = weight, city = city)

  mod_formula_adj = y ~ X1 + X2 + X3 + X4 + X5 | X6 + X7 + X8 + X9 | X10 + X11 + X12 + X13 + X14 / weight + i.city

  fgwqsr_fit_adj = fgwqsr::fgwqsr(formula = mod_formula_adj,
                                  data = data,
                                  quantiles = 5,
                                  n_mvn_sims = 10000,
                                  verbose = F)

  testthat::expect_success(
    testthat::expect_type(fgwqsr_fit_adj , "list")
  )

  testthat::expect_no_error(summary(fgwqsr_fit_adj))

})


testthat::test_that("Expectation: bgwqsr succeeds in fitting and summary works", {
  n = 10000

  gamma = log(c(.5,1,1.5)) # group index sizes in log odds scale

  ess = list(n = n,
             weights = list(w1 = c(1/3,1/3,1/3,0,0),
                            w2 = c(1/2,1/2, 0, 0),
                            w3 = c(1/3,1/3,1/3,0,0)
             )
  )

  ccs = c(.5,.1)

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

  corr = create_corr_mat(ESS = ess, CCS = ccs)

  corr %>% data.frame

  set.seed(1) # for reproducibility
  chem_data = mvtnorm::rmvnorm(n, mean = rep(0, 14), sigma = corr)
  chem_data = apply(chem_data, MARGIN = 2, statar::xtile, n = 5) - 1
  set.seed(1)

  # create adjusting covariates
  weight = stats::rnorm(n = n, mean = 68, sd = 2.5)
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

  pi = plogis(logit_pi)

  # transform to bernoulli outcome
  y = sapply(pi, FUN = function(p) rbinom(1,1,p))

  # create dataset

  data = data.frame(y = y, chem_data, weight = weight, city = city)

  mod_formula_adj = y ~ X1 + X2 + X3 + X4 + X5 | X6 + X7 + X8 + X9 | X10 + X11 + X12 + X13 + X14 / weight + i.city

  bgwqsr_fit_adj = fgwqsr::bgwqsr(formula = mod_formula_adj,
                                  data = data,
                                  quantiles = 5,
                                  n.iter = 1000,
                                  n.burnin = 40,
                                  n.thin = 1, n.chains = 3,
                                  n.adapt = 10,
                                  method = "parallel")

  testthat::expect_success(
    testthat::expect_type(bgwqsr_fit_adj , "list")
  )

  testthat::expect_no_error(plot_results(bgwqsr_fit_adj))
  testthat::expect_no_error(plot_betas(bgwqsr_fit_adj))
  testthat::expect_no_error(plot_weights(bgwqsr_fit_adj))

})


