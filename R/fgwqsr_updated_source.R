### Libraries ################################################################
library(fastDummies)
library(statar)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(rlang)
library(progressr)
library(future)
library(future.apply)
library(pracma)
library(sigmoid) # for softplus function
##############################################################################


clean_vars = function(vars)
{
  # separates formula into mixture component, continuous confounder,
  # and categorical confounder names
  sep = strsplit(vars, split = "/", fixed = TRUE)[[1]] # separate exposure vars from adjusting covariates

  mixture_string = sep[1] # exposure grouped vars

  # seperate into groups of exposures
  mixture_names = strsplit(strsplit(gsub(" ", "", mixture_string), "\\|")[[1]], "\\+")

  # find confounders/adjusting variables
  confounders = strsplit(gsub("\\+", " ", sep[2]), " ")[[1]]

  confounders = confounders[which(confounders!="")]

  # find indices of categorical variables
  cat_ind = which(substr(confounders, 1, 2) == "i.")

  if(length(cat_ind)>0) # if there are categorical variables
  {
    cat_confounders = gsub("i.", "", confounders[cat_ind], fixed = TRUE) # variable name with out i.

    cont_confounders = confounders[-cat_ind] # continuous confounders are left over variables

    variables = list(mixture = mixture_names, continuous = cont_confounders, categorical = cat_confounders)

  }else
  {
    if(length(confounders) > 0) # if there are confounders, but no categorical vars
    {

      cont_confounders = confounders

      variables = list(mixture = mixture_names, continuous = cont_confounders, categorical = NULL)

    }else # no adjusting variables
    {
      variables = list(mixture = mixture_names, continuous = NULL, categorical = NULL)
    }

  }

  return(variables)
}

quantize_vars = function(data, vars, quantiles)
{
  # quantize exposure grouped variables into quantile variables
  mixture_comps = unlist(vars$mixture)

  for(i in 1: length(mixture_comps))
  {
    data[[mixture_comps[i]]] = statar::xtile(data[[mixture_comps[i]]], n = quantiles) - 1
  }

  return(data)
}




get_formulas= function(vars,formula)
{
  # get formulas for LRTs
  weight_vars = unlist(vars$mixture)

  RHS = paste(" ", as.character(formula)[3], " ", sep = "") # this is for the edge cases of vector

  backslash_loc = gregexpr("/",RHS)[[1]] %>% as.numeric # need to add space where last covariate is before confounders

  if(backslash_loc != -1) # if there is a backslash (confounders present)
  {
    RHS = sub("/", " /", RHS)  # add space where backslash is
    backslash_loc = backslash_loc + 1
  }

  plus_locs = gregexpr("+", RHS, fixed = TRUE)[[1]]

  l = length(plus_locs)

  y = as.character(formula)[2] # string for dependent variable

  formulas = list()

  # generate single pollutant exluded formulas

  for(i in 1: length(weight_vars ))
  {
    loc_var = gregexpr(paste(" ", weight_vars[i]," ", sep = ""), RHS, fixed = TRUE)[[1]] # find loc of variable

    new_formula = sub(weight_vars[i], "", RHS, fixed = TRUE) # remove from formula

    which_plus1 = abs(plus_locs - rep(loc_var, l)) # find location dist of + from first index of var

    min_plus1 = min(which_plus1)

    which_plus1 = which(which_plus1 == min_plus1) # find which + is closest to 1st index

    which_plus2 = abs(plus_locs - rep((loc_var + nchar(weight_vars[i])-1), l)) # find loc dist of + from end of var

    min_plus2 = min(which_plus2)

    which_plus2 = which(which_plus2 == min(which_plus2))

    plus_remove = ifelse(min_plus1 < min_plus2,which_plus1, which_plus2) # remove + that is closest to word

    RHS_tokenize = strsplit(new_formula, "")[[1]]

    new_plus_locs = which(RHS_tokenize == "+")

    RHS_tokenize[new_plus_locs[plus_remove]] = ""

    new_formula = paste(y," ~ ",paste(RHS_tokenize, collapse = ""), sep = "")

    formulas[[i]] = new_formula
  }

  # get formulas for group lrt IF more than one group -- special case if one group

  num_groups = vars$mixture %>% length

  group_formulas = list()

  if(num_groups > 1) # if more than one group, perform the excluded group fg runs
  {

    bar_locs = gregexpr("|", RHS, fixed = T)[[1]] %>% as.numeric # where bars are between groups

    if(backslash_loc != -1)
    {
      bar_locs = c(bar_locs, backslash_loc)
    }else
    {
      bar_locs =  c(bar_locs, nchar(RHS) + 1)
    }

    # treat first group differently

    group_formulas[[1]] = paste(y, " ~ ", substring(RHS,bar_locs[1]+1), sep = "")

    for(i in 2:num_groups) # iterate over groups 2,3,...
    {

      group_formulas[[i]] = paste(y, " ~ ",substr(RHS, 1,bar_locs[i-1]-1),
                                  substr(RHS,bar_locs[i], nchar(RHS)), sep = "")

    }

  }else # get formula if we remove the only group in model
    # 1 group null ll will call glm %>% logLik
  {
    # find substring to remove
    end_loc = ifelse(backslash_loc != -1,backslash_loc+1, nchar(RHS))
    new_formula = substring(RHS, first = end_loc)

    # need to account for when there are categorical predictors

    if(vars$categorical %>% length > 0 )
    {
      new_cat_var_name = paste("factor(", vars$categorical, ")", sep = "")

      for(i in 1:length(vars$categorical)) # switch i. encoding to factor() encoding for glm
      {
        RHS = gsub(paste("i.", vars$categorical[i], sep = ""),
                   new_cat_var_name[i], RHS)
      }
    }

    group_formulas[[1]] = paste(y, "~",
                                ifelse(new_formula == "", "1", new_formula))
  }

  formulas = append(formulas, group_formulas)
  return(unlist(formulas))
}


create_likelihood_string = function(data, vars)
{
  ### Beta formulation##########################################################

  beta_string = "B[1]"

  beta_index = 2 # to hold value where group index effect is in vector

  num_mixes = length(vars$mixture) # number of mixtures

  beta_orders_mix = character(1+length(unlist(vars$mixture)))

  beta_orders_mix[1] = "constant"

  for(i in 1: length(vars$mixture)) # iterate over number of mixtures
  {
    num_elements = length(vars$mixture[[i]])

    group_string = ""

    if(num_elements>1) # if more than 1 element in mixture group
    {
      coef_denoms = paste("(1 + ", paste("exp(", paste("B[", (beta_index+1):(beta_index+num_elements -1),"]",
                                                       sep = ""), ")", sep = "", collapse = " + "),")", sep = "") # denominators of coefs

      coef_nums = c(paste("exp(B[", (beta_index+1):(beta_index+num_elements - 1),"])",
                          sep = ""), "1") # numerators of coefficients

      coeffs = paste("(",coef_nums, "/", coef_denoms, ")", sep = "") # vector of weight coefficient strings

      quant_vars = paste("data$", vars$mixture[[i]], "[i]", sep = "")

      group_string = paste("B[", beta_index, "]*(", paste(coeffs, quant_vars, sep = "*", collapse = " + "),
                           ")", sep = "")

      beta_string = paste(beta_string, group_string, sep = " + ")

      beta_orders_mix[beta_index:(beta_index+num_elements -1)] = c(
        paste("B_mix_", i, sep = ""), paste("alpha_", vars$mixture[[i]][-length(vars$mixture[[i]])], sep = ""))

    }else # if one element in mixture group, weight is 1
    {
      quant_vars = paste("data$", vars$mixture[[i]], "[i]", sep = "")

      group_string = paste("B[", beta_index, "]*",quant_vars, sep = "")

      beta_string = paste(beta_string, group_string, sep = " + ")

      beta_orders_mix[beta_index:(beta_index+num_elements -1)] = paste("B_mix_", i, sep = "")
    }

    beta_index = beta_index + length(vars$mixture[[i]])
  }

  ### For Phis #################################################################

  phi_start_index = beta_index

  length1 = length(vars$continuous)

  length2 = length(vars$categorical)

  phi1 = c()

  if(length1>0)
  {
    for(i in 1: length1)
    {
      phi1[i] = paste(vars$continuous[i], sep = "")
    }
  }

  phi2 = list()

  if(length2>0)
  {
    for(i in 1:length2)
    {
      varName = vars$categorical[i]

      names = names(data)[which(substr(names(data), 1, nchar(varName)+1) ==
                                  paste(varName, "_", sep = ""))]
      phi2[[i]] = paste(names, sep= "")

    }
  }

  phis = c(phi1, unlist(phi2))

  phi_betas = c()

  if(length(phis)>0) # if confounders, add betas
  {
    phi_betas = paste("B[", (phi_start_index):(phi_start_index + length(phis) -1), "]", sep = "")
  }


  order_of_betas = c(beta_orders_mix, phis)



  phiString = c()

  if(length(phis)>0)
    phiString = paste(" + ", paste(phi_betas, "*", "data$", phis,"[i]", sep = "", collapse = " + "))

  mixString = paste(beta_string, phiString, sep = "")

  return(list(mixString=mixString, order_of_betas=order_of_betas, numPhis = length(phis)))
}

make_beta_vec = function(B, vars)
{
  # we organize into a list, then perform the reparameterization

  mix_sizes = sapply(vars$mixture, length) # sizes of mixture groups
  num_mixes = vars$mixture %>% length # number of mixtures

  current_index = 2 # 1st index is constant

  for(i in 1:num_mixes) # iterate of number of mixtures
  {
    # if mix_size is 1, just leave param for unconstrained optimization
    if(mix_sizes[i] != 1)
    {
      # where is the group effect
      group_effect = B[current_index]
      # where are the alphas -- the 0 is so that exp(0) = 1 for last group
      alphas = c(B[(current_index+1):(current_index + mix_sizes[i] - 1)], 0)

      # reparameterize the beta vector
      B[current_index:(current_index + mix_sizes[i] - 1)] =
        (group_effect*exp(alphas)) / (exp(alphas) %>% sum)
    }

    # change the current index for next iteration
    current_index = current_index + mix_sizes[i]
  }
  return(B)
}

reparam_GI_weights = function(B, vars)
{
  # we organize into a list, then perform the reparameterization

  mix_sizes = sapply(vars$mixture, length) # sizes of mixture groups
  num_mixes = vars$mixture %>% length # number of mixtures

  current_index = 2 # 1st index is constant

  reparam = vector(mode = "list")

  for(i in 1:num_mixes) # iterate of number of mixtures
  {
    # group effect
    group_effect = sum(B[current_index:(current_index + mix_sizes[i] - 1)])
    reparam$group_index[i] = group_effect

    # weights
    reparam$weights[[paste("mix_", i, sep = "")]] =
      B[current_index:(current_index + mix_sizes[i] - 1)] /
      group_effect

    # change the current index for next iteration
    current_index = current_index + mix_sizes[i]
  }

  # if there are confounders
  if(length(union(vars$continuous, vars$categorical)) > 0)
  {
    reparam$confounders = B[current_index: length(B)]
  }

  reparam$constant = B[1]

  return(reparam)
}

get_optimization_region = function(logistic_param, vars)
{
  mix_sizes = sapply(vars$mixture, length) # sizes of mixture groups
  num_mixes = vars$mixture %>% length # number of mixtures

  current_index = 2 # 1st index is constant

  optimization_lower = rep(-Inf, logistic_param %>% length)
  optimization_upper = rep(Inf, logistic_param %>% length)

  for(i in 1:num_mixes) # iterate of number of mixtures
  {
    # find location of indices for each mixture group
    group_indices = current_index: (current_index + mix_sizes[i] - 1)
    # find the sign of the group
    group_sign = sign(logistic_param[group_indices]) %>% sum

    if(group_sign > 0)
    {
      optimization_lower[group_indices] = 0
      optimization_upper[group_indices] = Inf
    }else
    {
      optimization_lower[group_indices] = -Inf
      optimization_upper[group_indices] = 0
    }

    # change the current index for next iteration
    current_index = current_index + mix_sizes[i]
  }
  return(list(optimization_lower = optimization_lower,
              optimization_upper = optimization_upper))
}

fgwqsr_ll = function(B,design_matrix, y_vec, vars) # log likelihood
{
  B_logistic = make_beta_vec(B,vars) # to evaluate special Beta vector

  ll = logistic_neg_ll(B_logistic,design_matrix, y_vec, vars)
  # will return negative of log likelihood -- what we need for optim

  return(ll)
}

logistic_neg_ll = function(B,design_matrix, y_vec, vars) # log likelihood
{
  lin_pred = tcrossprod(design_matrix,matrix(B, nrow = 1))

  ll =ifelse(y_vec == 1, -lin_pred, lin_pred) %>% sigmoid::softplus %>% sum

  return(ll)
}
#
# fast_ll = function(B,design_matrix, y_vec, vars) # log likelihood
# {
#   # this function is for running one iteration of L-BFGS-B constrained optim
#   # split up log likelihood
#
#   ll_1 = (y_vec * crossprod(t(design_matrix), B)) %>% sum
#   ll_2 = (apply(design_matrix, MARGIN = 1, FUN = function(x) log(1+exp(crossprod(x, B))))) %>% sum
#
#   ll = ll_1 - ll_2
#
#   return(-1*ll) # -1 is for optim call
# }

logistic_neg_gr = function(B,design_matrix, y_vec, vars)
{

  gr = colSums(design_matrix * (c(y_vec - pracma::sigmoid(crossprod(t(design_matrix),B)))))

  return(-1*gr)
}

logistic_hessian = function(B,design_matrix)
{
  design_matrix = design_matrix %>% as.matrix

  z = as.matrix(design_matrix) %*% matrix(B, ncol = 1) # get z_i's for sigmoid function

  D = (pracma::sigmoid(z)*(1-pracma::sigmoid(z))) %>% as.vector

  hessian = -1*crossprod(design_matrix, (design_matrix * D))

  return(hessian)
}


reparam_to_weights = function(ML_sol, vars)
{
  # we need to recover the model weights
  num_mixes = length(vars$mixture)

  index_effects = ML_sol[which(substr(names(ML_sol), 1,2) == "B_")]

  constant = ML_sol[1]

  weight_list = list()

  for(i in 1:num_mixes)
  {
    group_weights = ML_sol[which(substring(names(ML_sol), 7) %in% vars$mixture[[i]])]

    num_elements = length(vars$mixture[[i]])

    group_index = 1:length(num_elements)

    denom = 1+ sum(sapply(group_weights, FUN = exp))

    weights = c(exp(group_weights) /denom, 1/denom )

    names(weights) = paste("w",i, 1:num_elements , sep = "")

    weight_list[[i]] = weights
  }

  length_phis = length(vars$continuous) + length(vars$categorical)

  return_list = NULL

  phis = c()

  if(length_phis >0)
  {
    phis = ML_sol[-(1:(length(unlist(vars$mixture)) + 1))] # find phi sols
    return_list = list(index_effects = index_effects, weight_list = weight_list, phis = phis, constant = constant)
  }else
  {
    return_list = list(index_effects = index_effects, weight_list = weight_list, constant = constant)
  }
  return(return_list)
}


fit_fgwqsr = function(formula, data, quantiles, output_hessian = F,
                      initial_cov_vals, return_y = F, return_data = F,
                      optim_control_list)
{
  f = as.character(formula)

  f[3] = gsub("\n", "", f[3])

  y = f[2]

  vars = clean_vars(f[3])

  data = data[, c(y, unlist(vars))] # keep only variables in model, order them for fast ll
  # NOTICE THAT DATA SET IS MODIFIED LOCALLY IN THIS FUNCTION CALL

  if(length(vars$categorical)>0) # if there are categorical variables
  {
    data = fastDummies::dummy_cols(data, select_columns = vars$categorical, remove_first_dummy = TRUE,
                      remove_selected_columns = TRUE) # add dummy variables
  }

  data = quantize_vars(data,vars, quantiles) # quantize mixture components

  # create likelihood model for ML

  num_confounders = ncol(data) - unlist(vars$mixture)%>%length -1 # -1 for outcome y column

  initial_vals =  rep(0, ncol(data))# also accounts for constant

  initial_vals[1] = initial_cov_vals$intercept_est # intercept estimate

  if(num_confounders > 0)
  {
    initial_vals[(ncol(data) - num_confounders + 1):ncol(data)] =
      initial_cov_vals$confounder_ests
  }
  # this is for likelihood ratio test

  y_vec = data[,1] # y outcome vector to send to likelihood, will always be first col
  design_matrix = cbind(rep(1, nrow(data)), data[, -1]) %>% as.matrix # the first column of data has y value

  # first perform a few initial iterations in reparameterization

  ML_sol = stats::optim(par = initial_vals,
                 fn = fgwqsr_ll,
                 design_matrix = design_matrix,
                 y_vec = y_vec,
                 vars = vars,
                 method = "BFGS",
                 control = list(maxit = 200,
                                factr = 1E-14,
                                reltol = 1E-14#,
                                #fnscale = fnscale
                 ))

  # now we extract the ML_sol and run using optim LBFGS for precise likelihood estimate

  # ML solution in logistic parameterization
  ML_sol_logistic_param = make_beta_vec(ML_sol$par, vars)

  optimization_region = get_optimization_region(ML_sol_logistic_param, vars)

  initial_vals = ML_sol_logistic_param

  final_fit = stats::optim(par = initial_vals,
                    fn = logistic_neg_ll,
                    gr = logistic_neg_gr,
                    design_matrix = design_matrix,
                    y_vec = y_vec,
                    vars = vars,
                    method = "L-BFGS-B",
                    lower = optimization_region$optimization_lower,
                    upper = optimization_region$optimization_upper,
                    control = optim_control_list)

  return_list = list(ML_sol= final_fit, vars = vars)

  # return y and new data only for fist FG formula.

  if(return_y == T)
  {
    return_list[["y"]] = y_vec
  }
  if(return_data == T)
  {
    return_list[["new_data"]] = data
  }
  return(return_list)
}

generate_optimization_regions = function(vars, num_confounders)
{
  perms = gtools::permutations(2, r = length(vars$mixture), v = c(-1,1),
                       repeats.allowed = T)
  num_mixes = length(vars$mixture)

  num_in_each_group = sapply(vars$mixture, length)

  opt_list = apply(perms, MARGIN = 1, FUN = function(x)
  {
    # +1 is for constant at beginning
    lower = rep(-Inf, sum(num_in_each_group) + num_confounders + 1)
    upper = rep(Inf, sum(num_in_each_group) + num_confounders+ 1)
    #perm_vec = rep(x, times = num_in_each_group))

    # iterate over groups for each perm
    current_index = 2
    for(i in 1:num_mixes)
    {
      # if there is not a single pollutant in the group -- we would want to leave
      # that single pollutant unconstrained
      if(num_in_each_group[i] != 1)
      {
        if(x[i] == 1)
        {
          lower[current_index:(current_index + num_in_each_group[i] - 1)] = 0
          upper[current_index:(current_index + num_in_each_group[i] - 1)] = Inf
        }else
        {
          lower[current_index:(current_index + num_in_each_group[i] - 1)] = -Inf
          upper[current_index:(current_index + num_in_each_group[i] - 1)] = 0
        }

      }

      current_index = current_index + num_in_each_group[i]

    }

    return(list(lower = lower, upper = upper))


  })
  return(opt_list)
}

get_cov_initial_vals = function(formula, data, quantiles)
{
  # get initial values for intercept and confounder estimates using quantile
  # logistic regression

  f = as.character(formula)

  f[3] = gsub("\n", "", f[3])

  y = f[2]

  vars = clean_vars(f[3])

  data = data[, c(y, unlist(vars))] # keep only variables in model, order them for fast ll
  # NOTICE THAT DATA SET IS MODIFIED LOCALLY IN THIS FUNCTION CALL

  if(length(vars$categorical)>0) # if there are categorical variables
  {
    data = fastDummies::dummy_cols(data, select_columns = vars$categorical, remove_first_dummy = TRUE,
                      remove_selected_columns = TRUE) # add dummy variables
  }

  num_confounders = ncol(data) - unlist(vars$mixture)%>%length -1 # -1 for outcome y column

  data = quantize_vars(data,vars, quantiles) # quantize mixture components

  # run GLM -- we use all variables in data matrix since we only kept used vars

  glm_formula = paste(y, "~ .") %>% stats::as.formula

  unconstr_glm = glm(formula =glm_formula, data = data, family = binomial)

  glm_coefs = coef(unconstr_glm)

  intercept_est = glm_coefs[1] # get intercept estimate

  confounder_ests = NA

  if(num_confounders > 0) # if there are confounders, get the estimates
  {
    confounder_ests = glm_coefs[-c(1:(length(glm_coefs) - num_confounders))]
  }

  initial_ests = list(intercept_est = intercept_est,
                      confounder_ests = confounder_ests)
  return(initial_ests)

}

fgwqsr_caller = function(formulas, data, quantiles,vars, verbose, cores, optim_control_list)
{
  if(verbose == TRUE) # call to progress
  {
    p <- progressr::progressor(along = formulas)
  }

  # get initial vals for intercept and adjusting covariates
  initial_cov_vals = get_cov_initial_vals(formula = stats::as.formula(formulas[[1]]),
                                          data = data,
                                          quantiles = quantiles)
  future::plan(future::multisession, workers = cores)

  fits = future.apply::future_lapply(1:length(formulas), FUN = function(x)
  {
    result = NULL
    if(x == 1) # if the original formula
    {
      result = fit_fgwqsr(stats::as.formula(formulas[x]), data, quantiles, return_y = T,
                          return_data = T, initial_cov_vals = initial_cov_vals,
                          optim_control_list = optim_control_list)
      if(verbose == TRUE)
      {
        p(paste("Fitted Model ", x, " of ", length(formulas), sep = ""))
      }
    }else # if LRT for weights, return LL
    {
      # the ifelse in "result" is to handle the ll for the 1 group null ll for a 1 group lrt
      result = ifelse( (length(vars$mixture) == 1) && (x == (vars$mixture %>% unlist %>% length %>% sum(2))),
                       stats::glm(stats::as.formula(paste(formulas[x], "1")), data = data, family = "binomial") %>% stats::logLik,
                       -1*fit_fgwqsr(stats::as.formula(formulas[x]), data, quantiles,
                                     initial_cov_vals = initial_cov_vals,
                                     optim_control_list = optim_control_list)$ML_sol$value)

      if(verbose == TRUE)
      {
        p(paste("Fitted Model ", x, " of ", length(formulas), sep = ""))
      }
    }
    return(result)

  })

  future::plan(future::sequential)

  return(fits)
}

#' Fit a FGWQSR Model
#'
#' @param formula A formula for model fitting of FGWQSR.  Please see description for formula construction
#' @param data R dataframe that contains all covariates and the outcome data
#' @param quantiles number of quantiles to quantize the exposure variables in the mixture portion of the model.  Default value is 5.
#' @param n_mvn_sims defines resolution for simulated null distribution for group index and single chemical LRTs.  Default is 10,000.
#' @param zero_threshold_cutoff defines a .
#' @param verbose Displays messages and progress bar while fitting FGWQSR model.  Default is TRUE.
#' @param cores number of cores to parallelize on for fitting nested models and simulated null LRT distributions.  Default is number of available cores on user device.
#' @param optim_control_list - option to supply control options to optim.
#' @return list with attributes from fgwqsr model fitting.
#' @import statar
#' @import progressr
#' @export

fgwqsr = function(formula, data, quantiles = 5, n_mvn_sims = 10000,
                  zero_threshold_cutoff = .5, verbose = T,
                  cores = future::availableCores(),
                  optim_control_list = list(maxit = 1000, factr = 1E-12, fnscale = 1))
{

  time_begin = proc.time()[3]

  # get formula and vars object that stores info about model
  f = as.character(formula); vars = clean_vars(gsub("\n", "", f[3]))

  # get list of formulas to run submodels for LRT
  formulas = c(c(paste(format(formula), collapse = ""),get_formulas(vars,formula)))

  # run models in parallel, verbose option
  if(verbose == TRUE)
  {
    progressr::handlers("progress")
    message("\nFitting nested models for Weight Inference:")
    progressr::with_progress(fits <- fgwqsr_caller(formulas, data, quantiles,vars, verbose,
                                        cores, optim_control_list))
    message("\nNow Performing Inference...")
  }else
  {
    fits <- fgwqsr_caller(formulas, data, quantiles,vars, verbose, cores, optim_control_list)
  }

  # full model
  fgwqsr_fit = fits[[1]]

  # ll from submodels
  ll_models = c(-1*fgwqsr_fit$ML_sol$value, unlist(fits[-1]))

  # model attributes of original formula
  vars = fgwqsr_fit$vars

  # outcome variable vector
  y = fgwqsr_fit$y

  # new data with categorical variables dummified and quantile exposures
  new_data = fgwqsr_fit$new_data

  # final ML solution in logistic regression formulation
  params_logistic_form = fgwqsr_fit$ML_sol$par

  # compute inverse of observed fisher info at constrained logistic regression
  # ML solution
  design_matrix = cbind(intercept = rep(1, nrow(new_data)),new_data[,-1])
  observed_fisher = -1*logistic_hessian(B = params_logistic_form,
                                        design_matrix = design_matrix)
  cov_mat = pracma::pinv(observed_fisher)

  # get ML sol in group index and weight formulation
  names(params_logistic_form) = names(design_matrix)
  reparam_reg = reparam_GI_weights(params_logistic_form, vars)

  # now we perform inference

  # call function to perform inference #####################################

  inference_frames = NULL # initialize

  if(verbose == T)
  {
    with_progress(inference_frames <- perform_inference(ll_models = ll_models,
                                                        params_logistic_form = params_logistic_form,
                                                        vars = vars, cov_mat = cov_mat,
                                                        zero_threshold_cutoff = zero_threshold_cutoff,
                                                        num_sims = n_mvn_sims,
                                                        reparam_reg = reparam_reg,
                                                        cores = cores,
                                                        verbose = verbose))
  }else
  {
    inference_frames = perform_inference(ll_models = ll_models,
                                         params_logistic_form = params_logistic_form,
                                         vars = vars, cov_mat = cov_mat,
                                         zero_threshold_cutoff = zero_threshold_cutoff,
                                         num_sims = n_mvn_sims,
                                         reparam_reg = reparam_reg,
                                         cores = cores,
                                         verbose = verbose)
  }

  ##########################################################################




  total_time = proc.time()[3]-time_begin

  # for summary method

  n = nrow(new_data)

  full_ll_val = ll_models[1]

  aic = -2*full_ll_val + 2*length(params_logistic_form)

  bic = -2*full_ll_val + log(n)*length(params_logistic_form)

  # when certain estimates are on the boundary, convergance code issues
  # message regarding LNSRCH == line search method.  Therefore, do not look at
  # convergence code in some instances

  L_BFGS_B_convergance = list(counts = fgwqsr_fit$ML_sol$counts,
                              convergence = fgwqsr_fit$ML_sol$convergence,
                              message = fgwqsr_fit$ML_sol$message)

  return(list(inference_frames = inference_frames,
              total_time = total_time,
              n = n,
              ll = full_ll_val,
              formula = formula,
              aic = aic,
              bic = bic,
              vars = vars,
              n_mvn_sims = n_mvn_sims,
              cores = cores,
              L_BFGS_B_convergance =L_BFGS_B_convergance,
              param_cov_mat = cov_mat))
}



#' Summarize a fgwqsr model fit
#' @param fgwqsr_sol a fitted object from a fgwqsr() call
#' @param digits the number of rounding digits to display in summary tables.
#' @export
fgwqsr_summary = function(fgwqsr_sol, digits  = 6)
{
  # rounding for mixture index frame
  fgwqsr_sol$inference_frames$group_index_frame[,1:2] =
    round(fgwqsr_sol$inference_frames$group_index_frame[,1:2], digits = digits)

  # fix pvalues
  fgwqsr_sol$inference_frames$group_index_frame[,3] =
    sapply(fgwqsr_sol$inference_frames$group_index_frame[,3],
           FUN = format_scientific, cutoff = 10^(-digits))

  # rounding for weight frame

  fgwqsr_sol$inference_frames$weight_frame[,1:2] =
    round(fgwqsr_sol$inference_frames$weight_frame[,1:2], digits = digits)

  # fix pvalues
  fgwqsr_sol$inference_frames$weight_frame[,3] =
    sapply(fgwqsr_sol$inference_frames$weight_frame[,3],
           FUN = format_scientific, cutoff = 10^(-digits))

  # fix adjusting covariates

  fgwqsr_sol$inference_frames$adj_param_frame[,1:3] =
    round(fgwqsr_sol$inference_frames$adj_param_frame[,1:3], digits = digits)

  # fix pvalues
  fgwqsr_sol$inference_frames$adj_param_frame[,4] =
    sapply(fgwqsr_sol$inference_frames$adj_param_frame[,4],
           FUN = format_scientific, cutoff = 10^(-digits))

  # make confidence interval rounded for adjusting covariates
  fgwqsr_sol$inference_frames$adj_param_frame[, 5:6] =
    round(fgwqsr_sol$inference_frames$adj_param_frame[, 5:6], digits = digits)

  ci_95 = paste("(",fgwqsr_sol$inference_frames$adj_param_frame[, 5], ", ",
                fgwqsr_sol$inference_frames$adj_param_frame[, 6], ")", sep = "")

  fgwqsr_sol$inference_frames$adj_param_frame[,5] = ci_95
  colnames(fgwqsr_sol$inference_frames$adj_param_frame)[5] = "95% CI"
  fgwqsr_sol$inference_frames$adj_param_frame =
    subset(fgwqsr_sol$inference_frames$adj_param_frame, select = -6) # remove column with upper 95% endpoint

  cat("\nCall: \nFGWQSR with formula '",gsub("  ", "",paste(format(fgwqsr_sol$formula), collapse = "")),"' on n = ",fgwqsr_sol$n," observations.", sep = "" )
  cat("\n\n", fgwqsr_sol$n_mvn_sims, " samples used for simulated LRT distirbution.",sep = "")
  cat("\n\nLog Likelihood:", fgwqsr_sol$ll, "| AIC:",fgwqsr_sol$aic, "| BIC:", fgwqsr_sol$bic)
  cat("\n\nEstimates and Inference for Group Index Effects\n", sep = "")
  print(fgwqsr_sol$inference_frames$group_index_frame, digits = digits)
  cat("\nEstimates and Inference for Weights\n")

  current_index = 1
  for(i in 1: length(fgwqsr_sol$vars$mixture))
  {
    group_size = fgwqsr_sol$vars$mixture[[i]] %>% length
    output = fgwqsr_sol$inference_frames$weight_frame[current_index: (current_index + group_size - 1), ]
    print(output, digits = digits)
    current_index = current_index + group_size
    cat("-------------------------------------------------\n")
  }
  cat("\nEstimates and Inference for Intercept and Adjusting Covariates\n")
  print(fgwqsr_sol$inference_frames$adj_param_frame)
  cat("\nSignificance Codes: <0.001 '***' <0.01 '**' <0.05 '*' <0.10 '.' \n")
  cat("\nTotal runtime for FGWQSR: ",
      ifelse(fgwqsr_sol$total_time < 60,paste(round(fgwqsr_sol$total_time,2), "seconds"),
             paste(round(fgwqsr_sol$total_time/60, 2), "minutes")), "on",fgwqsr_sol$cores, "cores." )
}


