# Main functions for FGWQSR

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


# create_likelihood_string = function(data, vars)
# {
#   ### Beta formulation##########################################################
#
#   beta_string = "B[1]"
#
#   beta_index = 2 # to hold value where group index effect is in vector
#
#   num_mixes = length(vars$mixture) # number of mixtures
#
#   beta_orders_mix = character(1+length(unlist(vars$mixture)))
#
#   beta_orders_mix[1] = "constant"
#
#   for(i in 1: length(vars$mixture)) # iterate over number of mixtures
#   {
#     num_elements = length(vars$mixture[[i]])
#
#     group_string = ""
#
#     if(num_elements>1) # if more than 1 element in mixture group
#     {
#       coef_denoms = paste("(1 + ", paste("exp(", paste("B[", (beta_index+1):(beta_index+num_elements -1),"]",
#                                                        sep = ""), ")", sep = "", collapse = " + "),")", sep = "") # denominators of coefs
#
#       coef_nums = c(paste("exp(B[", (beta_index+1):(beta_index+num_elements - 1),"])",
#                           sep = ""), "1") # numerators of coefficients
#
#       coeffs = paste("(",coef_nums, "/", coef_denoms, ")", sep = "") # vector of weight coefficient strings
#
#       quant_vars = paste("data$", vars$mixture[[i]], "[i]", sep = "")
#
#       group_string = paste("B[", beta_index, "]*(", paste(coeffs, quant_vars, sep = "*", collapse = " + "),
#                            ")", sep = "")
#
#       beta_string = paste(beta_string, group_string, sep = " + ")
#
#       beta_orders_mix[beta_index:(beta_index+num_elements -1)] = c(
#         paste("B_mix_", i, sep = ""), paste("alpha_", vars$mixture[[i]][-length(vars$mixture[[i]])], sep = ""))
#
#     }else # if one element in mixture group, weight is 1
#     {
#       quant_vars = paste("data$", vars$mixture[[i]], "[i]", sep = "")
#
#       group_string = paste("B[", beta_index, "]*",quant_vars, sep = "")
#
#       beta_string = paste(beta_string, group_string, sep = " + ")
#
#       beta_orders_mix[beta_index:(beta_index+num_elements -1)] = paste("B_mix_", i, sep = "")
#     }
#
#     beta_index = beta_index + length(vars$mixture[[i]])
#   }
#
#   ### For Phis #################################################################
#
#   phi_start_index = beta_index
#
#   length1 = length(vars$continuous)
#
#   length2 = length(vars$categorical)
#
#   phi1 = c()
#
#   if(length1>0)
#   {
#     for(i in 1: length1)
#     {
#       phi1[i] = paste(vars$continuous[i], sep = "")
#     }
#   }
#
#   phi2 = list()
#
#   if(length2>0)
#   {
#     for(i in 1:length2)
#     {
#       varName = vars$categorical[i]
#
#       names = names(data)[which(substr(names(data), 1, nchar(varName)+1) ==
#                                   paste(varName, "_", sep = ""))]
#       phi2[[i]] = paste(names, sep= "")
#
#     }
#   }
#
#   phis = c(phi1, unlist(phi2))
#
#   phi_betas = c()
#
#   if(length(phis)>0) # if confounders, add betas
#   {
#     phi_betas = paste("B[", (phi_start_index):(phi_start_index + length(phis) -1), "]", sep = "")
#   }
#
#
#   order_of_betas = c(beta_orders_mix, phis)
#
#
#
#   phiString = c()
#
#   if(length(phis)>0)
#     phiString = paste(" + ", paste(phi_betas, "*", "data$", phis,"[i]", sep = "", collapse = " + "))
#
#   mixString = paste(beta_string, phiString, sep = "")
#
#   return(list(mixString=mixString, order_of_betas=order_of_betas, numPhis = length(phis)))
# }
#
make_beta_vec = function(B, vars)
{
  # we organize into a list, then perform the reparameterization

  mix_sizes = sapply(vars$mixture, length) # sizes of mixture groups
  num_mixes = vars$mixture %>% length # number of mixtures

  current_index = 2 # 1st index is constant

  for(i in 1:num_mixes) # iterate of number of mixtures
  {
    # if mix_size is 1, just leave param for unconstrained optimization
    if(mix_sizes[i] > 1)
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

# generate initial values for nested models
generate_nested_initial_vals = function(full_mod_sol, vars)
{
  # need to handle when there is a 1 group LRT

  # full_mod_sol is organized as (intercept, chems, ..., adjusting covs)

  # store initial values
  num_formulas = length(vars$mixture) + sum(sapply(vars$mixture, length))
  reparam_initial_vals = vector(mode = "list", length = num_formulas)

  # first single pollutant initial values
  run_ind = 1
  for(i in 1: length(vars$mixture))
  {
    temp_vars = lapply(1:length(vars$mixture[[i]]),
                       FUN = function(j)
                       {
                         v_temp = vars
                         v_temp$mixture[[i]] = v_temp$mixture[[i]][-j]
                         return(v_temp)
                       })

    exclude_index = (run_ind+1):(run_ind + length(vars$mixture[[i]]))

    initial_vals= lapply(1:length(temp_vars), FUN = function(j)
        {
        return(full_mod_sol[-exclude_index[j]])
      })

    # reparameterize to multinomial parameterization

    reparam_initial_vals[run_ind: (run_ind + length(vars$mixture[[i]]) -1)] =
      lapply(1:length(temp_vars), FUN = function(j)
      return(reparam_to_multinomial(initial_vals[[j]], temp_vars[[j]])))

    # increment running indices
    run_ind = run_ind + length(vars$mixture[[i]])
  }

  # add initial values for group LRTs
  if(length(vars$mixture) > 1) # if more than 1 mixture group, then we can put initial values
  {
    counter = 2
    for(i in 1:length(vars$mixture))
    {
      temp_vars = vars
      temp_vars$mixture = temp_vars$mixture[-i]
      exclude_index = counter:(counter + length(vars$mixture[[i]]) - 1)
      reparam_initial_vals[[run_ind]] = reparam_to_multinomial(full_mod_sol[-exclude_index], temp_vars)
      counter = counter + length(vars$mixture[[i]])
      run_ind = run_ind + 1
    }
  }
  return(reparam_initial_vals)
}


# this function is for the initial values of the nested models, for the hybridized fitting
reparam_to_multinomial = function(B, temp_vars)
{
  num_groups = length(temp_vars$mixture)
  num_in_each_group = sapply(temp_vars$mixture, length)

  reparam = B # initialize

  current_ind = 2
  for(i in 1: num_groups)
  {
    # if only 1 in the group, leave for unconstrained optimziation
    if(num_in_each_group[i] > 1)
    {
      # betas for group
      betas = B[current_ind:(current_ind + num_in_each_group[i] -1)]
      gamma = sum(betas)

      # if there are numerically 0 values, need them to not be zero for reparameterization
      # ALSO -- we set small betas to 1E-3, so that solving system is numerically stable!
      # these are just to solve for initial values so okay if off.
      if(gamma > 0) # if group effect is positive
      {
        betas = ifelse(betas <= 1E-3, 1E-3, betas)
      }else if(gamma< 0) # if group effect is negative
      {
        betas = ifelse(betas >= -1E-3 , -1E-3, betas)
      }else
        # if single pollutant with all the weight from the group is removed,
        # choose at random to make effects 1E-3 or -1E-3 -- all betas have same sign though.
      {
        betas = rep( sample(c(-1,1), 1) * 1E-3, length(betas))
      }

      # recompute gamma
      gamma = sum(betas)

      # need to solve for alphas
      # create matrix to solve system of equations
      mat = matrix(data = 1, nrow = num_in_each_group[i] - 1,
                   ncol = num_in_each_group[i] - 1)

      # add some numerical jitter to diagonal to avoid computational singularity
      diag(mat) = 1 - gamma/betas[-length(betas)] + abs(stats::rnorm(n = length(betas)-1, mean = 1E-6, sd = 0.0001))

      # solve for alphas
      alphas = solve(mat, rep(-1,num_in_each_group[i] - 1)) %>% log

      # return (gamma, alphas) for group
      reparam[current_ind:(current_ind + num_in_each_group[i] -1)] = c(gamma, alphas)
    }
    # increment counter
    current_ind = current_ind + num_in_each_group[i]
  }
  return(reparam)
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

fgwqsr_logistic_ll = function(B,design_matrix, y_vec, vars) # log likelihood
{
  B_logistic = make_beta_vec(B,vars) # to evaluate special Beta vector

  ll = logistic_neg_ll(B_logistic,design_matrix, y_vec)
  # will return negative of log likelihood -- what we need for optim

  return(ll)
}

fgwqsr_poisson_ll = function(B,design_matrix, y_vec, vars) # log likelihood
{
  new_B = make_beta_vec(B,vars) # to evaluate special Beta vector

  ll = poisson_neg_ll_WO_factorial(new_B,design_matrix, y_vec)
  # will return negative of log likelihood -- what we need for optim
  return(ll)
}

fgwqsr_ols_loss = function(B,design_matrix, y_vec, vars) # log likelihood
{
  new_B = make_beta_vec(B,vars) # to evaluate special Beta vector

  loss = ols_loss(new_B,design_matrix, y_vec)

  return(loss)
}

logistic_neg_ll = function(B,design_matrix, y_vec) # log likelihood
{
  lin_pred = tcrossprod(design_matrix,matrix(B, nrow = 1))

  ll =ifelse(y_vec == 1, -lin_pred, lin_pred) %>% sigmoid::softplus() %>% sum

  return(ll)
}

poisson_neg_ll_WO_factorial = function(B,design_matrix, y_vec) # log likelihood
{
  lin_pred = tcrossprod(design_matrix,matrix(B, nrow = 1))

  ll = sum(exp(lin_pred) - y_vec * lin_pred)

  return(ll)
}

poisson_neg_ll = function(B,design_matrix, y_vec) # log likelihood
{
  lin_pred = tcrossprod(design_matrix,matrix(B, nrow = 1))

  ll = sum(exp(lin_pred) - y_vec * lin_pred + lfactorial(y_vec))

  return(ll)
}

ols_loss = function(B,design_matrix, y_vec) # log likelihood
{
  lin_pred = tcrossprod(design_matrix,matrix(B, nrow = 1))
  loss = as.numeric(crossprod(y_vec - lin_pred, y_vec - lin_pred))
  return(loss)
}

ols_gr = function(B,design_matrix, y_vec) # log likelihood
{
  lin_pred = tcrossprod(design_matrix,matrix(B, nrow = 1))
  gr =as.numeric(-2*crossprod(y_vec - lin_pred, design_matrix))
  return(gr)
}

mvn_neg_ll = function(B,design_matrix, y_vec, sigma_2) # log likelihood
{
  n = nrow(design_matrix)

  lin_pred = tcrossprod(design_matrix,matrix(B, nrow = 1))

  ll = -(n/2)*log(2 * pi) - (n/2)*log(sigma_2) -
    (1/(2*sigma_2)) * crossprod(y_vec -lin_pred, y_vec -lin_pred)

  return(-1*ll)
}

logistic_neg_gr = function(B,design_matrix, y_vec)
{
  gr = colSums(design_matrix * (c(y_vec - pracma::sigmoid(crossprod(t(design_matrix),B)))))
  return(-1*gr)
}

poisson_neg_gr = function(B,design_matrix, y_vec)
{
  lin_pred = tcrossprod(design_matrix,matrix(B, nrow = 1))

  gr = as.vector(crossprod(design_matrix, matrix(y_vec, ncol = 1)) - crossprod(design_matrix, exp(lin_pred)))

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

poisson_hessian = function(B,design_matrix)
{
  design_matrix = design_matrix %>% as.matrix

  lin_pred = tcrossprod(design_matrix,matrix(B, nrow = 1))

  W = lin_pred %>% exp %>% as.numeric %>% diag

  hessian = -1*crossprod(design_matrix, W) %*% design_matrix

  return(hessian)
}

gaussian_hessian = function(B,design_matrix, sigma_2)
{
  design_matrix = design_matrix %>% as.matrix

  hessian = -(1/sigma_2) * crossprod(design_matrix, design_matrix)

  return(hessian)
}

#
# reparam_to_weights = function(ML_sol, vars)
# {
#   # we need to recover the model weights
#   num_mixes = length(vars$mixture)
#
#   index_effects = ML_sol[which(substr(names(ML_sol), 1,2) == "B_")]
#
#   constant = ML_sol[1]
#
#   weight_list = list()
#
#   for(i in 1:num_mixes)
#   {
#     group_weights = ML_sol[which(substring(names(ML_sol), 7) %in% vars$mixture[[i]])]
#
#     num_elements = length(vars$mixture[[i]])
#
#     group_index = 1:length(num_elements)
#
#     denom = 1+ sum(sapply(group_weights, FUN = exp))
#
#     weights = c(exp(group_weights) /denom, 1/denom )
#
#     names(weights) = paste("w",i, 1:num_elements , sep = "")
#
#     weight_list[[i]] = weights
#   }
#
#   length_phis = length(vars$continuous) + length(vars$categorical)
#
#   return_list = NULL
#
#   phis = c()
#
#   if(length_phis >0)
#   {
#     phis = ML_sol[-(1:(length(unlist(vars$mixture)) + 1))] # find phi sols
#     return_list = list(index_effects = index_effects, weight_list = weight_list, phis = phis, constant = constant)
#   }else
#   {
#     return_list = list(index_effects = index_effects, weight_list = weight_list, constant = constant)
#   }
#   return(return_list)
# }

# Hybrdized approach had problem -- weights that were estimated to be 0
# had LRTS that were nonzero.  The reason this occurred was because
# the initial optimization, in the reparameterization, was not converging
# in a particular subregion for a group, and the LLs between full and nested
# models when exluding a pollutant with weight 0 was nonzero.

# this happens on few occasions, we can manually make the LRT
fit_fgwqsr_hybrid = function(formula, data, quantiles,family, output_hessian = F,
                      initial_vals, return_y = F, return_data = F,
                      optim_control_list)
{
  f = as.character(formula)

  f[3] = gsub("\n", "", f[3])

  y = f[2]

  vars = clean_vars(f[3])

  # this dataset is modified locally.  Different LRT formulas specify exclusion
  # of different variables from the dataframe.  Each local data is built to
  # have the data organized such that it can be called in likelihood function.
  data = data[, c(y, unlist(vars))]

  if(length(vars$categorical)>0) # if there are categorical variables
  {
    data = fastDummies::dummy_cols(data, select_columns = vars$categorical, remove_first_dummy = TRUE,
                      remove_selected_columns = TRUE) # add dummy variables
  }

  data = quantize_vars(data,vars, quantiles) # quantize mixture components

  # create likelihood model for ML

  num_confounders = ncol(data) - unlist(vars$mixture)%>%length -1 # -1 for outcome y column

  # initial_vals =  rep(0, ncol(data))# also accounts for constant
  #
  # initial_vals[1] = initial_cov_vals$intercept_est # intercept estimate
  #
  # if(num_confounders > 0)
  # {
  #   initial_vals[(ncol(data) - num_confounders + 1):ncol(data)] =
  #     initial_cov_vals$confounder_ests
  # }
  # this is for likelihood ratio test

  y_vec = data[,1] # y outcome vector to send to likelihood, will always be first col
  design_matrix = cbind(rep(1, nrow(data)), data[, -1]) %>% as.matrix # the first column of data has y value

  fn_hyb = fgwqsr_logistic_ll
  if(family == "poisson")
  {
    fn_hyb = fgwqsr_poisson_ll
  }else if(family == "gaussian") # NEED TO MODIFY THIS!!!
  {
    fn_hyb = fgwqsr_ols_loss
  }


  # first perform a few initial iterations in reparameterization

  ML_sol = stats::optim(par = initial_vals,
                 fn = fn_hyb,
                 design_matrix = design_matrix,
                 y_vec = y_vec,
                 vars = vars,
                 method = "BFGS",
                 control = list(maxit = 1000,
                                factr = 1E-14,
                                reltol = 1E-20#,
                                #fnscale = fnscale
                 ))

  # now we extract the ML_sol and run using optim LBFGS for precise likelihood estimate

  # ML solution in logistic parameterization
  ML_sol_logistic_param = make_beta_vec(ML_sol$par, vars)

  optimization_region = get_optimization_region(ML_sol_logistic_param, vars)

  initial_vals = ML_sol_logistic_param

  # adjust log likelihood and gradient by family
  # assume first logistic regression
  fn = logistic_neg_ll
  gr = logistic_neg_gr
  if(family == "poisson")
  {
    fn = poisson_neg_ll_WO_factorial
    gr = poisson_neg_gr
  }else if(family == "gaussian") # NEED TO MODIFY THIS!!!
  {
    fn = ols_loss
    gr = ols_gr
  }

  final_fit = NULL

  # set these -- adjusted in loop
  fn_scale = 1
  counter = 1

  max_iter = 30

  # loop in case solution does not converge
  # if does not converge,
  while(is.null(final_fit))
  {

    final_fit = tryCatch(expr =
                           {
                             stats::optim(par = initial_vals,
                                          fn = fn,
                                          gr = gr,
                                          design_matrix = design_matrix,
                                          y_vec = y_vec,
                                          method = "L-BFGS-B",
                                          lower = optimization_region$optimization_lower,
                                          upper = optimization_region$optimization_upper,
                                          control = optim_control_list)
                           },
                         error = function(err)
                         {
                           NULL
                         })
    if(counter == max_iter)
    {
      stop("FGWQSR failed to converge.")
    }

    # increment scale and counter
    fn_scale = 1*10^(counter)
    counter = counter + 1
  }




  # if gaussian, compute sigma^2 and calculate log likelihood
  if(family == "gaussian")
  {
    sigma_2 = sum((y_vec - tcrossprod(design_matrix,matrix(final_fit$par, nrow = 1)))^2) /
      (nrow(design_matrix) - ncol(design_matrix))

    # change loss value to log likelihood value
    final_fit$value = mvn_neg_ll(final_fit$par, design_matrix, y_vec, sigma_2)
  }else if(family == "poisson") # evaluate ll with factorial term
  {
    final_fit$value = poisson_neg_ll(final_fit$par, design_matrix, y_vec)
  }

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

# MO stands for multiple optim
fit_fgwqsr_MO = function(formula, data, quantiles,family, output_hessian = F,
                      initial_cov_vals, return_y = F, return_data = F,
                      optim_control_list, cores)
{
  f = as.character(formula)

  f[3] = gsub("\n", "", f[3])

  y = f[2]

  vars = clean_vars(f[3])

  # this dataset is modified locally.  Different LRT formulas specify exclusion
  # of different variables from the dataframe.  Each local data is built to
  # have the data organized such that it can be called in likelihood function.
  data = data[, c(y, unlist(vars))]

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


  optimization_regions = generate_optimization_regions(vars = vars,
                                                       num_confounders = num_confounders)

  # adjust log likelihood and gradient by family
  # assume first logistic regression
  fn = logistic_neg_ll
  gr = logistic_neg_gr
  if(family == "poisson")
  {
    fn = poisson_neg_ll_WO_factorial
    gr = poisson_neg_gr
  }else if(family == "gaussian")
  {
    fn = ols_loss
    gr = ols_gr
  }

  future::plan(future::multisession, workers = cores)

  fits = future.apply::future_lapply(1: length(optimization_regions),
                                     FUN = function(i)
                                     {
                                       initial_vals[2:(ncol(data) - num_confounders)] = generate_initial_vals(optimization_regions[[i]]$lower,
                                                                                                              optimization_regions[[i]]$upper)[2:(ncol(data) - num_confounders)]
                                       sol = NULL

                                       # set these -- adjusted in loop
                                       fn_scale = 1
                                       counter = 1

                                       max_iter = 30

                                       while(is.null(sol) && counter <= max_iter )
                                       {

                                         sol = tryCatch(expr =
                                                          {
                                                            stats::optim(par = initial_vals,
                                                                         fn = fn,
                                                                         gr = gr,
                                                                         design_matrix = design_matrix,
                                                                         y_vec = y_vec,
                                                                         method = "L-BFGS-B",
                                                                         lower = optimization_regions[[i]]$lower,
                                                                         upper = optimization_regions[[i]]$upper,
                                                                         control = c(optim_control_list, fnscale = fn_scale))
                                                          },
                                                        error = function(err)
                                                        {
                                                          NULL
                                                        })

                                         fn_scale = 1*10^(counter)
                                         counter = counter + 1
                                       }

                                       if(counter == max_iter)
                                       {
                                         warning("FGWQSR failed to converge in a particular subregion.
                                                   \nProceeding without evaluating likelihood in this region.")
                                       }
                                       return(sol)
                                     })

  future::plan(future::sequential)

  fits = fits[lengths(fits) != 0]

  final_fit = fits[[which.min(sapply(fits, "[[", 2))]]

  # if gaussian, compute sigma^2 and calculate log likelihood
  if(family == "gaussian")
  {
    sigma_2 = sum((y_vec - tcrossprod(design_matrix,matrix(final_fit$par, nrow = 1)))^2) /
      (nrow(design_matrix) - ncol(design_matrix))

    # change loss value to log likelihood value
    final_fit$value = mvn_neg_ll(final_fit$par, design_matrix, y_vec, sigma_2)
  }else if(family == "poisson") # evaluate ll with factorial term
  {
    final_fit$value = poisson_neg_ll(final_fit$par, design_matrix, y_vec)
  }

  return_list = list(ML_sol= final_fit, vars = vars)

  if(family == "gaussian")
  {
    return_list$residual_var = sum((y_vec - tcrossprod(design_matrix,matrix(final_fit$par, nrow = 1)))^2) /
      (nrow(design_matrix) - ncol(design_matrix))
  }

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

get_cov_initial_vals = function(formula, data, quantiles, family)
{
  # get initial values for intercept and confounder estimates using quantile
  # logistic regression

  f = as.character(formula)

  f[3] = gsub("\n", "", f[3])

  y = f[2]

  vars = clean_vars(f[3])

  data = data[, c(y, unlist(vars))]

  if(length(vars$categorical)>0) # if there are categorical variables
  {
    data = fastDummies::dummy_cols(data, select_columns = vars$categorical, remove_first_dummy = TRUE,
                      remove_selected_columns = TRUE) # add dummy variables
  }

  num_confounders = ncol(data) - unlist(vars$mixture)%>%length -1 # -1 for outcome y column

  data = quantize_vars(data,vars, quantiles) # quantize mixture components

  # run GLM -- we use all variables in data matrix since we only kept used vars

  glm_formula = paste(y, "~ .") %>% stats::formula()

  unconstr_glm = stats::glm(formula =glm_formula, data = data, family = family)

  glm_coefs = stats::coef(unconstr_glm)

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

fgwqsr_caller = function(formulas, data, quantiles, family, vars, verbose, cores, optim_control_list)
{
  if(verbose == TRUE) # call to progress
  {
    p <- progressr::progressor(steps = length(formulas))
  }

  # get initial vals for intercept and adjusting covariates
  initial_cov_vals = get_cov_initial_vals(formula = stats::formula(formulas[[1]]),
                                          data = data,
                                          quantiles = quantiles,
                                          family = family)
  # initialize for storing results
  fits = vector(mode = "list", length = length(formulas))

  # perform full model fit using multiple optim -- not hybridized!

  fits[[1]] = fit_fgwqsr_MO(formula = stats::formula(formulas[[1]]),
                            data = data,
                            quantiles = quantiles,
                            family = family,
                            return_y = T,
                            return_data = T, initial_cov_vals = initial_cov_vals,
                            optim_control_list = optim_control_list,
                            cores = cores)

  if(verbose == T)
  {
    p()
  }

  # From this first fit, get solution -- generate initial values for nested models

  full_mod_sol = fits[[1]]$ML_sol$par

  # find which weights were estimated to be 0 -- do NOT waste time fitting
  # SPLRT for chemicals with weights of 0

  single_chem_effects = full_mod_sol[2: (2+sum(sapply(vars$mixture, length)) - 1)]

  formulas[2:(2+sum(sapply(vars$mixture, length)) - 1)] = ifelse(single_chem_effects == 0, "0 weight",
                                                                 formulas[2:(2+sum(sapply(vars$mixture, length)) - 1)] )

  # generate initial values for nested models -- initial values for chemicals too!

  initial_vals = generate_nested_initial_vals(full_mod_sol = full_mod_sol,
                                              vars = vars)

  future::plan(future::multisession, workers = cores)

  fits[2:length(formulas)] = future.apply::future_lapply(2:length(formulas), FUN = function(x)
  {
    result = 0 # initialize
    # if only a 1 group LRT with only 1 mixture group in formula
    if((length(vars$mixture) == 1) && (x == (vars$mixture %>% unlist %>% length %>% sum(2))))
    {
      result = stats::logLik(stats::glm(stats::formula(paste(formulas[x], "1")), data = data, family = family))[1]
    }else # if either single pollutant LRT or 1 group LRT (more than 1 mixture group)
    {
      if(formulas[[x]] == "0 weight")
      {
        result = -1*fits[[1]]$ML_sol$value
      }else
      {
        result = -1*fit_fgwqsr_hybrid(stats::formula(formulas[x]), data, quantiles,
                                      family = family,
                                      initial_vals = initial_vals[[x-1]],
                                      optim_control_list = optim_control_list)$ML_sol$value
      }
    }

    if(verbose == TRUE)
    {
      p()
    }

    return(result)
  })

  future::plan(future::sequential)

  return(fits)
}

#' Fit a FGWQSR Model
#'
#' @param formula A formula for model fitting of FGWQSR.  Please see details for formula construction.
#' @param data  Dataframe that contains all covariates and the outcome data.
#' Column names of dataframe should match those referenced in the model formula.
#' @param quantiles Number of quantiles to quantize the exposure variables in the
#' mixture portion of the model.  Default value is 5.
#' @param family String, one of either \code{'binomial'}, \code{'gaussian'}, or \code{'poisson'} for binary,
#' continuous, or count outcomes respectively.  The link function for the binomial outcome is the logit link.
#' @param n_mvn_sims Defines resolution for simulated null distribution for group index and single
#' chemical LRTs.  Default is 10,000.
#' @param zero_threshold_cutoff Value within (0,.5] that defines how often
#' parameters estimated close to the boundary of the parameter
#' space are assigned a boundary cone in the constrained multivariate normal
#' monte carlo inference procedure.  A \code{zero_tolerance_threshold} value of 0.5 will
#' assign parameters with FGWQSR maximum likelihood estimates of precisely 0 the
#' boundary cone while a \code{zero_tolerance_threshold} value of 0 will assign all
#' parameters a boundary cone.  Reasonable values may be within [0.05, 0.5] --
#' all choices of \code{zero_tolerance_threshold} are asymptotically equivalent.
#' The default is set to zero_tolerance_threshold = 0.5.
#' @param verbose Displays messages and progress bar while fitting FGWQSR model.
#'  Default is TRUE.
#' @param cores Number of cores to parallelize on for fitting nested models and
#' simulated null LRT distributions.
#'  Default is number of available cores on user device.
#' @param optim_control_list Option to supply control options to optim.
#' @return List with attributes from fgwqsr model fitting.
#' @description
#' Fits a group signed constrained logistic regression with
#' inference for both group indices and single chemical effects on an entire
#' dataset (without requiring any data splitting).
#' @details
#' The model formula that is passed to fgwqsr() is different than traditional
#' formulas in lm and glm, as it needs to denote mixture group arrangement.
#' Three special characters are used in fgwqsr formulas:
#' \itemize{  \item \code{|} - denotes the boundary of a mixture group, used to
#' separate chemicals within a mixture group.
#'    \item  \code{/} - denotes the end of the mixture group specification,
#'    adjusting covariates can be added to the formula after
#' this character.  If no adjusting covariates, do not need to specify.
#'    \item \code{i.} - precedes categorical variables to denote a categorical variable.
#' For example, if we have a categorical variable cat_var, we would denote
#' this in the model formula by i.cat_var.
#' This is similar to the stata syntax to declare categorical variables.}
#' An example formula may look like \code{formula = y ~ x1 + x2 | x3 + x4/ height + i.city}.
#' For examples, please see the github package vignette or
#' \url{https://github.com/Daniel-Rud/fgwqsr/blob/main/README.md}
#' @seealso \code{\link{summary.fgwqsr}}
#' @export

fgwqsr = function(formula, data, quantiles = 5,
                  family = "binomial",
                  n_mvn_sims = 10000,
                  zero_threshold_cutoff = .5, verbose = T,
                  cores = future::availableCores(),
                  optim_control_list = list(maxit = 1000, factr = 1E-12, fnscale = 1))
{

  time_begin = proc.time()[3]

  # cast to numeric -- just in case!
  quantiles = as.integer(quantiles)
  n_mvn_sims = as.integer(n_mvn_sims)
  zero_threshold_cutoff = as.numeric(zero_threshold_cutoff)
  cores = as.integer(cores)

  # perform initial checks ##################################################
  # check if formula
  if(!inherits(formula,"formula"))
  {
    stop("The formula argument must be of type formula.
         If using a string formula, consider using as.formula(formula).")
  }

  # get formula and vars object that stores info about model
  f = as.character(formula); vars = clean_vars(gsub("\n", "", f[3]))

  # check if data is a dataframe
  if(!is.data.frame(data))
  {
    stop("The data argument must be a dataframe object.
         Please ensure that it is a dataframe, where the
         columnames of the dataframe correspond to the variable
         names referenced in the model formula.")
  }

  # check if variables in formula are in dataframe
  all_vars = c(vars$mixture %>% unlist, vars$continuous, vars$categorical)

  if( sum(all_vars %in% colnames(data)) != length(all_vars))
  {
    missing_vars = all_vars[which(!(all_vars %in% colnames(data)))]
    stop(paste0("The following variable names included in the model formula are not included
         as columnames in the passed `data` argument: ", paste(missing_vars, collapse = ", "), "."))
  }

  # keep only relevant covariates in dataset -- this call is so we can check for
  # complete cases
  data = data[c(f[2], vars %>% unlist)]

  # check to see if all observations have complete cases
  if(sum(stats::complete.cases(data)) != nrow(data))
  {
    warning("Dataframe contains observations with missing values. Please consider
    using na.omit(data).  Only using complete observations from passed data argument.")

    data = stats::na.omit(data)

    if(nrow(data) == 0)
    {
      stop("All observations are incomplete.  Consider using imputation, fgwqsr()
           does not allow for missing data.")
    }
  }

  # check family argument
  if(!(family %in% c("binomial", "gaussian", "poisson")))
  {
    stop("The `family` argument must be one of 'binomial', 'poisson', or 'gaussian'.")
  }



  # check if the outcome has both cases and controls -- make sure percentage
  # of cases is at least greater than 1%
  if(family == "binomial")
  {
    # check that outcome variable is a 0 1 numeric variable
    if(!is.numeric(data[f[2]] %>% unlist))
      # check that the vector is numeric
    {
      stop("The outcome variable must be coded as a numeric vector with
         0 denoting controls and 1 denoting cases")
    }else if((sum(unique(data[f[2]] %>% unlist) %in% c(0,1)) != 2) && (family == "binomial"))
      # if the vector is numeric, make sure its elements are only 0s or 1s
    {
      stop("The outcome variable must be coded as a numeric vector with
         0 denoting controls and 1 denoting cases")
    }

    if(sum(data[f[2]]) == 0) # if only 0s in outcome vector (controls)
    {
      stop("Outcome variable only contains controls, please check outcome variable coding.")

    } else if (sum(data[f[2]]) == nrow(data))# if only 1s in outcome vector (cases)
    {
      stop("Outcome variable only contains cases, please check outcome variable coding.")

    } else if(sum(data[f[2]]) / nrow(data) < .05)
    {
      message("Alert: Less than 5% of observations are cases.")
    }
  }

  if(family == "poisson")
  {
    min_y = min(data[f[2]])

    if(min_y < 0)
    {
      stop("Count outcome for `poisson` family must be nonnegative!  Negative counts are present in the outcome.")
    }

    # cast outcome to integer to ensure counts
    if(sapply(unlist(data[f[2]]), FUN = function(x) !is.integer(x)) %>% sum !=0)
    {
      message("Casting outcome variable to nonnegative integer type (counts).")
      data[f[2]] = data[f[2]] %>% unlist %>% as.integer
    }

  }


  # check if quantiles variable is within reasonable range

  if(quantiles <= 1)
  {
    message("Quantiles variable should be greater than 1.  Resetting to default value of 5.")
    quantiles = 5
  }else if(quantiles >20)
  {
    message("Quantiles argument should be less than 20.  Resetting to default value of 5.")
    quantiles = 5

  }else if(quantiles > 10)
  {
    message("Consider setting quantiles argument to no larger than deciles (quantiles = 10)")
  }

  # check if n_mvn_sims >= 100
  if(n_mvn_sims < 100)
  {
    message("The n_mvn_sims argument must be greater than 100.  Resetting to default value of 10,000")
    n_mvn_sims = 10000
  }

  # check if zero_threshold_cutoff is in (0,.5]
  if(zero_threshold_cutoff <= 0 || zero_threshold_cutoff > .5)
  {
    message("The zero_threshold_cutoff argument must be within (0,0.5].  Resetting to the default value of 0.5.")
    zero_threshold_cutoff = .5
  }

  # check that verbose is of type boolean
  if(!is.logical(verbose))
  {
    message("The verbose argument must be of type logical (either T or F).  Resetting to the default value of T.")
    verbose = T
  }

  # check that cores > 0, if cores > available cores,
  if(cores <= 0 || cores > future::availableCores())
  {
    message(paste0("The cores argument must be between 1 and ", future::availableCores(),
                  ".  Resetting to the default value which is the result of the availableCores() call, ", future::availableCores(), "."))
    cores = future::availableCores()
  }

  # check arguments in optim control list

  if(!setequal(optim_control_list, list(maxit = 1000, factr = 1E-12, fnscale = 1)))
  {
    optim_control_opts = c("trace", "fnscale", "parscale", "ndeps", "maxit",
                           "abstol", "reltol", "alpha", "beta", "gamma", "REPORT",
                           "warn.1d.NelderMead", "type", "lmm", "factr", "pgtol",
                           "temp", "tmax") # all optim control parameters

    # if there are names in the optim list that are not optim control parameters
    if( sum(names(optim_control_list) %in% optim_control_opts) != length(optim_control_list))
    {
      failed_names = names(optim_control_list)[which(!(names(optim_control_list) %in% optim_control_opts))]

      message("The following optim control parameters passed in optim_control_list are not optim control parameters: ",
              paste0(failed_names, collapse = ", "), ".  Resetting to the default optim control options list")

      optim_control_list = list(maxit = 1000, factr = 1E-12, fnscale = 1)
    }
  }

  ###########################################################################

  # get list of formulas to run submodels for LRT
  formulas = c(c(paste(format(formula), collapse = ""),get_formulas(vars,formula)))

  # run models in parallel, verbose option
  if(verbose == TRUE)
  {
    progressr::handlers("progress")
    message("Fitting full and nested FGWQSR models...")
    progressr::with_progress(fits <- fgwqsr_caller(formulas, data, quantiles,family,vars, verbose,
                                        cores, optim_control_list))
    message("\nGenerating LRT distributions under H0...")
  }else
  {
    fits <- fgwqsr_caller(formulas, data, quantiles,family, vars, verbose, cores, optim_control_list)
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

  # compute fisher information based on proper family
  observed_fisher = NULL
  if(family == "binomial")
  {
    observed_fisher = -1*logistic_hessian(B = params_logistic_form,
                                          design_matrix = design_matrix)
  }else if(family == "poisson")
  {
    observed_fisher = -1*poisson_hessian(B = params_logistic_form,
                                         design_matrix = design_matrix)
  }else # if gaussian
  {
    observed_fisher = -1*gaussian_hessian(B = params_logistic_form,
                                          design_matrix = design_matrix,
                                          sigma_2 = fgwqsr_fit$residual_var)
  }
  # compute inverse of observed fisher information
  cov_mat = pracma::pinv(observed_fisher)

  # get ML sol in group index and weight formulation
  names(params_logistic_form) = names(design_matrix)
  reparam_reg = reparam_GI_weights(params_logistic_form, vars)

  # now we perform inference

  # call function to perform inference #####################################

  inference_frames = NULL # initialize

  if(verbose == T)
  {
    progressr::with_progress(inference_frames <- perform_inference(ll_models = ll_models,
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

  # if gaussian, residual variance is estimated
  aic = ifelse(family != "gaussian",
               -2*full_ll_val + 2*length(params_logistic_form),
               -2*full_ll_val + 2*(length(params_logistic_form) + 1))

  bic = ifelse(family != "gaussian",
               -2*full_ll_val + log(n)*length(params_logistic_form),
               -2*full_ll_val + log(n)*(length(params_logistic_form)+1))

  # when certain estimates are on the boundary, convergance code issues
  # message regarding LNSRCH == line search method.  Therefore, do not look at
  # convergence code in some instances

  L_BFGS_B_convergance = list(counts = fgwqsr_fit$ML_sol$counts,
                              convergence = fgwqsr_fit$ML_sol$convergence,
                              message = fgwqsr_fit$ML_sol$message)

  return_list = list(inference_frames = inference_frames,
                     total_time = total_time,
                     n = n,
                     ll = full_ll_val,
                     formula = formula,
                     family = family,
                     aic = aic,
                     bic = bic,
                     vars = vars,
                     n_mvn_sims = n_mvn_sims,
                     cores = cores,
                     L_BFGS_B_convergance =L_BFGS_B_convergance,
                     param_cov_mat = cov_mat,
                     data = data,
                     quantiles = quantiles,
                     zero_threshold_cutoff = zero_threshold_cutoff,
                     optim_control_list = optim_control_list)
  if(family == "gaussian")
  {
    return_list$residual_var = fgwqsr_fit$residual_var
  }

  class(return_list) = "fgwqsr"

  return(return_list)
}



#' Summarize a fgwqsr model fit
#' @param object a fitted object from a fgwqsr() call
#' @param ... optional options -- see details
#' @description
#'  Summarizes a fgwqsr object.
#'@details
#'  Optional option `digits` to set the number of rounding digits to display in summary tables.
#' @export
summary.fgwqsr = function(object, ...)
{

  # check for fgwqsr object
  if(!methods::is(object, "fgwqsr"))
  {
    stop("Must pass the object output of the fgwqsr() function.  For example...
         some_model = fgwqsr(formula, data)
         summary(some_model)")
  }

  print.fgwqsr(object, ...)

}

print.fgwqsr = function(object,...)
{


  passed_args = list(...)

  digits = 6 # default
  if(length(passed_args) == 1)
  {
    digits = passed_args[[1]]
  }

  cat("\nCall: \nFGWQSR with formula '",gsub("  ", "",paste(format(object$formula), collapse = "")),"' on n = ",object$n,
      " observations and family = '", object$family, "'.", sep = "" )
  cat("\n\n", object$n_mvn_sims, " samples used for simulated LRT distirbution.",sep = "")
  cat("\n\nLog Likelihood:", object$ll, "| AIC:",object$aic, "| BIC:", object$bic)
  cat("\n\nEstimates and Inference for Group Index Effects\n\n", sep = "")
  stats::printCoefmat(object$inference_frames$group_index_frame[,-4], digits = digits,
                      signif.stars = T, signif.legend = F,
                      cs.ind = 1,
                      tst.ind = 2,
                      P.values = T,
                      has.Pvalue = T)
  cat("\nEstimates and Inference for Weights\n\n")

  current_index = 1
  for(i in 1: length(object$vars$mixture))
  {
    group_size = object$vars$mixture[[i]] %>% length
    output = object$inference_frames$weight_frame[current_index: (current_index + group_size - 1), ]
    stats::printCoefmat(output[,-4], digits = digits,
                        signif.stars = T, signif.legend = F,
                        cs.ind = 1,
                        tst.ind = 2,
                        P.values = T,
                        has.Pvalue = T)
    current_index = current_index + group_size
    cat("-------------------------------------------------\n")
  }
  cat("\nEstimates and Inference for Intercept and Adjusting Covariates\n\n")
  stats::printCoefmat(object$inference_frames$adj_param_frame[,1:4], digits = digits,
                      signif.stars = T, signif.legend = F,
                      cs.ind = 1:2,
                      tst.ind = 2,
                      P.values = T,
                      has.Pvalue = T)
  cat("\nSignificance Codes: <0.001 '***' <0.01 '**' <0.05 '*' <0.10 '.' \n")

  if(object$family == "gaussian")
  {
    cat(paste("\n(Dispersion parameter for gaussian family taken to be ", round(object$residual_var,digits), ")\n", sep = ""))
  }
  cat("\nTotal runtime for FGWQSR: ",
      ifelse(object$total_time < 60,paste(round(object$total_time,2), "seconds"),
             paste(round(object$total_time/60, 2), "minutes")), "on",object$cores, "cores." )

}

# Nonparametric bootstrap doesnt work in this situation -- appeals to central limit theorem,
# requires some smoothness of the functional statistic fo the bootstrap

#' #' Generate confidence intervals for group and single pollutant effects through nonparametric bootstrapping.
#' #' @param object a fitted object from a fgwqsr() call
#' #' @param parm should be left unspecified
#' #' @param level level of confidence interval to produce -- default is 0.95 for 95\% bootstrap confidence intervals
#' #' @param ... optional options -- see details
#' #' @description
#' #'  Provides confidence intervals for group index, single chemical, and adjusting covariate effects through
#' #'  nonparameteric bootstrapping.  Bootstrap samples are stratified by cases and controls, with each bootstrap sample keeping
#' #'  same proportions of cases and controls as in the original data.  Each bootstrap sample has the same number of observations
#' #'  as the original dataset.  NOTE: Bootstrapped samples are sampled from the original data with replacement!
#' #' @details
#' #' One can provide the following two arguments:
#' #'`boot_reps` - number of bootstrapped replicates
#' #'`verbose` - allows messages and progress bars to track progress of bootstrapping.
#' #'
#' #'
#' #' @export
#' #'
#' confint.fgwqsr = function(object, parm, level = 0.95, ...)
#' {
#'   passed_args = list(...)
#'
#'   if(!methods::is(object, "fgwqsr"))
#'   {
#'     stop("Must pass the object output of the fgwqsr() function.  For example...
#'          some_model = fgwqsr(formula, data)
#'          summary(some_model)")
#'   }
#'
#'   boot_reps = 1000 # initialize
#'   if("boot_reps" %in% names(passed_args))
#'   {
#'     if(!is.numeric(passed_args$boot_reps))
#'     {
#'       stop("The passed argument `boot_reps` must be an integer.")
#'     }else if(passed_args$boot_reps <= 10)
#'     {
#'       stop("The passed argument `boot_reps` should be a positive integer greater than 10")
#'     }else
#'     {
#'       boot_reps = passed_args$boot_reps
#'     }
#'   }
#'
#'   verbose = T
#'   if("verbose" %in% names(passed_args))
#'   {
#'     if(!is.logical(passed_args$verbose))
#'     {
#'       stop("The passed argument `verbose` must either be a logical -- either `TRUE` or `FALSE`. ")
#'     }else
#'     {
#'       verbose = passed_args$verbose
#'     }
#'   }
#'
#'   parm = NA
#'
#'   # stratified bootstrap procedure here -- we do not use `boot` because we want progress bars!
#'   # we sample cases and controls based on their proportion int he dataset.
#'
#'   if(verbose == T)
#'   {
#'     message("Generating bootstrap replicates...")
#'   }
#'
#'   progressr::with_progress(bootstrap_results <- generate_bootstrap(object, level,boot_reps, verbose))
#'
#'
#'   return(bootstrap_results)
#'
#' }
#'
#' generate_bootstrap = function(object, level, boot_reps, verbose)
#' {
#'   # stratified bootstrap procedure here -- we do not use `boot` because we want progress bars!
#'   # we sample cases and controls based on their proportion int he dataset.
#'
#'   total_data = object$data
#'
#'   n = nrow(total_data)
#'
#'   f = as.character(object$formula)
#'
#'   f[3] = gsub("\n", "", f[3])
#'
#'   outcome = f[2]
#'
#'   candidate_indices_cases = (1:n)[object$data[, outcome] ==1]
#'   candidate_indices_controls = (1:n)[object$data[, outcome] ==0]
#'
#'   percent_cases = sum(object$data[, outcome]) / n
#'   percent_controls = 1-percent_cases
#'
#'   if(verbose == T)
#'   {
#'     p <- progressr::progressor(steps = boot_reps)
#'   }
#'
#'
#'   future::plan(future::multisession, workers = object$cores)
#'
#'   boot_results = future.apply::future_sapply(1:boot_reps, FUN = function(x)
#'   {
#'     case_indices_sample = sample(x = candidate_indices_cases,
#'                                  size = round(n * percent_cases),
#'                                  replace = T)
#'
#'     control_indices_sample = sample(x = candidate_indices_controls,
#'                                     size = round(n * percent_controls),
#'                                     replace = T)
#'
#'     boot_data = total_data[c(case_indices_sample,control_indices_sample), ]
#'
#'
#'     initial_cov_vals = get_cov_initial_vals(formula = object$formula,
#'                                             data = boot_data,
#'                                             quantiles = object$quantiles)
#'     # we fit on one core -- parallelize through bootstrap iterations
#'     fg_fit = fit_fgwqsr_MO(formula = object$formula,
#'                            data = boot_data,
#'                            quantiles = object$quantiles,
#'                            return_y = F,
#'                            return_data = F,
#'                            initial_cov_vals = initial_cov_vals,
#'                            optim_control_list = object$optim_control_list,
#'                            cores = 1)
#'
#'     sol_logistic = fg_fit$ML_sol$par
#'     names(sol_logistic) = c("Intercept", names(boot_data)[-1])
#'
#'     # get group effects
#'     group_effects = vector(mode = "numeric", length = length(vars$mixture))
#'
#'     run_ind = 2 # first element is intercept
#'     for(i in 1: length(vars$mixture))
#'     {
#'       group_effects[i] = sum(sol_logistic[run_ind: (run_ind + length(vars$mixture[[i]]) -1)])
#'       run_ind = run_ind + length(vars$mixture[[i]])
#'     }
#'
#'     names(group_effects) = paste("Mixture Effect", 1:length(vars$mixture))
#'
#'     if(verbose == T)
#'     {
#'       p()
#'     }
#'
#'     return(c(group_effects,sol_logistic))
#'   }, future.seed = 2023) %>% t
#'
#'   future::plan(future::sequential)
#'
#'   result_frame = matrix(nrow = ncol(boot_results), ncol = 4)
#'
#'   for(i in 1:ncol(boot_results))
#'   {
#'     result_frame[i,2] = boot_results[,i] %>% mean
#'     result_frame[i,3:4] = stats::quantile(boot_results[,i],
#'                                           probs = c((1-level)/2, 1-(1-level)/2))
#'   }
#'
#'   result_frame[,1] = c(object$inference_frames$group_index_frame$Estimate,
#'                        object$inference_frames$adj_param_frame$Estimate[1],
#'                        rep(object$inference_frames$group_index_frame$Estimate,
#'                            times = sapply(object$vars$mixture, length)) *
#'                          object$inference_frames$weight_frame$`Weight Estimate`,
#'                        object$inference_frames$adj_param_frame$Estimate[-1])
#'
#'
#'   rownames(result_frame) = colnames(boot_results)
#'   colnames(result_frame) = c("FGWQSR Estimate","Bootstrap Mean", paste0((1-level)/2, " %"),
#'                              paste0( 1-(1-level)/2, " %"))
#'
#'
#'   # bootstrap_results = boot::boot(data = object$data,
#'   #                                statistic = function(data,i)
#'   #                                  {
#'   #                                  boot_data = data[i,]
#'   #
#'   #                                  f = as.character(formula)
#'   #                                  outcome = f[2]
#'   #                                  vars = clean_vars(gsub("\n", "", f[3]))
#'   #
#'   #                                  initial_cov_vals = get_cov_initial_vals(formula = object$formula,
#'   #                                                                          data = boot_data,
#'   #                                                                          quantiles = object$quantiles)
#'   #                                  # we fit on one core -- parallelize through bootstrap iterations
#'   #                                  fg_fit = fit_fgwqsr_MO(formula = object$formula,
#'   #                                                         data = boot_data,
#'   #                                                         quantiles = object$quantiles,
#'   #                                                         return_y = F,
#'   #                                                         return_data = F,
#'   #                                                         initial_cov_vals = initial_cov_vals,
#'   #                                                         optim_control_list = object$optim_control_list,
#'   #                                                         cores = 1)
#'   #
#'   #                                  sol_logistic = fg_fit$ML_sol$par
#'   #                                  names(sol_logistic) = c("Intercept", names(boot_data)[-1])
#'   #
#'   #                                  # get group effects
#'   #                                  group_effects = vector(mode = "numeric", length = length(vars$mixture))
#'   #
#'   #                                  run_ind = 2 # first element is intercept
#'   #                                  for(i in 1: length(vars$mixture))
#'   #                                  {
#'   #                                    group_effects[i] = sum(sol_logistic[run_ind: (run_ind + length(vars$mixture[[i]]) -1)])
#'   #                                    run_ind = run_ind + length(vars$mixture[[i]])
#'   #                                  }
#'   #
#'   #                                  names(group_effects) = paste("Mixture Effect", 1:length(vars$mixture))
#'   #
#'   #                                  return(c(group_effects,sol_logistic ))
#'   #                                },
#'   #                                strata = object$data[, outcome],
#'   #                                R = 10, )
#'
#'   # result_frame = matrix(nrow = ncol(bootstrap_results$t), ncol = 3)
#'   # for(i in 1:ncol(bootstrap_results$t))
#'   # {
#'   #   result_frame[i,1] = bootstrap_results$t0[i]
#'   #   result_frame[i,2] = bootstrap_results$t[,i] %>% mean
#'   #   result_frame[i,3:4] = boot::boot.ci(bootstrap_results,type = "perc",
#'   #                                       index = i, conf = level)$percent[1,4:5]
#'   # }
#'   #
#'   # rownames(result_frame) = names(bootstrap_results$t0)
#'   # colnames(result_frame) = c("FGWQSR Estimate","Bootstrap Mean", paste0((1-level)/2, " %"),
#'   #                            paste0( 1-(1-level)/2, " %"))
#'
#'   return(list(result_frame = result_frame,
#'               bootstrap_results = boot_results))
#' }

