

# 12/13/2021
# Author: Daniel Rud
# We implement the BGWQSR in the form of Wheeler et al (2021).
# the formula should be in the form of y ~ u + v | w + x / y + i.z
# where '|' denotes the group separations and  '/' denotes
# the separation between mixture components and confounders
# i. denotes creating a categorical response


clean_vars = function(vars) # separates formula into mixture component, continuous confounder,
  # and categorical confounder names
{
  sep = strsplit(vars, split = "/", fixed = TRUE)[[1]]

  mixtureString = sep[1]

  mixtureNames = strsplit(strsplit(gsub(" ", "", mixtureString), "\\|")[[1]], "\\+")

  confounders = strsplit(gsub("\\+", " ", sep[2]), " ")[[1]]

  confounders = confounders[which(confounders!="")]

  catInd = which(substr(confounders, 1, 2) == "i.")


  if(length(catInd)>0)
  {
    catConfounders = gsub("i.", "", confounders[catInd], fixed = TRUE)

    contConfounders = confounders[-catInd]

    variables = list(mixture = mixtureNames, continuous = contConfounders, categorical = catConfounders)

  }else
  {
    if(length(confounders) > 0)
    {

      contConfounders = confounders

      variables = list(mixture = mixtureNames, continuous = contConfounders, categorical = NULL)

    }else
    {
      variables = list(mixture = mixtureNames, continuous = NULL, categorical = NULL)
    }

  }

  return(variables)
}

quantize_vars = function(data, vars, quantiles)
{

  mixtureComps = unlist(vars$mixture)

  for(i in 1: length(mixtureComps))
  {
    data[[mixtureComps[i]]] = statar::xtile(data[[mixtureComps[i]]], n = quantiles) - 1
  }

  return(data)
}

create_jags_model_string = function(data,y, vars)
{
  numBetas = length(vars$mixture)

  betaNames = paste("B",seq(0,numBetas,1), sep = "")

  weights = list()

  for(i in 1: numBetas) # weight arrangement
  {
    numWeights = length(vars$mixture[[i]])
    weights[[i]] = paste("w",i,"[", seq(1:numWeights),"]",sep = "")
  }

  phi1 = c()

  length1 = length(vars$continuous)

  length2 = length(vars$categorical)

  if(length1>0)
  {
    for(i in 1: length1)
    {
      phi1[i] = paste("phi_", vars$continuous[i], sep = "")
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
      phi2[[i]] = paste("phi_", names, sep= "")

    }
  }


  mixString = "B0 "

  for(i in 1: length(vars$mixture))
  {
    weightString = paste("(",
                         paste(paste(weights[[i]], "*", vars$mixture[[i]],"[i]", sep = ""), collapse = " + "),
                         ")", sep = "")
    mixString = paste(mixString," + ", betaNames[i+1], "*", weightString, sep = "" )

  }

  phi1String = ""

  if(length1>0)
    phi1String = paste(" + ", paste(paste(phi1, "*", vars$continuous,"[i]", sep = ""), collapse = " + "), sep = "")

  phi2String = ""

  if(length2>0)
  {

    for(i in 1: length2)
    {
      catName = paste(vars$categorical[i], "_", sep = "")
      varNames = names(data)[which(substr(names(data), 1, nchar(catName)) == catName)]
      phi2String = paste(paste(phi2String," + ", paste(paste(phi2[[i]], "*", varNames,"[i]", sep = ""), collapse = " + ") ))
    }
  }

  ystring = paste(y,"[i] ~ dbern(pi[i]) \n logit(pi[i]) = ", mixString, phi1String, phi2String, sep = "")

  weightPriors = ""

  for(i in 1:length(weights))
  {
    numWeights = length(weights[[i]])
    alpha = toString(rep(1/numWeights, numWeights))
    weightPriors = paste(weightPriors, " \n w",i, "[1:", numWeights, "] ~ ddirch( c(", alpha, ") )", sep = "")
  }



  betaPriors = ""

  sigmaBetas = paste("sigma_", betaNames, sep = "")

  for(i in 1:length(betaNames))
  {
    betaPriors = paste(betaPriors, " \n ", betaNames[i], " ~ dnorm(0, 1/(", sigmaBetas[i],"^2))", sep = "")
  }

  phis = unlist(c(phi1, phi2))

  sigmaPhis = ""

  if(length(phis)>0)
  {
    sigmaPhis = paste("sigma_", substr(phis, 5, sapply(phis, nchar)), sep = "")
  }

  phiPriors = ""

  if(length(phis)>0)
  {
    for(i in 1: length(phis))
    {
      phiPriors = paste(phiPriors, " \n ", phis[i], " ~ dnorm(0, (", sigmaPhis[i],"^2))", sep = "")
    }
  }

  sigmas = c()

  if(length(phis)>0)
  {
    sigmas =  c(sigmaBetas, sigmaPhis)
  }else
  {
    sigmas = sigmaBetas
  }

  sigmaPriors = ""

  for(i in 1:length(sigmas))
  {
    sigmaPriors = paste(sigmaPriors, " \n ", sigmas[i], " ~ dunif(0,100)", sep = "")
  }


  mod_string = paste("model{ \n for(i in 1: length(", y, ")) \n { \n",
                     ystring, " \n } \n ",
                     weightPriors,
                     betaPriors,
                     phiPriors,
                     sigmaPriors,
                     "}", sep = "")

  paramlist = c()

  if(length(phis)>0)
  {
    paramList = list(modelString = mod_string, betaNames = betaNames, phis = phis,
                     weights = weights, sigmaBetas=sigmaBetas, sigmaPhis = sigmaPhis)
  }else{
    paramList = list(modelString = mod_string, betaNames = betaNames,
                     weights = weights, sigmaBetas=sigmaBetas)
  }

  return(paramList)
}


#' Fit a Bayesian Grouped Weighted Quantile Sum Regression Model
#' @description
#' Fits a BGWQSR using the runjags package in parallel over multiple cores
#' @param `formula`  A formula for model fitting of BGWQSR.  Please see description for formula construction
#' @param `data` dataframe that contains all covariates and the outcome data.  Column names of dataframe should match those referenced int he model formula.
#' @param `quantiles` number of quantiles to quantize the exposure variables in the mixture portion of the model.  Default value is 5.
#' @param `n.iter` number of mcmc iterations after burnin and adapt iterations PER mcmc chain.
#' @param `n.burnin` number of mcmc burnin samples PER mcmc chain
#' @param `n.thin` thinning interval for mcmc samples PER mcmc chain
#' @param `n.chains` number of separate independent mcmc chains to use.
#' @param `n.adapt` number of mcmc samples to perform mcmc adaption PER mcmc chain.
#' @param `inits` initial values to provide for prior distributions.
#' @param `method` method for mcmc fitting, a passed argument to run.jags function.  Can be one of: `rjags`, `simple`, `interruptible`,
#' `parallel`, `rjparallel`, `background`, `bgparallel`, or `snow`.
#' @return list with attributes from fgwqsr model fitting.
#' @export

bgwqsr= function(formula, data, quantiles = 5,  n.iter = 10000 / n.chains, n.burnin = 5000,
                 n.thin = 1, n.chains=3, n.adapt= 1000, inits = NA, method = "parallel")
{
  if(!inherits(formula,"formula"))
  {
    stop("The formula argument must be of type formula.
         If using a string formula, consider using as.formula(formula).")
  }

  f = as.character(formula)

  f[3] = gsub("\n", "", f[3])

  y = f[2]

  vars = clean_vars(f[3])

  # performing initial checks ##################################################

  # check if data is a dataframe
  if(!is.data.frame(data))
  {
    stop("The data argument must be a dataframe object.
         Please ensure that it is a dataframe, where the
         columnames of the dataframe correspond to the variable
         names referenced in the model formula.")
  }

  # check that outcome variable is a 0 1 numeric variable
  if(!is.numeric(data[f[2]] %>% unlist))
    # check that the vector is numeric
  {
    stop("The outcome variable must be coded as a numeric vector with
         0 denoting controls and 1 denoting cases")
  }else if(sum(unique(data[f[2]] %>% unlist) %in% c(0,1)) != 2)
    # if the vector is numeric, make sure its elements are only 0s or 1s
  {
    stop("The outcome variable must be coded as a numeric vector with
         0 denoting controls and 1 denoting cases")
  }

  # check if the outcome has both cases and controls -- make sure percentage
  # of cases is at least greater than 1%

  if(sum(data[f[2]]) == 0) # if only 0s in outcome vector (controls)
  {
    stop("Outcome variable only contains controls, please check outcome variable coding.")

  } else if (sum(data[f[2]]) == nrow(data))# if only 1s in outcome vector (cases)
  {
    stop("Outcome variable only contains cases, please check outcome variable coding.")

  } else if(sum(data[f[2]]) / nrow(data) < .05)
  {
    warning("Case control ratio is extremely unproportional.  Less than 5% of observations are cases.  Proceed with caution.")
  }

  # check if variables in formula are in dataframe
  all_vars = c(vars$mixture %>% unlist, vars$continuous, vars$categorical)

  if( sum(all_vars %in% colnames(data)) != length(all_vars))
  {
    missing_vars = all_vars[which(!(all_vars %in% colnames(data)))]
    stop(paste0("The following variable names included in the model formula are not included
         as columnames in the passed `data` argument: ", paste(all_vars, collapse = ", "), "."))
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

  # check n.chains > 0
  if(n.chains <=0)
  {
    message("The `n.chains` argument should be greater than 0.
            Setting n.chains = 1")
    n.chains = 1
  }

  # check that n.iter is greater than 100
  if(n.iter < 100)
  {
    message("The `n.iter` argument should be greater than 100 -- we suggest >= 1,000
            and ideally 10,000 (in total, n.iter is the number of iterations per mcmc chain).
            Resetting to default value of 10,000 / n.chain")
    n.iter = 10000 / n.chains
  }

  # check that n.burnin is greater than 0
  if(n.burnin <=0)
  {
    message("The `n.burnin` argument should be greater than 0 -- we suggest >= 1,000.
            Resetting to default value of 5,000")
    n.burnin = 5000
  }

  # check that n.thin > 0
  if( n.thin <=0)
  {
    message("The `n.thin` argument should be greater than 0.
            Resetting to default value of 1")
    n.thin = 1
  }

  if(n.adapt < 0)
  {
    message("The `n.adapt` argument should be greater than or equal to 0.
            Setting value to 0.  ")
    n.adapt = 0
  }

  # check is method is a run.jags method option
  if(method != "parallel")
  {
    runjags_methods = c("rjags", "simple", "interruptible",
                        "parallel", "rjparallel", "background",
                        "bgparallel", "snow")

    if(!(method %in% runjags_methods))
    {
      message(paste0("The argument `method` should be one of the following options: ", paste(runjags_methods, collapse = ", "),
                     ". Resetting to default value of 'parallel'.") )

      method = "parallel"
    }
  }

  ##############################################################################

  if(length(vars$categorical)>0) # if there are categorical variables
  {
    data = fastDummies::dummy_cols(data, select_columns = vars$categorical,
                                   remove_first_dummy = TRUE, remove_selected_columns = TRUE) # add dummy variables
  }

  data = quantize_vars(data,vars, quantiles) # quantize mixture components

  model = create_jags_model_string(data,y, vars) # create form of Jags model

  model_string = model[[1]]

  save_params = unlist(c(model$betaNames, model$phis, model$weights))

  suppressWarnings(jags.model <- runjags::run.jags(model = model_string, data = data,
                                                   monitor = save_params,
                                                   n.chains = n.chains, sample = n.iter, burnin = n.burnin,
                                                   adapt = n.adapt, thin = n.thin, method = method,
                                                   inits = inits, modules = "glm"))

  jags.model = runjags::add.summary(jags.model)

  model_params = unlist(model[-1])

  return_list = list(model = jags.model,vars = vars, params =  model_params)

  return(return_list)
}


#' Plot Results of Beta and Weight Estimates
#' @description
#' A short description...
#'
#' @param model a BGWQSR model object fit from bgwqsr_parallel
#' @param filename optional file path to save results of plot in a pdf format
#' @param weight_axis_pos choice to put legend for weight plots on the left or right side of plot.  Options are "left" or "right".  Default is "left"
#' @param beta_axis_pos choice to put legend for group index plots on the left or right side of plot. Options are "left" or "right".  Default is "left"
#' @return void, prints plots, optionally saves plot in pdf format to filename
#' @import ggplot2
#' @export
plot_results = function(model, filename = NULL, weight_axis_pos = "left", beta_axis_pos = "left")
{

  # initial checks #############################################################

  # check for bgwqsr object
  if(!setequal(names(model), c("model", "vars", "params")))
  {
    stop("Must pass the object output of the bgwqsr() function.  For example...
         some_model = bgwqsr(formula, data)
         plot_results(some_model)")
  }

  if(!(weight_axis_pos %in% c("left", "right")))
  {
    message("The argument `weight_axis_pos` can only be set to `left` or `right`.
            Resetting to default value of `left`. ")

    weight_axis_pos = "left"
  }

  if(!(beta_axis_pos %in% c("left", "right")))
  {
    message("The argument `beta_axis_pos` can only be set to `left` or `right`.
            Resetting to default value of `left`. ")

    beta_axis_pos = "left"
  }

  ##############################################################################

  results = model$model$summaries

  #Weight Plotting

  numMixes = length(model$vars$mixture)

  group = c() # for group category

  for(i in 1: numMixes) # create group vector
  {
    group = c(group, rep(i, length(model$vars$mixture[[i]])))
  }

  weightNames = unlist(model$vars$mixture)

  weightData = results[which(substr(row.names(results), 1,1)=="w"), c(1,3,4)]

  names = row.names(weightData)

  for(i in 1: numMixes) # order weights for when >= 10 weights in a group
  {
    searchW = paste("w", i, sep = "")

    locs = which(substr(names, 1, nchar(searchW))==searchW)

    minLoc = min(locs)

    maxLoc = max(locs)

    tempData = weightData[minLoc:maxLoc, ]

    weightNumbers = as.numeric(substr(names[minLoc:maxLoc],
                                      regexpr("[", names[minLoc:maxLoc], fixed = TRUE)+1, nchar(names[minLoc:maxLoc])-1))

    newOrder = order(weightNumbers)

    tempData = tempData[newOrder, ]

    weightData[minLoc:maxLoc, ] = tempData

    row.names(weightData)[minLoc:maxLoc] = row.names(tempData)

  }

  weightPlot = data.frame(name = factor(weightNames, levels = weightNames), group = factor(group), mean = weightData[, 3],
                          L2.5 = weightData[, 1], U97.5 = weightData[, 2]  )


  plotWeights = ggplot(data = weightPlot, mapping = aes(x = .data$mean, y = .data$name)) +
    geom_point(mapping = aes(color = .data$group), size = 3) +
    geom_errorbar(aes(xmin =.data$L2.5, xmax = .data$U97.5), linewidth = .8) + xlab("2.5% quantile, Mean, 97.5% quantile")+
    ylab("Pollutant Name") + scale_y_discrete(position = weight_axis_pos) +
    theme(panel.grid.major.y = element_line(colour = "grey50", linewidth = 0.2))

  #scale_color_manual(values = c("1" = "deepskyblue1", "2" = "red")) +

  #Beta / Phi Plotting

  betas = results[which(substr(row.names(results), 1,1)=="B"), c(1,3,4)]
  phis = results[which(substr(row.names(results), 1,1)=="p"), c(1,3,4)]

  effectEstimates = rbind(betas, phis)

  significant = ifelse(effectEstimates[,1]>0 | effectEstimates[,3]<0, "sig", "insig")

  effectData = data.frame(name = factor(row.names(effectEstimates)), mean = effectEstimates[, 3],
                          L2.5 = effectEstimates[,1], U97.5 = effectEstimates[, 2],
                          significant = factor(significant))

  plotBetas = ggplot(data = effectData, mapping = aes(x = .data$mean, y = .data$name)) +
    geom_vline(xintercept = 0, linetype = "longdash", colour = "red") +
    geom_point(size = 3, aes(color = .data$significant)) +
    scale_color_manual(values = c("sig" = "red", "insig" = "black")) +
    geom_errorbar(aes(xmin =.data$L2.5, xmax = .data$U97.5), linewidth = .8) + ylab("Effect Name") +
    xlab("2.5% quantile, Mean, 97.5% quantile") +
    scale_y_discrete(position = beta_axis_pos) +
    theme(panel.grid.major.y = element_line(colour = "grey50", linewidth = 0.2))

  gridExtra::grid.arrange(plotWeights, plotBetas, ncol = 2)

  if(!is.null(filename))
  {
    grDevices::pdf(paste(filename, ".pdf", sep = ""))
    gridExtra::grid.arrange(plotWeights, plotBetas, ncol = 1)
    grDevices::dev.off()
  }
}

#' Plot Weight Results
#' @description
#' Plots the 95% credible intervals for each single chemical weight, along with posterior mean weight estimates, from the BGWQSR weight posterior distributions.
#' @param model a BGWQSR model object fit from bgwqsr_parallel
#' @param filename optional file path to save results of plot in a pdf format
#' @param weight_axis_pos choice to put legend for weight plots on the left or right side of plot.  Options are "left" or "right".  Default is "left"
#' @return void, prints weight plot, optionally saves plot in pdf format to filename
#' @import ggplot2
#' @export
plot_weights = function(model, filename = NULL, weight_axis_pos = "left")
{

  # initial checks #############################################################

  # check for bgwqsr object
  if(!setequal(names(model), c("model", "vars", "params")))
  {
    stop("Must pass the object output of the bgwqsr() function.  For example...
         some_model = bgwqsr(formula, data)
         plot_results(some_model)")
  }

  if(!(weight_axis_pos %in% c("left", "right")))
  {
    message("The argument `weight_axis_pos` can only be set to `left` or `right`.
            Resetting to default value of `left`. ")

    weight_axis_pos = "left"
  }

  ##############################################################################
  results = model$model$summaries

  #Weight Plotting

  numMixes = length(model$vars$mixture)

  group = c() # for group category

  for(i in 1: numMixes) # create group vector
  {
    group = c(group, rep(i, length(model$vars$mixture[[i]])))
  }

  weightNames = unlist(model$vars$mixture)

  weightData = results[which(substr(row.names(results), 1,1)=="w"), c(1,3,4)]

  names = row.names(weightData)

  for(i in 1: numMixes) # order weights for when >= 10 weights in a group
  {
    searchW = paste("w", i, sep = "")

    locs = which(substr(names, 1, nchar(searchW))==searchW)

    minLoc = min(locs)

    maxLoc = max(locs)

    tempData = weightData[minLoc:maxLoc, ]

    weightNumbers = as.numeric(substr(names[minLoc:maxLoc],
                                      regexpr("[", names[minLoc:maxLoc], fixed = TRUE)+1, nchar(names[minLoc:maxLoc])-1))

    newOrder = order(weightNumbers)

    tempData = tempData[newOrder, ]

    weightData[minLoc:maxLoc, ] = tempData

    row.names(weightData)[minLoc:maxLoc] = row.names(tempData)

  }

  weightPlot = data.frame(name = factor(weightNames, levels = weightNames), group = factor(group), mean = weightData[, 3],
                          L2.5 = weightData[, 1], U97.5 = weightData[, 2]  )


  plotWeights = ggplot(data = weightPlot, mapping = aes(x = .data$mean, y = .data$name)) +
    geom_point(mapping = aes(color = .data$group), size = 3) +
    geom_errorbar(aes(xmin =.data$L2.5, xmax = .data$U97.5), linewidth = .8) + xlab("2.5% quantile, Mean, 97.5% quantile")+
    ylab("Pollutant Name") + scale_y_discrete(position = weight_axis_pos) +
    theme(panel.grid.major.y = element_line(colour = "grey50", linewidth = 0.2))


  print(plotWeights)

  if(!is.null(filename))
  {
    grDevices::pdf(paste(filename, ".pdf", sep = ""))
    print(plotWeights)
    grDevices::dev.off()
  }

}

#' Plot Results of Betas
#' @description
#' Plots the posterior 95% credible interval for group index posteriors from BGWQSR, as well as posterior mean.
#' @param model a BGWQSR model object fit from bgwqsr_parallel
#' @param filename optional file path to save results of plot in a pdf format
#' @param beta_axis_pos choice to put legend for group index plots on the left or right side of plot. Options are "left" or "right".  Default is "left"
#' @return void, prints plots, optionally saves plot in pdf format to filename
#' @import ggplot2
#' @export

plot_betas = function(model, filename = NULL, beta_axis_pos = "left") # enter just name of file, not with .pdf
{

  # initial checks #############################################################

  # check for bgwqsr object
  if(!setequal(names(model), c("model", "vars", "params")))
  {
    stop("Must pass the object output of the bgwqsr() function.  For example...
         some_model = bgwqsr(formula, data)
         plot_results(some_model)")
  }

  if(!(beta_axis_pos %in% c("left", "right")))
  {
    message("The argument `beta_axis_pos` can only be set to `left` or `right`.
            Resetting to default value of `left`. ")

    beta_axis_pos = "left"
  }

  ##############################################################################

  results = model$model$summaries


  betas = results[which(substr(row.names(results), 1,1)=="B"), c(1,3,4)]
  phis = results[which(substr(row.names(results), 1,1)=="p"), c(1,3,4)]

  effectEstimates = rbind(betas, phis)

  significant = ifelse(effectEstimates[,1]>0 | effectEstimates[,3]<0, "sig", "insig")

  effectData = data.frame(name = factor(row.names(effectEstimates)), mean = effectEstimates[, 3],
                          L2.5 = effectEstimates[,1], U97.5 = effectEstimates[, 2],
                          significant = factor(significant))

  plotBetas = ggplot(data = effectData, mapping = aes(x = .data$mean, y = .data$name)) +
    geom_vline(xintercept = 0, linetype = "longdash", colour = "red") +
    geom_point(size = 3, aes(color = .data$significant)) +
    scale_color_manual(values = c("sig" = "red", "insig" = "black")) +
    geom_errorbar(aes(xmin =.data$L2.5, xmax = .data$U97.5), linewidth = .8) + ylab("Effect Name") +
    xlab("2.5% quantile, Mean, 97.5% quantile") +
    scale_y_discrete(position = beta_axis_pos) +
    theme(panel.grid.major.y = element_line(colour = "grey50", linewidth = 0.2))


  print(plotBetas)

  if(!is.null(filename))
  {
    grDevices::pdf(paste(filename, ".pdf", sep = ""))
    print(plotBetas)
    grDevices::dev.off()
  }
}


