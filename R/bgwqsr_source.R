

# 12/13/2021
# Author: Daniel Rud
# We implement the BGWQSR in the form of Wheeler et al (2021).
# the formula should be in the form of y ~ u + v | w + x / y + i.z
# where '|' denotes the group separations and  '/' denotes
# the separation between mixture components and confounders
# i. denotes creating a categorical response


cleanVars = function(vars) # separates formula into mixture component, continuous confounder,
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

quantizeVars = function(data, vars, quantiles)
{

  mixtureComps = unlist(vars$mixture)

  for(i in 1: length(mixtureComps))
  {
    data[[mixtureComps[i]]] = statar::xtile(data[[mixtureComps[i]]], n = quantiles) - 1
  }

  return(data)
}

createJagsModelString = function(data,y, vars)
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

bgwqsr_parallel = function(formula, data, quantiles = 4,  n.iter = 1000, n.burnin = 5000,
                           n.thin = 1, n.chains=3, n.adapt= 1000, inits = NA)
{
  f = as.character(formula)

  f[3] = gsub("\n", "", f[3])

  y = f[2]

  vars = cleanVars(f[3])

  if(length(vars$categorical)>0) # if there are categorical variables
  {
    data = fastDummies::dummy_cols(data, select_columns = vars$categorical,
                      remove_first_dummy = TRUE, remove_selected_columns = TRUE) # add dummy variables
  }

  data = quantizeVars(data,vars, quantiles) # quantize mixture components

  model = createJagsModelString(data,y, vars) # create form of Jags model

  modelString = model[[1]]

  save_params = unlist(c(model$betaNames, model$phis, model$weights))

  jags.model = runjags::run.jags(model = modelString, data = data,
                        monitor = save_params,
                        n.chains = n.chains, sample = n.iter, burnin = n.burnin,
                        adapt = n.adapt, thin = n.thin, method = "parallel",
                        inits = inits, modules = "glm")
  jags.model = runjags::add.summary(jags.model)

  modelParams = unlist(model[-1])

  returnlist = list(model = jags.model,vars = vars, params =  modelParams)

  return(returnlist)
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
plotResults = function(model, filename = NULL, weight_axis_pos = "left", beta_axis_pos = "left")
{
  results = results = model$results[[2]]

  if(length(names(model))==3) # if it was run in parallel
  {
    results = model$model$summary[[2]]
  }

  #Weight Plotting

  numMixes = length(model$vars$mixture)

  group = c() # for group category

  for(i in 1: numMixes) # create group vector
  {
    group = c(group, rep(i, length(model$vars$mixture[[i]])))
  }

  weightNames = unlist(model$vars$mixture)

  weightData = results[which(substr(row.names(results), 1,1)=="w"), c(1,3,5)]

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

  weightPlot = data.frame(name = factor(weightNames, levels = weightNames), group = factor(group), median = weightData[, 2],
                          L2.5 = weightData[, 1], U97.5 = weightData[, 3]  )


  plotWeights = ggplot(data = weightPlot, mapping = aes(x = .data$median, y = .data$name)) +
    geom_point(mapping = aes(color = .data$group), size = 3) +
    geom_errorbar(aes(xmin =.data$L2.5, xmax = .data$U97.5), size = .8) + xlab("2.5% quantile, Median, 97.5% quantile")+
    ylab("Pollutant Name") + scale_y_discrete(position = weight_axis_pos) +
    theme(panel.grid.major.y = element_line(colour = "grey50", size = 0.2))

  #scale_color_manual(values = c("1" = "deepskyblue1", "2" = "red")) +

  #Beta / Phi Plotting

  betas = results[which(substr(row.names(results), 1,1)=="B"), c(1,3,5)]
  phis = results[which(substr(row.names(results), 1,1)=="p"), c(1,3,5)]

  effectEstimates = rbind(betas, phis)

  significant = ifelse(effectEstimates[,1]>0 | effectEstimates[,3]<0, "sig", "insig")

  effectData = data.frame(name = factor(row.names(effectEstimates)), median = effectEstimates[, 2],
                          L2.5 = effectEstimates[,1], U97.5 = effectEstimates[, 3],
                          significant = factor(significant))

  plotBetas = ggplot(data = effectData, mapping = aes(x = .data$median, y = .data$name)) +
    geom_vline(xintercept = 0, linetype = "longdash", colour = "red") +
    geom_point(size = 3, aes(color = .data$significant)) +
    scale_color_manual(values = c("sig" = "red", "insig" = "black")) +
    geom_errorbar(aes(xmin =.data$L2.5, xmax = .data$U97.5), size = .8) + ylab("Effect Name") +
    xlab("2.5% quantile, Median, 97.5% quantile") +
    scale_y_discrete(position = beta_axis_pos) +
    theme(panel.grid.major.y = element_line(colour = "grey50", size = 0.2))

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
plotWeights = function(model, filename = NULL, weight_axis_pos = "left")
{
  numMixes = length(model$vars$mixture)

  results = model$results[[2]]

  if(length(names(model))==3) # if it was run in parallel
  {
    results = model$model$summary[[2]]
  }

  group = c() # for group category

  for(i in 1: numMixes) # create group vector
  {
    group = c(group, rep(i, length(model$vars$mixture[[i]])))
  }

  weightNames = unlist(model$vars$mixture)

  weightData = results[which(substr(row.names(results), 1,1)=="w"), c(1,3,5)]

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

  weightPlot = data.frame(name = factor(weightNames, levels = weightNames), group = factor(group), median = weightData[, 2],
                          L2.5 = weightData[, 1], U97.5 = weightData[, 3] )


  plotWeights = ggplot(data = weightPlot, mapping = aes(x = .data$median, y = .data$name)) +
    geom_point(mapping = aes(color = .data$group), size = 3) +
    geom_errorbar(aes(xmin =.data$L2.5, xmax = .data$U97.5), size = .8) + xlab("2.5% quantile, Median, 97.5% quantile")+
    ylab("Pollutant Name") + scale_y_discrete(position = weight_axis_pos) +
    theme(panel.grid.major.y = element_line(colour = "grey50", size = 0.2))


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

plotBetas = function(model, filename = NULL, beta_axis_pos = "left") # enter just name of file, not with .pdf
{
  results = model$results[[2]]

  if(length(names(model))==3) # if it was run in parallel
  {
    results = model$model$summary[[2]]
  }

  betas = results[which(substr(row.names(results), 1,1)=="B"), c(1,3,5)]
  phis = results[which(substr(row.names(results), 1,1)=="p"), c(1,3,5)]

  effectEstimates = rbind(betas, phis)

  significant = ifelse(effectEstimates[,1]>0 | effectEstimates[,3]<0, "sig", "insig")

  effectData = data.frame(name = factor(row.names(effectEstimates)), median = effectEstimates[, 2],
                          L2.5 = effectEstimates[,1], U97.5 = effectEstimates[, 3],
                          significant = factor(significant))

  plotBetas = ggplot(data = effectData, mapping = aes(x = .data$median, y = .data$name)) +
    geom_vline(xintercept = 0, linetype = "longdash", colour = "red")+
    geom_point(size = 3, aes(color = .data$significant)) +
    scale_color_manual(values = c("sig" = "red", "insig" = "black")) +
    geom_errorbar(aes(xmin =.data$L2.5, xmax = .data$U97.5), size = .8) + ylab("Effect Name") +
    xlab("2.5% quantile, Median, 97.5% quantile")  +
    scale_y_discrete(position = beta_axis_pos) +
    theme(panel.grid.major.y = element_line(colour = "grey50", size = 0.2))


  print(plotBetas)

  if(!is.null(filename))
  {
    grDevices::pdf(paste(filename, ".pdf", sep = ""))
    print(plotBetas)
    grDevices::dev.off()
  }
}


