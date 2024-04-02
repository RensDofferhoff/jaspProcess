#
# Copyright (C) 2023 University of Amsterdam and Netherlands eScience Center
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

# Main function ----
BayesianProcess <- function(jaspResults, dataset = NULL, options) {
  # Set title
  jaspResults$title <- gettext("Process Analysis")
  # Check if analysis has any variables to read in
  if (options[["dependent"]] == "" || (length(options[["covariates"]]) == 0 && length(options[["factors"]]) == 0)) {
    # create empty summary table
    .procBayesModelSummaryTable(jaspResults, options, NULL)
    return()
  }
  options$naAction <- "listwise"
  # Read dataset
  dataset <- .procReadData(options)
  # Check for errors in dataset
  .procErrorHandling(dataset, options)
  # Create a container for each model
  modelsContainer <- .procContainerModels(jaspResults, options)
  # Check if all models are ready to compute something
  if (.procIsReady(options)) {
    # Transform input for each model into a graph for further processing
    .procModelGraph(jaspResults, options)
    # Add factor dummy variables manually to dataset
    # because lavaan does not do it automatically;
    # also add three-way interaction variables manually to dataset 
    # because lavaan does not allow threeway interaction regression terms
    dataset <- .procAddFactorDummyIntVars(jaspResults, dataset, options)
    # Add parameter names to graph of each model
    .procGraphAddParNames(jaspResults, options)
    # Compute quantiles at which to probe moderators for each model
    .procModProbes(jaspResults, dataset, options)
    # Add undirected graph for variances and covariances
    .procResCovGraph(jaspResults, options)
    # Create lavaan syntax for each model from graph
    .procModelSyntax(jaspResults, options)
    # Fit lavaan models based on syntax and dataset and update models container
    modelsContainer <- .procBayesComputeResults(jaspResults, dataset, options)
  }
  # Create container for path plots for each model
  pathPlotContainer <- .procContainerPathPlots(jaspResults, options)
  # Create path plots for each model and add to container
  .procPathPlots(pathPlotContainer, options, modelsContainer)
  # Create table with model fit indices (AIC, ...)
  .procBayesModelSummaryTable(jaspResults, options, modelsContainer)
  # Create container for parameter estimates for each model
  parEstContainer <- .procContainerParameterEstimates(jaspResults, options, modelsContainer)
  # Create tables for parameter estimates
  .procBayesParameterEstimateTables(parEstContainer, options, modelsContainer)

  .procBayesMcmcPlots(jaspResults, options, modelsContainer)

  # Create html output with lavaan syntax for each model
  .procPlotSyntax(jaspResults, options, modelsContainer)

  return()
}

.procBayesComputeResults <- function(jaspResults, dataset, options) {
  modelsContainer <- jaspResults[["modelsContainer"]]
  nModels <- length(options[["processModels"]])

  for (i in 1:nModels) {
    modelOptions <- options[["processModels"]][[i]]
    modelName <- modelOptions[["name"]]

    if (is.null(modelsContainer[[modelName]][["fittedModel"]])) {
      if (modelOptions[["inputType"]] == "inputModelNumber" && !modelOptions[["modelNumber"]] %in% .procHardCodedModelNumbers()) {
        fittedModel <- gettextf("%1$s: Hayes model %2$s not implemented", modelName, modelOptions[["modelNumber"]])
      } else {
        fittedModel <- .procBayesResultsFitModel(
          modelsContainer[[modelName]],
          dataset,
          options
        )
      }
      state <- createJaspState(object = fittedModel)
      modelsContainer[[modelName]][["fittedModel"]] <- state
    }
  }

  return(modelsContainer)
}

.procBayesResultsFitModel <- function(container, dataset, options) {
  # Should model be fitted?
  doFit <- .procCheckFitModel(container[["graph"]]$object)

  if (!doFit) {
    dataset <- NULL
  }

  # Necessary for JASP to find function blavaan
  blavaan <- blavaan::blavaan

  fittedModel <- try(blavaan::bsem(
    model           = container[["syntax"]]$object,
    data            = dataset,
    n.chains        = options$mcmcChains,
    burnin          = options$mcmcBurnin,
    sample          = options$mcmcSamples,
    do.fit          = doFit,
    target          = "stan"
  ))

  if (jaspBase::isTryError(fittedModel)) {
    return(gettextf("Estimation failed: %s", gsub("lavaan ERROR:", "", jaspBase::.extractErrorMessage(fittedModel))))
  }

  if (doFit) {
    container[["graph"]]$object <- .procGraphAddEstimates(container[["graph"]]$object, fittedModel)
    container[["resCovGraph"]]$object <- .procGraphAddEstimates(container[["resCovGraph"]]$object, fittedModel, type = "variances")
  }

  return(fittedModel)
}

# Output functions ----

.procBayesModelSummaryTable <- function(jaspResults, options, modelsContainer) {
  if (!is.null(jaspResults[["modelSummaryTable"]])) return()

  modelNumbers <- lapply(options[["processModels"]], function(mod) {
    graph <- modelsContainer[[mod[["name"]]]][["graph"]]$object
    if (!igraph::is.igraph(graph)) return(NULL)
    return(.procRecognizeModelNumber(graph))
  })

  modelNumberIsValid <- !sapply(modelNumbers, is.null)

  modelNames <- sapply(options[["processModels"]], function(mod) mod[["name"]])

  procResults <- lapply(options[["processModels"]], function(mod) modelsContainer[[mod[["name"]]]][["fittedModel"]]$object)

  # Remove invalid models
  resultIsValid <- sapply(procResults, function(mod) inherits(mod, "lavaan") && mod@Options[["do.fit"]])
  procResults <- procResults[resultIsValid]

  tableRowIsValid <- modelNumberIsValid & resultIsValid

  summaryTable <- createJaspTable(title = gettext("Model summary"), rowNames = modelNames[tableRowIsValid])
  summaryTable$dependOn(c(.procGetDependencies(), "processModels", "naAction"))
  summaryTable$position <- 1

  ovtWaic <- gettext("WAIC")
  ovtLoo  <- gettext("LOO")

  summaryTable$addColumnInfo(name = "Model",        title = "",                         type = "string" )
  summaryTable$addColumnInfo(name = "modelNumber",  title = "Hayes model number",       type = "integer" )
  summaryTable$addColumnInfo(name = "waicEst",      title = gettext("Estimate"),        type = "number", overtitle = ovtWaic )
  summaryTable$addColumnInfo(name = "waicSE",       title = gettext("SE"),              type = "number", overtitle = ovtWaic )
  summaryTable$addColumnInfo(name = "looEst",       title = gettext("Estimate"),        type = "number", overtitle = ovtLoo )
  summaryTable$addColumnInfo(name = "looSE",        title = gettext("SE"),              type = "number", overtitle = ovtLoo )
  
  summaryTable$addColumnInfo(name = "N",            title = gettext("n"),               type = "integer")

  jaspResults[["modelSummaryTable"]] <- summaryTable

  summaryTable[["Model"]]       <- modelNames[tableRowIsValid]
  summaryTable[["modelNumber"]] <- modelNumbers[tableRowIsValid]

  converged <- sapply(procResults, function(mod) mod@Fit@converged)

  if (length(procResults) == 0) {
    summaryTable$addFootnote(message = gettext("At least one model is incomplete or no model is specified. Please add at least one model and complete specified models."))
    return()
  }

  # Use lavaan::fitMeasures instead
  looResults <- lapply(procResults, function(mod) loo::loo(mod@external$mcmcout))
  waicResults <- lapply(procResults, function(mod) loo::waic(loo::extract_log_lik(mod@external$mcmcout)))
  # df  <- sapply(procResults, lavaan::fitMeasures, fit.measures = "df")

  ### Add footnotes about warnings in WAIC and LOO

  summaryTable[["waicEst"]]  <- sapply(waicResults, function(mod) mod$estimates["waic", "Estimate"])
  summaryTable[["waicSE"]]  <- sapply(waicResults, function(mod) mod$estimates["waic", "SE"])
  summaryTable[["looEst"]]  <- sapply(looResults, function(mod) mod$estimates["looic", "Estimate"])
  summaryTable[["looSE"]]  <- sapply(looResults, function(mod) mod$estimates["looic", "SE"])
  summaryTable[["N"]]        <- sapply(procResults, lavaan::lavInspect, what = "nobs")
}

.procBayesParameterEstimateTables <- function(container, options, modelsContainer) {
  if (is.null(modelsContainer)) return()
  
  procResults <- lapply(options[["processModels"]], function(mod) modelsContainer[[mod[["name"]]]][["fittedModel"]]$object)
  modelNames <- sapply(options[["processModels"]], function(mod) mod[["name"]])

  for (i in 1:length(procResults)) {
    if (is.null(container[[modelNames[i]]])) {
      modelContainer <- createJaspContainer(title = modelNames[i], , initCollapsed = TRUE)
      modelContainer$dependOn(
        nestedOptions = .procGetSingleModelsDependencies(as.character(i))
      )
      container[[modelNames[i]]] <- modelContainer
    } else {
      modelContainer <- container[[modelNames[i]]]
    }
    
    .procSetContainerError(modelContainer, procResults[[i]])

    if (options[["processModels"]][[i]][["pathCoefficients"]])
      .procBayesPathCoefficientsTable(modelContainer, options, procResults[[i]], i)

    if (options[["processModels"]][[i]][["mediationEffects"]])
      .procBayesPathMediationEffectsTable(modelContainer, options, procResults[[i]], i)

    # if (options[["processModels"]][[i]][["totalEffects"]])
    #   .procPathTotalEffectsTable(modelContainer, options, procResults[[i]], i)

    # if (options[["processModels"]][[i]][["residualCovariances"]])
    #   .procCovariancesTable(modelContainer, options, procResults[[i]], i)
  }
}

.procBayesCoefficientsTable <- function(tbl, options, coefs) {
  titlePosterior <- gettext("Posterior")
  titleCI <- gettextf("%s%% Credible Interval", options$ciLevel * 100)

  tbl$addColumnInfo(name = "mean",     title = gettext("Estimate"),   type = "number", format = "sf:4;dp:3", overtitle = titlePosterior)
  tbl$addColumnInfo(name = "median",   title = gettext("Median"),     type = "number", format = "sf:4;dp:3", overtitle = titlePosterior)
  tbl$addColumnInfo(name = "sd",       title = gettext("SD"), type = "number", format = "sf:4;dp:3", overtitle = titlePosterior)
  tbl$addColumnInfo(name = "ci.lower", title = gettext("Lower"),      type = "number", format = "sf:4;dp:3",
                    overtitle = titleCI)
  tbl$addColumnInfo(name = "ci.upper", title = gettext("Upper"),      type = "number", format = "sf:4;dp:3",
                    overtitle = titleCI)
  tbl$addColumnInfo(name = "rhat", title = gettext("R-hat"), type = "number")
  tbl$addColumnInfo(name = "bulkEss", title = gettext("ESS (bulk)"), type = "integer")
  tbl$addColumnInfo(name = "tailEss", title = gettext("ESS (tail)"), type = "integer")

  tbl[["mean"]]      <- coefs$mean
  tbl[["median"]]   <- coefs$median
  tbl[["sd"]]       <- coefs$sd
  tbl[["ci.lower"]] <- coefs$ci.lower
  tbl[["ci.upper"]] <- coefs$ci.upper
  tbl[["rhat"]]     <- coefs$Rhat
  tbl[["bulkEss"]]  <- coefs$Bulk_ESS
  tbl[["tailEss"]]  <- coefs$Tail_ESS

  if (options$standardizedEstimates != "unstandardized") {
    txt <- switch(options$standardizedEstimates,
      centered = gettext("mean-centered"),
      standardized = gettext("standardized")
    )
    tbl$addFootnote(gettextf("Summary statistics are based on %s estimates.", txt))
  }
}

.procBayesGetParStats <- function(stanFit, parTable, options) {
  ci <- c((1-options$ciLevel)/2, 1-(1-options$ciLevel)/2)
  ciNames <- paste0(100*ci, "%")

  includePars <- parTable$stanpnum
  
  if (all(is.na(includePars))) {
    includePars <- 1:nrow(parTable)
  }

  parStats <- as.data.frame(
    rstan::monitor(
      as.array(stanFit)[, , na.omit(includePars), drop = FALSE],
      probs = c(0.5, ci), # also compute median
      # This should be zero because as.array does not include warmup samples
      warmup = 0
    )
  )[, c("mean", "50%", "sd", ciNames, "Rhat", "Bulk_ESS", "Tail_ESS")]
  
  names(parStats) <- c("mean", "median", "sd", "ci.lower", "ci.upper", "Rhat", "Bulk_ESS", "Tail_ESS")

  return(parStats)
}

.procBayesPathCoefficientsTable <- function(container, options, procResults, modelIdx) {
  if (!is.null(container[["pathCoefficientsTable"]])) return()

  pathCoefTable <- createJaspTable(title = gettext("Path coefficients"))
  pathCoefTable$dependOn(
    options = "parameterLabels",
    nestedOptions = list(c("processModels", as.character(modelIdx), "pathCoefficients"))
  )
  container[["pathCoefficientsTable"]] <- pathCoefTable

  if (!.procIsValidModel(container, procResults)) return()

  parTable <- lavaan::parameterTable(procResults)
  parTable <- parTable[parTable$op == "~",]

  parStats <- .procBayesGetParStats(as.array(procResults@external$mcmcout), parTable, options)

  pathCoefTable$addColumnInfo(name = "lhs", title = "", type = "string")
  pathCoefTable$addColumnInfo(name = "op",  title = "", type = "string")
  pathCoefTable$addColumnInfo(name = "rhs", title = "", type = "string")

  pathCoefTable[["lhs"]] <- gsub("__", ":", parTable$rhs)
  pathCoefTable[["op"]]  <- rep("\u2192", nrow(parTable))
  pathCoefTable[["rhs"]] <- gsub("__", ":", parTable$lhs)

  if (options$parameterLabels) {
    pathCoefTable$addColumnInfo(name = "label", title = gettext("Label"), type = "string")
    pathCoefTable[["label"]] <- parTable$label
  }

  .procBayesCoefficientsTable(pathCoefTable, options, parStats)
}

.procBayesAddMedSamples <- function(stanFit, parTable) {
  rhs <- parTable$rhs[parTable$op == ":="]

  samples <- lapply(parTable$pxnames[parTable$op == "~"], function(nm) {
    return(as.array(stanFit)[, , nm])
  })
  
  # rstan::extract(stanFit, pars = parTable$pxnames[parTable$op == "~"], inc_warmup = TRUE)

  names(samples) <- parTable$label[parTable$op == "~"]

  medEffectSamples <- sapply(rhs, function(trm) {
    return(eval(str2lang(trm), envir = samples))
  })

  return(array(medEffectSamples, dim = c(dim(as.array(stanFit))[1:2], length(rhs)), dimnames = list(NULL, NULL, parameters = parTable$label[parTable$op == ":="])))
}

.procBayesPathMediationEffectsTable <- function(container, options, procResults, modelIdx) {
  if (!is.null(container[["mediationEffectsTable"]])) return()

  medEffectsTable <- createJaspTable(title = gettext("Mediation effects"))
  medEffectsTable$dependOn(
    options = c("parameterLabels", "moderationProbes"),
    nestedOptions = list(c("processModels", as.character(modelIdx), "mediationEffects"))
  )

  container[["mediationEffectsTable"]] <- medEffectsTable

  if (!.procIsValidModel(container, procResults)) return()

  parTable <- lavaan::parameterTable(procResults)
  medEffects <- parTable[parTable$op == ":=",]

  labelSplit <- lapply(strsplit(medEffects$lhs, "\\."), strsplit, split = "__")

  # Only use label splits of length > 1 to omit total effects
  isMediation <- sapply(labelSplit, function(path) length(path[[1]]) > 1)

  labelSplit <- labelSplit[isMediation]
  medEffects <- medEffects[isMediation,]

  parStats <- .procBayesGetParStats(.procBayesAddMedSamples(procResults@external$mcmcout, parTable), medEffects, options)

  # Get paths from label of mediation effect
  medPaths <- lapply(labelSplit, function(path) path[[1]])

  # Get path lengths
  medPathLengths <- sapply(medPaths, length)

  # Sort paths to incresaing length
  medLengthSortIdx <- sort(medPathLengths, index.return = TRUE)$ix
  medEffects <- medEffects[medLengthSortIdx, ]

  # Add a column for each step of longest path
  for (i in 1:max(medPathLengths)) {
    # If path has step add var name otherwise empty
    medEffect <- sapply(medPaths[medLengthSortIdx], function(path) ifelse(length(path) >= i, path[i], ""))

    # Add operator columns
    if (i > 1) {
      # Add operator for non-empty path steps otherwise empty
      medOp <- ifelse(medEffect == "", "", "\u2192")
      medEffectsTable$addColumnInfo(name = paste0("op_", i), title = "", type = "string")
      medEffectsTable[[paste0("op_", i)]] <- medOp
    }

    medEffectsTable$addColumnInfo(name = paste0("lhs_", i), title = "", type = "string")
    medEffectsTable[[paste0("lhs_", i)]] <- medEffect
  }

  medEffectIsConditional <- sapply(labelSplit, function(path) length(path) > 1)

  uniqueMods <- unique(unlist(lapply(labelSplit[medEffectIsConditional], function(path) lapply(path[-1], function(row) row[1]))))

  modProbes <- .procEffectsTablesGetConditionalLabels(labelSplit[medEffectIsConditional], uniqueMods)

  for (mod in uniqueMods) {
    medEffectsTable$addColumnInfo(name = mod, title = mod, type = "string", combine = FALSE) # combine = F because empty cells indicate no moderation
    modLabels <- vector("character", length(medEffectIsConditional))
    modLabels[medEffectIsConditional] <- modProbes[[mod]]
    medEffectsTable[[mod]] <- modLabels
  }

  # Add column with parameter labels
  if (options$parameterLabels) {
    medEffectsTable <- .procEffectsTablesParameterLabels(medEffectsTable, medEffects)
  }

  .procBayesCoefficientsTable(medEffectsTable, options, parStats)
}

.procBayesMcmcPlots <- function(container, options, modelsContainer) {

  if (is.null(container[["mcmcPlotContainer"]])) {
    plotContainer <- createJaspContainer(
      title = gettext("MCMC Plots"),
      dependencies = c("monitoredParametersShown", "useColorPalette", "colorPalette")
    )
    container[["mcmcPlotContainer"]] <- plotContainer
  } else {
    plotContainer <- container[["mcmcPlotContainer"]]
  }

  procResults <- lapply(options[["processModels"]], function(mod) modelsContainer[[mod[["name"]]]][["fittedModel"]]$object)
  modelNames <- sapply(options[["processModels"]], function(mod) mod[["name"]])

  for (i in 1:length(procResults)) {
    if (is.null(plotContainer[[modelNames[i]]])) {
      modelContainer <- createJaspContainer(title = modelNames[i], , initCollapsed = TRUE)
      modelContainer$dependOn(
        nestedOptions = .procGetSingleModelsDependencies(as.character(i))
      )
      plotContainer[[modelNames[i]]] <- modelContainer
    } else {
      modelContainer <- plotContainer[[modelNames[i]]]
    }
    
    if (!.procIsValidModel(modelContainer, procResults[[i]])) next

    parTbl <- lavaan::parTable(procResults[[i]])
    # Only plot free parameters but not fixed
    parTbl <- parTbl[parTbl$free > 0, ]
    # Names of parameters in MCMC output
    mcmcParams <- parTbl[, "pxnames"]
    # Names of parameters to display
    dispParams <- paste(parTbl$rhs, ifelse(parTbl$op == "~", "\u2192", "\u2194"), parTbl$lhs)
    # Get MCMC samples
    mcmcArray <- as.array(procResults[[i]]@external$mcmcout)
    # Replace parameter names
    dimnames(mcmcArray)[[3]][dimnames(mcmcArray)[[3]] %in% mcmcParams] <- dispParams

    # mcmcArray <- abind::abind(mcmcArray, .procBayesAddMedSamples(procResults[[i]]@external$mcmcout, parTbl), along = 3)
    
    # Create dummy mcmcResult for JAGS functions
    mcmcResult <- list(
      samples = lapply(seq(dim(mcmcArray)[2]), function(i) mcmcArray[, i, ])
    )

    # Plotting functions from jaspJags module
    containerObj <- jaspJags:::.JAGSInitPlotsContainers(modelContainer, options, dispParams)

    if (options[["useColorPalette"]]) {
      colorpalette <- options[["colorPalette"]]
      oldColorpalette <- jaspGraphs::getGraphOption("palette")
      on.exit(jaspGraphs::setGraphOption("palette", oldColorpalette))
      jaspGraphs::setGraphOption("palette", colorpalette)
    }

    jaspJags:::.JAGSFillPlotContainers(containerObj, options, mcmcResult, dispParams)

    jaspJags:::.JAGSPlotBivariateScatter(modelContainer, options, mcmcResult, dispParams)
  }
}
