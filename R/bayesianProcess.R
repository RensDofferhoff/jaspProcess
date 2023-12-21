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
    # .procModelSummaryTable(jaspResults, options, NULL)
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
  # .procModelSummaryTable(jaspResults, options, modelsContainer)
  # Create container for parameter estimates for each model
  parEstContainer <- .procContainerParameterEstimates(jaspResults, options, modelsContainer)
  # Create tables for parameter estimates
  .procBayesParameterEstimateTables(parEstContainer, options, modelsContainer)
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

  blavaan <- blavaan::blavaan

  fittedModel <- try(blavaan::bsem(
    model           = container[["syntax"]]$object,
    data            = dataset,
    n.chains        = options$mcmcChains,
    burnin          = options$mcmcBurnin,
    sample          = options$mcmcSamples,
    do.fit          = doFit
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
  tbl$addColumnInfo(name = "est",      title = gettext("Estimate"),   type = "number", format = "sf:4;dp:3")
  tbl$addColumnInfo(name = "sd",       title = gettext("SD"), type = "number", format = "sf:4;dp:3")
  tbl$addColumnInfo(name = "ci.lower", title = gettext("Lower"),      type = "number", format = "sf:4;dp:3",
                    overtitle = gettextf("%s%% Credible Interval", options$ciLevel * 100))
  tbl$addColumnInfo(name = "ci.upper", title = gettext("Upper"),      type = "number", format = "sf:4;dp:3",
                    overtitle = gettextf("%s%% Credible Interval", options$ciLevel * 100))
  tbl$addColumnInfo(name = "rhat", title = gettext("R-hat"), type = "number")
  tbl$addColumnInfo(name = "bulkEss", title = gettext("ESS (bulk)"), type = "integer")
  tbl$addColumnInfo(name = "tailEss", title = gettext("ESS (tail)"), type = "integer")

  tbl[["est"]]      <- coefs$mean
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
      probs = ci,
      # This should be zero because as.array does not include warmup samples
      warmup = 0
    )
  )[, c("mean", "sd", ciNames, "Rhat", "Bulk_ESS", "Tail_ESS")]
  
  names(parStats) <- c("mean", "sd", "ci.lower", "ci.upper", "Rhat", "Bulk_ESS", "Tail_ESS")

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
