################################################################################
################################################################################
## author Till Junge <till.junge@unil.ch>                                     ##
##                                                                            ##
## Copyright (c) UNIL (Universite de Lausanne)                                ##
## NCCR - LIVES (National Centre of Competence in Research "LIVES –           ##
## Overcoming vulnerability: life course perspectives",                       ##
## <http://www.lives-nccr.ch/>)                                               ##
##                                                                            ##
## spacom is free software: you can redistribute it and/or modify it under    ##
## the terms of the GNU General Public License as published by the Free       ##
## Software Foundation, either version 2 of the License or any later version. ##
##                                                                            ##
## spacom is distributed in the hope that it will be useful, but WITHOUT ANY  ##
## WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS  ##
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more      ##
## details, see <http://www.gnu.org/licenses/>.                               ##
################################################################################
################################################################################


## provides a ResampleExactSpawMLObject needed for the
## analysis and performs all consistency checks
MakeResampleExactSpawMLObject <-function(individual.level.data,
                                                 contextual.names,
                                                 context.id,
                                                 formula,
                                                 contextual.data=NULL,
                                                 precise.data=NULL,
                                                 contextual.weight.matrices=NULL,
                                                 individual.weight.names=NULL,
                                                 aggregation.functions="mean",
                                                 
                                                 confidence.intervals,
                                                 nb.resamples=1000,
                                                 individual.sample.seed=NULL,
                                                 obj=NULL) {

  ## create an empty ResampleExactSpawMLObject to fill, using the
  ## parent class SpawMLObject
  if (is.null(obj)){
    obj <- new("ResampleExactSpawMLObject")
  }
  obj <-
    MakeSpawMLObject(individual.level.data=individual.level.data,
                                 contextual.names=contextual.names,
                                 context.id=context.id,
                                 formula=formula,
                                 contextual.data=contextual.data,
                                 precise.data=precise.data,
                                 contextual.weight.matrices=contextual.weight.matrices,
                                 individual.weight.names=individual.weight.names,
                                 aggregation.functions=aggregation.functions,

                                 obj)
  
  ## makes sure the number of resamples is a positive integer
  obj@nb.resamples <- checkNbResamples(nb.resamples)
  
  ## make sure the sample seed is ok
  obj@individual.sample.seed <- checkSeed(individual.sample.seed)

  ## make sure the confidence intervals are correctly defined and contain the
  ## median
  obj@percentiles <-
    checkConfidenceIntervals(confidence.intervals, obj@nb.resamples)

  
  return(obj)
}

## ## Definition of an output object for describeResampledIndividualContext
printSraweObject <- function(x, is.print=TRUE) {
  cat("New model result \n")
  if (!("contextual.sample.seed" %in% slotNames(x))) {
    cat("Stratified resampling for precise contextual values\n")
  } else {
    cat("Stratified resampling for contextual variables aggregated with error")
  }
  cat(paste("\nNumber of resamples: ", x@nb.resamples, "\n"))
  cat("\nLinear mixed model fit by REML:\n")
  if (is.print) {
    print(x@model.fit)
  } else {
    show(x@model.fit)
  }
  cat("\nFixed effects:\n")
  if (is.print) {
    print(x@fixed)
  } else {
    show(x@fixed)
  }
  cat("\nRandom effects:\n")
  if (is.print) {
    print(x@random.var)
  } else {
    show(x@random.var)
  }
  cat("\nStandardised fixed effects:\n")
  if (TRUE==all.equal(dim(x@betas), c(0,0))) {
    cat("None\n")
  } else {  
    if (is.print) {
      print(x@betas)
    } else {
      show(x@betas)
    }
  }
  cat("\nslots for direct acces to data:\n  usage: <object.name>@<slot.name>\n")
  for (name in slotNames(x)) {
    if (is.print) {
      print(name)
    } else {
      show(name)
    }
  }
}

setMethod("print", signature="ResampleExactSpawMLOutput",
          definition=printSraweObject)
setMethod("show", signature="ResampleExactSpawMLOutput",
          definition=function(object) { printSraweObject(object, is.print=FALSE)})

PerformResampledMultilevel <- function(obj, srawe=TRUE, exploratory=FALSE){
  ## creation of an output data to bundle all the return data
  if (srawe) {
    output.obj <- new("ResampledSpawMLOutput")
  } else {
    output.obj <- new("ResampleExactSpawMLOutput")
  }
  output.obj@nb.resamples <- obj@nb.resamples

  ## temporary data structures to gather the results
  aic.results <- matrix(nrow=5, ncol=obj@nb.resamples)
  row.names(aic.results) <- c('AIC', 'BIC', 'logLik', 'deviance', 'REMLdev')
  fixed.effect.results <- NULL # will have to be filled at first occurence
  random.var.results <- NULL # will have to be filled at first occurence
  output.obj@ranefs <- matrix(nrow=obj@nb.area, ncol=obj@nb.resamples)
  beta.results <- NULL # will have to be filled at first occurence

  ##
  n <- 0
  nn <- 0
  prepareRandomEffects <- function(lme.obj, column) {
    vc <- lme4::VarCorr(lme.obj)
    if (is.null(random.var.results)) {
      n <<- nrow(vc[[1]])
      nn <<- n*(n-1)/2

      colnames <- character((1+n)*n+1) # to store 2 sd, (n*n-n)/2+n cov, (n*n-n)/2 corr and sc

      colnames[1] <- 'sd.e'
      rand.names <-row.names(vc[[1]])
      colnames[2:(n+1)] <-
        lapply(rand.names,
               FUN=function(x){return(paste("VAR.", x, sep=""))})
      colnames[(n+2):(2*n+1)] <- 
        lapply(rand.names,
               FUN=function(x){return(paste("SD.", x, sep=""))})

      if (n>1) {
        ind <- 2*n+1
        for (i in 2:n) {
          for (j in 1:(i-1)) {
            ind <- ind + 1
            colnames[ind] <-
              paste( "COV.", rand.names[i], ".", rand.names[j], sep="")
            colnames[ind+nn] <-
              paste("CORR.", rand.names[i], ".", rand.names[j], sep="")
          }
        }
      }

      random.var.results <<- matrix(nrow=(1+n)*n+1, ncol=obj@nb.resamples)
      row.names(random.var.results) <<- colnames
    }

    random.var.results[1, column] <<- attr(vc, "sc")
    mat <- vc[[1]]
    random.var.results[2:(n+1), column] <<- diag(mat)
    random.var.results[(n+2):(2*n+1), column] <<- attr(mat,"stddev")

    if (n>1) {
      ind <- 2*n+1
      for (i in 2:n) {
        for (j in 1:(i-1)) {
          ind <- ind + 1
          random.var.results[ind, column] <<- mat[i,j]
          random.var.results[ind+nn, column] <<- attr(mat, "correlation")[i,j]
        }
      }
    }      
  }

  ## inserts the analyses at the right positions in the result matrices
  fillMatrices <- function(lme.obj, beta, index) {
    if (is.null(fixed.effect.results)) {
      ##create the fixed.effect.results matrix of the right dims and name rows
      fix <- lme4::fixef(lme.obj)
      fixed.effect.results <<- matrix(nrow=length(fix),
                                     ncol=obj@nb.resamples)
      row.names(fixed.effect.results) <<- names(fix)

      ## TODO check whether we can make beta conditional when used with empty model
      if (length(beta)>0) {
        beta.results <<- matrix(nrow=length(beta), ncol=obj@nb.resamples)
      }
    }
    fixed.effect.results[, index] <<- lme4::fixef(lme.obj)
    prepareRandomEffects(lme.obj, index)

    aic.results[1, index] <<- lme4::AIC(lme.obj)
    aic.results[2, index] <<- lme4::BIC(lme.obj)
    aic.results[3, index] <<- as.numeric(logLik(lme.obj))
    aic.results[4, index] <<- deviance(lme.obj, REML=FALSE)
    aic.results[5, index] <<- deviance(lme.obj, REML=TRUE)
    ## TODO check whether we can make beta conditional when used with empty model
    if (!is.null(beta.results)){
      beta.results[, index] <<- beta
      rownames(beta.results) <<- names(beta)
    }
    
    
    output.obj@ranefs[ ,index] <<- lme4::ranef(lme.obj)[[1]][,1]
  }
  
  ## custom stateful iterator object which yields a different resample everytime
  ## it is used in the loop
  individual.iterator <-
    sampleIterator(obj@individual.sample.seed,
                   obj@individual.level.data[[obj@context.id]],
                   obj@nb.resamples)
  
  output.obj@individual.sample.seed <- obj@individual.sample.seed
  if (exploratory){
    output.obj@contextual.sample.seed <- obj@contextual.sample.seed
  } else {
    if (srawe) {
      if (obj@same.survey) {
        output.obj@contextual.sample.seed <- output.obj@individual.sample.seed
      } else {
        contextual.iterator <- sampleIterator(obj@contextual.sample.seed,
                                              GetContext(obj)[[obj@context.id]],
                                              obj@nb.resamples)
        
        output.obj@contextual.sample.seed <- obj@contextual.sample.seed
      } 
    }
  }
  
  ## loops through all resamples of the context data and aggregates them at the
  ## upper level. The results are stored in a list of matrices, each containing
  ## the resampled data for one variable
  start.time <- proc.time()[3]
  if (srawe & !exploratory) {
    
    analyses <- foreach(individual.slice=individual.iterator,
                        contextual.slice=contextual.iterator,
                        column=1:obj@nb.resamples,
                        .inorder=TRUE) %do% {
      sml.obj <-
        new("SpawMLObject",
            individual.level.data = obj@individual.level.data[individual.slice,],
            contextual.data = GetContext(obj)[contextual.slice,],
            precise.data = obj@precise.data,
            context.id = obj@context.id,
            contextual.names = obj@contextual.names,
            aggregation.names = obj@aggregation.names,
            precise.names = obj@precise.names,
            individual.weight.names = obj@individual.weight.names,
            aggregation.functions = obj@aggregation.functions,
            contextual.weight.matrices = obj@contextual.weight.matrices,
            formula = obj@formula,
            nb.area = obj@nb.area,
            nb.analyses = obj@nb.analyses,
            nb.aggregations= obj@nb.aggregations,
            nb.precise.weightings=obj@nb.precise.weightings)
      lme <- PerformSpawML(sml.obj)
      fillMatrices(lme@lme, lme@beta, column)
      elapsed.time <- proc.time()[3]-start.time
      cat("\rcomputed step ", column, " of ", obj@nb.resamples,
              ". ETA = ", (obj@nb.resamples/column-1)*elapsed.time)
      NULL
    }
  } else if (!srawe & !exploratory) {
    analyses <- foreach(individual.slice=individual.iterator,
                        column=1:obj@nb.resamples,
                        .inorder=TRUE) %do% {
      sml.obj <-
        new("SpawMLObject",
            individual.level.data = obj@individual.level.data[individual.slice,],
            contextual.data = GetContext(obj),
            precise.data = obj@precise.data,
            context.id = obj@context.id,
            contextual.names = obj@contextual.names,
            aggregation.names = obj@aggregation.names,
            precise.names = obj@precise.names,
            individual.weight.names = obj@individual.weight.names,
            aggregation.functions = obj@aggregation.functions,
            contextual.weight.matrices = obj@contextual.weight.matrices,
            formula = obj@formula,
            nb.area = obj@nb.area,
            nb.analyses = obj@nb.analyses,
            nb.aggregations= obj@nb.aggregations,
            nb.precise.weightings=obj@nb.precise.weightings)
      lme <- PerformSpawML(sml.obj)
      fillMatrices(lme@lme, lme@beta, column)
      elapsed.time <- proc.time()[3] - start.time
      cat("\rcomputed step ", column, " of ", obj@nb.resamples,
              ". ETA = ", (obj@nb.resamples/column-1)*elapsed.time)
      NULL
    }
  } else {
    analyses <- foreach(individual.slice=individual.iterator,
                        column=1:obj@nb.resamples,
                        .inorder=TRUE) %do%
    {
      ## build a data.frame out of the sample matrix
      contextual.data <- data.frame(obj@contextual.data[, column])
      names(contextual.data) <- substr(obj@contextual.names[[1]], 1,nchar(obj@contextual.names[[1]])-5)
      contextual.data[[obj@context.id]] <- 1:obj@nb.area
      sml.obj <-
        new("SpawMLObject",
            individual.level.data = obj@individual.level.data[individual.slice,],
            contextual.data = NULL,
            precise.data = contextual.data,
            context.id = obj@context.id,
            contextual.names = obj@contextual.names,
            aggregation.names = list(),
            precise.names = obj@contextual.names,
            individual.weight.names = obj@individual.weight.names,
            aggregation.functions = obj@aggregation.functions,
            contextual.weight.matrices = obj@contextual.weight.matrices,
            formula = obj@formula,
            nb.area = obj@nb.area,
            nb.analyses = obj@nb.analyses,
            nb.aggregations= as.integer(0),
            nb.precise.weightings=obj@nb.precise.weightings)
      lme <- PerformSpawML(sml.obj)
      fillMatrices(lme@lme, lme@beta, column)
      elapsed.time <- proc.time()[3] - start.time
      cat("\rcomputed step ", column, " of ", obj@nb.resamples,
              ". ETA = ", (obj@nb.resamples/column-1)*elapsed.time)
      NULL
    }
  }
  cat("\rspacom done                                                        \n")
  tmp.list <- list(fixed=fixed.effect.results,
                   model.fit=aic.results,
                   random.var=random.var.results)
  if (!is.null(beta.results)) {
    tmp.list[["betas"]] <- beta.results
  }
  
  for (name in names(tmp.list)) {
    ## compute mean (mandatory)
    mu <- rowMeans(tmp.list[[name]])
    frame <- data.frame(mean=mu)
    ## compute standard deviation (mandatory)
    frame$sd <- apply(tmp.list[[name]], 1, sd)
    ## compute percentiles
    perc.frame <- t(apply(tmp.list[[name]], 1, quantile, probs=obj@percentiles))
    slot(output.obj, name) <- cbind(frame, perc.frame)
  }

  return(output.obj)  
}

## PerformResampleExactSpawML <- function(obj){
##   ## creation of an output data to bundle all the return data
##   output.obj <- new("ResampleExactSpawMLOutput")
## 
##   ## custom stateful iterator object which yields a different resample everytime
##   ## it is used in the loop
##   iterator <- sampleIterator(obj@individual.sample.seed,
##                              obj@individual.level.data[[obj@context.id]],
##                              obj@nb.resamples)
## 
##   output.obj@samples <- as.list(iterator)
## 
##   
##   ## loops through all resamples of the context data and aggregates them at the
##   ## upper level. The results are stored in a list of matrices, each containing
##   ## the resampled data for one variable
##   analyses <- foreach(slice=output.obj@samples,
##                       column=1:obj@nb.resamples,
##                       .inorder=TRUE) %do% {
##     sml.obj <- new("SpawMLObject",
##                    individual.level.data = obj@individual.level.data[slice,],
##                    contextual.data = GetContext(obj),
##                    precise.data = obj@precise.data,
##                    context.id = obj@context.id,
##                    contextual.names = obj@contextual.names,
##                    individual.weight.names = obj@individual.weight.names,
##                    aggregation.functions = obj@aggregation.functions,
##                    contextual.weight.matrices = obj@contextual.weight.matrices,
##                    formula = obj@formula,
##                    nb.area = obj@nb.area,
##                    nb.analyses = obj@nb.analyses)
##    PerformSpawML(sml.obj)
##   }
##   
##   return(analyses)  
## }

ResampleExactSpawML <-function(
           individual.level.data,
           contextual.names,
           context.id,
           formula,
           ## contextual.data=NULL,
           precise.data,
           contextual.weight.matrices=NULL,
           ## individual.weight.names=NULL,
           ## aggregation.functions="mean",

           confidence.intervals=c(.95),
           nb.resamples=1000,
           individual.sample.seed=NULL) {
  obj <-
    MakeResampleExactSpawMLObject(individual.level.data,
                                                    contextual.names,
                                                    context.id,
                                                    formula,
                                                    contextual.data=NULL,
                                                    precise.data,
                                                    contextual.weight.matrices,
                                                    individual.weight.names=NULL,
                                                    aggregation.functions="mean",

                                                    confidence.intervals,
                                                    nb.resamples,
                                                    individual.sample.seed)
  model <- PerformResampledMultilevel(obj, srawe=FALSE)
}
