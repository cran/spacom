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


## provides a DescribeExactObject with basic consistency checks
MakeDescribeExactObject <- function(contextual.data,
                                           context.id,
                                           contextual.names,
                                           contextual.weight.matrices,
                                           obj=NULL){
  ## create an empty DescribeExactObject
  pedantic.check <- FALSE
  if (is.null(obj)){
    obj <- new("DescribeExactObject")
    pedantic.check <- TRUE
  }

  ## make sure contextual.data is a data.frame
  obj@contextual.data <- checkType(contextual.data, c("data.frame"),
                                   "contextual.data")

  ## extract number of upper level units
  obj@nb.area <- length(levels(as.factor(contextual.data[[context.id]])))
  if (pedantic.check && nrow(obj@contextual.data) != obj@nb.area) {
    stop("The contextual data has to have as many rows as there are areas, ",
         .obj@nb.area,". ",
         "You specified a frame with ", nrow(obj@contextual.data), " rows.")}
  
  ## make sure context.id is a name in contextual.data
  obj@context.id <- checkContextId(context.id, NULL,
                                   names(obj@contextual.data),
                                   NULL, TRUE)

  ## make sure the contextual names are in a list
  obj <- checkContextualNames(obj, contextual.names, names(contextual.data))
  
  ## make sure the weight matrix list is suitable
  obj@contextual.weight.matrices <-
    checkAllWeightMatrices(contextual.weight.matrices, obj@contextual.names,
                           obj@nb.area, obj@nb.analyses)

  return(obj)
}

performDescribeExact <- function(obj){
return(performAllSpatialWeighting(obj@contextual.names,
                                  obj@context.id,
                                  obj@nb.area,
                                  obj@nb.analyses,
                                  obj@contextual.data,
                                  obj@contextual.weight.matrices))
}

DescribeExact <- function(contextual.data,
                                 context.id,
                                 contextual.names,
                                 contextual.weight.matrices){
  obj <-MakeDescribeExactObject(contextual.data,
                                       context.id,
                                       contextual.names,
                                       contextual.weight.matrices)
  output <- performDescribeExact(obj)
  return(output)
}



## provides a DescribeExactObject with basic consistency checks
MakeDescribeAggregateObject <- function(contextual.data,
                                               context.id,
                                               contextual.names,
                                               contextual.weight.matrices,
                                               
                                               nb.resamples,
                                               aggregation.functions,
                                               confidence.intervals,
                                               individual.weight.names,
                                               sample.seed=NULL){
  ## use the parent class' maker function
  obj <- MakeDescribeExactObject(contextual.data,
                                        context.id,
                                        contextual.names,
                                        contextual.weight.matrices,
                                        new("DescribeAggregateObject"))
  ##
  obj@nb.resamples <- as.integer(nb.resamples)

  ## make sure the contextual names are in a list
  obj <- checkContextualNames(obj, contextual.names, names(contextual.data))
  

  ## make sure the aggregation function list is suitable
  obj@aggregation.functions <- checkAggregationFunctions(aggregation.functions,
                                                         obj@nb.analyses,
                                                         obj@aggregation.names)

  ## make sure the confidence intervals are correctly defined and contain the
  ## median
  obj@percentiles <-
    checkConfidenceIntervals(confidence.intervals, obj@nb.resamples)

  ## make sure the individual.weight.names are suitable
  obj@individual.weight.names <-
    checkAllIndividualWeights(individual.weight.names,
                              obj@contextual.names,
                              names(obj@contextual.data),
                              obj@nb.analyses)

  ## make sure the sample seed is ok
  obj@sample.seed <- checkSeed(sample.seed)

  return(obj)
}

print.descr <- function(object) {
  cat("aggregated data with descriptives:\n")
  print(object@frames)
  cat("\nslots for direct acces to data:\n  usage: <object.name>@<slot.name>\n")
  is.print <- TRUE
  for (name in slotNames(object)) {
    if (is.print) {
      print(name)
    } else {
      show(name)
    }
  }
}
setMethod("show", signature="DescribeAggregateOutput",definition=print.descr)


performDescribeAggregate <- function(obj, exploratory=FALSE) {
  ## creation of an output data to bundle all the return data
  output.obj <- new("DescribeAggregateOutput")

  ## temporary array object to gather the results
  aggregation.results <- array(dim=c(obj@nb.area, obj@nb.resamples, obj@nb.analyses),
                               dimnames=c('row', 'col', 'var'))

  ## inserts the analyses at the right position in the results
  fillAggregationResults <- function(frame, index){
    for (column in 2:ncol(frame)) {
      aggregation.results[,index, column-1] <<- frame[[column]]
    }
  }

  ## custom stateful iterator object which yields a different resample everytime
  ## it is used in the loop
  iterator <- sampleIterator(obj@sample.seed,
                             obj@contextual.data[[obj@context.id]],
                             obj@nb.resamples)
  output.obj@seed <- obj@sample.seed


  ## loops through all resamples of the context data and aggregates them at the
  ## upper level. The results are stored in a list of matrices, each containing
  ## the resampled data for one variable
  gc(FALSE)
  start.time <- proc.time()[3]
  analyses <- foreach(slice=iterator,
                      column=1:obj@nb.resamples,
                      .inorder=TRUE) %do% {
    

    frame <- 
      performAggregation(obj@contextual.names,
                         obj@context.id,
                         obj@nb.area,
                         obj@nb.analyses,
                         obj@individual.weight.names,
                         obj@contextual.data[slice,],
                         obj@aggregation.functions,
                         obj@contextual.weight.matrices)
    fillAggregationResults(frame, column)
    elapsed.time <- proc.time()[3] - start.time
    cat("\rcomputed step ", column, " of ", obj@nb.resamples,
        ". ETA = ", (obj@nb.resamples/column-1)*elapsed.time)

    NULL
  }
  gc(FALSE)
  cat("\rdescription done                                                   \n")
  if (exploratory) {
    return(matrix(aggregation.results, nrow=obj@nb.area, ncol=obj@nb.resamples))
  }
  variable.names <- names(frame)[1:obj@nb.analyses+1]

  for (i in 1:length(variable.names)) {
    name <- variable.names[i]
    ## add the computed matrices to the output object
    output.obj@aggregated.samples[[name]] <-
      matrix(aggregation.results[,,i],
              nrow=obj@nb.area)
    
    ## here, loop through the upper level aggregation functions (mean, sd,
    ## median, percentiles)
    frame = data.frame(1:obj@nb.area)
    names(frame) <- obj@context.id
    ## add mean(mandatory)
    frame[['mean']] <-
      rowMeans(aggregation.results[,,i])
    ## add standard deviation(mandatory)
    frame[['sd']] <-
      apply(aggregation.results[,,i],1,sd)
    ## add percentiles
    perc.frame <- t(apply(aggregation.results[,,i], 1, quantile, probs=obj@percentiles))
    output.obj@frames[[name]] <- cbind(frame, perc.frame)
  }
  return(output.obj)
}

DescribeAggregate <- function(contextual.data,
                                     context.id,
                                     contextual.names,
                                     contextual.weight.matrices,

                                     nb.resamples=1000,
                                     aggregation.functions='mean',
                                     confidence.intervals=.95,
                                     individual.weight.names=NULL,
                                     sample.seed=NULL){
  input.obj <- MakeDescribeAggregateObject(contextual.data,
                                                  context.id,
                                                  contextual.names,
                                                  contextual.weight.matrices,
                                                  
                                                  nb.resamples,
                                                  aggregation.functions,
                                                  confidence.intervals,
                                                  individual.weight.names,
                                                  sample.seed)
  return(performDescribeAggregate(input.obj))
}
  
