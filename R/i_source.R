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


MakeResampleExploreSpawMLObject <-
  function(individual.level.data,
           contextual.name,
           contextual.data,
           context.id,
           nb.resamples,
           formula,
           distance.matrix,
           multilevel.bandwidths,
           individual.weight.names=NULL,
           aggregation.function="mean",
           confidence.intervals=c(.95),
           individual.sample.seed=NULL,
           contextual.sample.seed=NULL,
           kernel=NULL)
{
  obj <- new("ResampleExploreSpawML")
  ## make sure individual.level.data is a data.frame
  obj@individual.level.data <- checkType(individual.level.data, c("data.frame"),
                                         "individual.level.data")

  ## make sure the number of resample is an integer > 1
  obj@nb.resamples <- checkNbResamples(nb.resamples)
  
  ## make sure contextual.data is either NULL or a data.frame
  obj@contextual.data <- checkType(obj=contextual.data,
                                   classes=c("NULL", "data.frame"),
                                   name="contextual.data")
  contextual.variables <- names(GetContextMethod(obj))
  obj <- checkContextualNames(obj, contextual.name, contextual.variables)

  obj@individual.weight.name <-
    checkIndividualWeights(individual.weight.names, obj@contextual.names,
                           contextual.variables)


  ## make sure context.id is a name in contextual.data
  obj@context.id <- checkContextId(context.id, names(obj@individual.level.data),
                                   contextual.variables,
                                   NULL)
  ## extract number of upper level units
  obj@nb.area <- length(levels(as.factor(GetContextMethod(obj)[[context.id]])))

  ## make sure the formula isn't a string anymore
  obj@formula <- checkFormula(formula)

  ## check the distance matrix for consistency
  if (is(distance.matrix, "data.frame")) {
    distance.matrix <- as.matrix(distance.matrix)
  }
  if (!is(distance.matrix, "matrix")) {
    stop("The distance matrix has to be of class 'matrix'. You specified an ",
         "object of class '", class(distance.matrix), "'.")
  }

  if (!nrow(distance.matrix) == ncol(distance.matrix)) {
    stop("The distance matrix has to be square of size ", obj@nb.area,
         " (number of areas, you specified a matrix of size ",
         nrow(distance.matrix),"x", ncol(distance.matrix))
  }
  if (!nrow(distance.matrix)==obj@nb.area) {
    stop("The distance matrix has to be square of size ", obj@nb.area,
         " (number of areas), you specified a square matrix of size ",
         nrow(distance.matrix))
  }
  obj@distance.matrix <- distance.matrix

  ## check multilevel.bandwidths for consistency
  obj@multilevel.bandwidths <- checkBandwidths(multilevel.bandwidths)


  ## check aggregation.function for consistency
  aggregation.function <- checkType(aggregation.function, "character",
                                    "aggregation.function")

  obj@aggregation.function <- checkAggregationFunctions(aggregation.function,
                                                        1,
                                                        obj@contextual.names
                                                        )

  ## check confidence intervals for consistency
  obj@percentiles <- checkConfidenceIntervals(confidence.intervals,
                                              obj@nb.resamples)

  ## check the individual.sample.seed for consistency
  obj@individual.sample.seed <- checkSeed(individual.sample.seed)

  ## check the contextual.sample.seed for consistency
  obj@contextual.sample.seed <- checkSeed(contextual.sample.seed)

  ## deal with kernel function
  obj@kernel <- checkKernel(kernel)

  return(obj)
}
performResampleExploreSpawML <- function(obj){

  ## prepare a weights object to compute weight.matrices on the fly
  weight.object <- new("weightsObject",
                       distance.matrix=obj@distance.matrix,
                       kernel=obj@kernel,
                       moran=FALSE)

  message("computing spatially unweighted contextual indicators")
  ## generate unweighted contextual data
  descr.obj <- new("DescribeAggregateObject",
                   contextual.data=obj@contextual.data,
                   context.id=obj@context.id,
                   contextual.names=obj@contextual.names,
                   nb.area=obj@nb.area,
                   nb.resamples=obj@nb.resamples,
                   nb.analyses=obj@nb.analyses,
                   sample.seed=obj@contextual.sample.seed,
                   individual.weight.names=getHash(obj@contextual.names,
                     obj@individual.weight.name),
                   aggregation.functions=obj@aggregation.function,
                   contextual.weight.matrices=getHash(obj@contextual.names,
                     NULL))
  unweighted.context.sampled <-
    performDescribeAggregate(descr.obj, exploratory=TRUE)
  colnames(unweighted.context.sampled) <-
    makeList(ncol(unweighted.context.sampled), obj@contextual.names[[1]])

  output.list=list()
  ## loop through multilevel.bandwidths
  message("Performing spacom computations")
  for (i in 1:length(obj@multilevel.bandwidths)) {
    multilevel.bandwidth = obj@multilevel.bandwidths[i]
    
    ## generate weight matrix

    message("computing spatial weights for bandwidth = ", multilevel.bandwidth)
    weight.matrix <- performWeights(weight.object, multilevel.bandwidth)
    message("performing spacom")
    srawe.obj <-
      new("ResampleAggregateSpawMLObject",
          contextual.sample.seed=obj@contextual.sample.seed,
          individual.sample.seed=obj@individual.sample.seed,
          individual.level.data=obj@individual.level.data,
          contextual.data=unweighted.context.sampled,
          precise.data=NULL,
          context.id=obj@context.id,
          contextual.names=obj@contextual.names,
          aggregation.names=list(),
          precise.names=list(),
          individual.weight.names=getHash(obj@contextual.names, NULL),
          aggregation.functions=list("mean"),
          contextual.weight.matrices=getHash(obj@contextual.names, weight.matrix),
          formula=obj@formula,
          nb.area=obj@nb.area,
          nb.analyses=as.integer(1),
          nb.aggregations=as.integer(0),
          nb.precise.weightings=as.integer(1),
          nb.resamples=obj@nb.resamples,
          same.survey=FALSE,
          percentiles=obj@percentiles)

    output.list[[paste("bandwidth = ", multilevel.bandwidth)]] <-
      PerformResampledMultilevel(srawe.obj, srawe=TRUE, exploratory=TRUE)
  }
  return(output.list)
}
ResampleExploreSpawML <-
  function(individual.level.data,
           contextual.name,
           contextual.data,
           context.id,
           nb.resamples,
           formula,
           distance.matrix,
           multilevel.bandwidths,
           individual.weight.name=NULL,
           aggregation.function="mean",
           confidence.intervals=c(.95),
           individual.sample.seed=NULL,
           contextual.sample.seed=NULL,
           kernel=NULL)
{
obj <-
    MakeResampleExploreSpawMLObject(individual.level.data,
                                               contextual.name,
                                               contextual.data,
                                               context.id,
                                               nb.resamples,
                                               formula,
                                               distance.matrix,
                                               multilevel.bandwidths,
                                               individual.weight.name,
                                               aggregation.function,
                                               confidence.intervals,
                                               individual.sample.seed,
                                               contextual.sample.seed,
                                               kernel=NULL)
  output.obj <- performResampleExploreSpawML(obj)
}
