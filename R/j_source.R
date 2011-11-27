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

MakeExploreSpawMLObject <-
  function(individual.level.data,
           contextual.name,
           contextual.data,
           precise.data,
           context.id,
           formula,
           distance.matrix,
           multilevel.bandwidths,
           individual.weight.names=NULL,
           aggregation.function="mean",
           kernel=NULL)
{
  obj <- new("ExploreSpawML")

  ## make sure individual.level.data is a data.frame
  obj@individual.level.data <- checkType(individual.level.data, c("data.frame"),
                                         "individual.level.data")

  ## make sure contextual.data is either NULL or a data.frame
  obj@contextual.data <- checkType(obj=contextual.data,
                                   classes=c("NULL", "data.frame"),
                                   name="contextual.data")
  obj@precise.data <- checkType(obj=precise.data,
                                   classes=c("NULL", "data.frame"),
                                   name="precise.data")
  contextual.variables <- names(GetContextMethod(obj))
  if (!is.null(precise.data)) {
    contextual.variables <- c(contextual.variables,
                              names(precise.data))
  }
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

  ## deal with kernel function
  obj@kernel <- checkKernel(kernel)

  return(obj)
}

performExploreSpawML <- function(obj){
  ## prepare a weights object to compute weight.matrices on the fly
  weight.object <- new("weightsObject",
                       distance.matrix=obj@distance.matrix,
                       kernel=obj@kernel,
                       moran=FALSE)
  coded.name <- obj@contextual.names[[1]]
  name <- substr(coded.name, 1, nchar(coded.name)-5)

  message("computing spatially unweighted contextual indicators")

  ## generate unweighted contextual data
  if (length(obj@aggregation.names) == 0){
    merge.data <- obj@precise.data
  } else {
    if (is.null(obj@individual.weight.name)) {
      aggregated.context <-
        matrix(aggregate(obj@contextual.data[[name]],
                         list(GetContextMethod(obj)[[obj@context.id]]),
                         obj@aggregation.function[[1]])[,2],
               ncol=1)
    } else {
      columns.extract <- c(name,
                           obj@individual.weight.name)
      context.and.weight <- subset(obj@contextual.data, select=columns.extract)
      names(context.and.weight) <- c('x', 'w')
      
      split.context <- split(context.and.weight,
                             obj@contextual.data[[obj@context.id]])
      
      aggregated.context <- matrix(sapply(X=split.context,
                                          FUN=function(x,y) {do.call(y,x)},
                                          obj@aggregation.function[[1]]),
                                   ncol=1)
      
    }
    merge.data <- data.frame(1:obj@nb.area)
    names(merge.data) <- obj@context.id
    merge.data[[name]] <- aggregated.context
    if(!is.null(obj@precise.data)) {
      merge.data <- merge(merge.data, obj@precise.data, by=obj@context.id)
    }
  }

  output.list=list()

  ## loop through multilevel.bandwidths
  message("performing spacom computations")
  for (i in 1:length(obj@multilevel.bandwidths)) {
    multilevel.bandwidth = obj@multilevel.bandwidths[i]
    
    ## generate weight matrix

    message("computing spatial weights for bandwidth = ", multilevel.bandwidth)
    weight.matrix <- performWeights(weight.object, multilevel.bandwidth)
    message("performing spacom")
    sml.obj <-
      new("ResampleAggregateSpawMLObject",
          individual.level.data=obj@individual.level.data,
          contextual.data=NULL,
          precise.data=merge.data,
          context.id=obj@context.id,
          contextual.names=obj@contextual.names,
          aggregation.names=list(),
          precise.names=obj@contextual.names,
          individual.weight.names=getHash(obj@contextual.names, NULL),
          aggregation.functions=list("mean"),
          contextual.weight.matrices=getHash(obj@contextual.names, weight.matrix),
          formula=obj@formula,
          nb.area=obj@nb.area,
          nb.analyses=as.integer(1),
          nb.aggregations=as.integer(0),
          nb.precise.weightings=as.integer(1))
    output.list[[paste("bandwidth = ", multilevel.bandwidth)]] <-
      PerformSpawML(sml.obj)
  }
  return(output.list)
}

ExploreSpawML <-
  function(individual.level.data,
           contextual.name,
           contextual.data,
           context.id,
           formula,
           distance.matrix,
           multilevel.bandwidths,
           precise.data=NULL,
           individual.weight.names=NULL,
           aggregation.function="mean",
           kernel=NULL)
{
obj <-
  MakeExploreSpawMLObject(individual.level.data,
                                contextual.name,
                                contextual.data,
                                precise.data,
                                context.id,
                                formula,
                                distance.matrix,
                                multilevel.bandwidths,
                                individual.weight.names,
                                aggregation.function,
                                kernel)
  output.obj <- performExploreSpawML(obj)
}
