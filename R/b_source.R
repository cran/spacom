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


getHash <- function(keys, data) {
  indices <- list()

  if (!is(data, "list")){
    for (key in keys) {
      indices[[key]] <- 1
    }
    ret.data <- list()
    if (is.null(data)){
      data <- "RETARDEDNULL"
    }
    ret.data[[1]] <- data
    if (is.null(data)) {
      ret.data[[2]] <- 1
    }
    return(new("MultiKeyHash", indices=indices, data=ret.data))
  }
  for (i in seq(length(keys))) {
    indices[[keys[[i]]]] <- i
    if(is.null(data[[i]])) {
      data[[i]] <- "RETARDEDNULL"
    }
  }
  return (new("MultiKeyHash", indices=indices, data=data))
}

slc <- function(x, i, j, ..., drop) {
  item <- x@data[[x@indices[[i]]]]
  if(is(item, "character") && item == "RETARDEDNULL") {
    return (NULL)
  }
  return(item)
}
setMethod("[",
          signature=c("MultiKeyHash"),
          definition=slc)
setMethod("[[",
          signature=c("MultiKeyHash"),
          definition=slc)

assgnmt <- function(x, i, j, value) {
  if (is.null(value)) {
    value <- "RETARDEDNULL"
  }
  index <- length(x@data)+1
  x@data[[index]] <- value
  x@indices[[i]] <- index
  return(x)
}
setMethod(`[<-`,
          signature=c("MultiKeyHash"),
          definition=assgnmt)

setMethod(`[[<-`,
          signature=c("MultiKeyHash"),
          definition=assgnmt)



## replace the replicate function but for a list output
makeList <- function(nb.elem, object){
  ret.list <- list()
  for (i in seq (nb.elem)){
    ret.list[[i]] <- object
  }
  return(ret.list)
}

## checks whether a variable name is NULL or denotes a column of obj
checkIndividualWeights <- function(weight.name, contextual.names,
                                   context.variables) {
  if (is.null(weight.name)) {
    return(weight.name)
  }
  else if (is(weight.name,"character")) {
    if (weight.name %in% context.variables) {
      return(weight.name)
    } else {
      stop("The individual weight name '", weight.name,
           "' is not in the contextual names {", paste(contextual.variables, collapse=", "), "}")
    }
  } else {
    stop("The individual weight names have to be of type 'character' or 'NULL'",
         " or a list thereof. You gave an object of type ", class(weight.name))
  }
}

## checks whether all variable names chosen for individual weights are suitable
checkAllIndividualWeights <- function(individual.weight.names, contextual.names,
                                      context.variables,
                                      nb.analyses){
  if (!is(individual.weight.names, "list")) {
    individual.weight.names <- checkIndividualWeights(individual.weight.names,
                                                      contextual.names,
                                                      context.variables)
    return(getHash(contextual.names, individual.weight.names))
  } else {
    len.mat <- length(individual.weight.names)
    if (len.mat != nb.analyses) {
      stop("You specified ", len.mat, " individual weight variables and ",
           nb.analyses, " contextual names. ",
           "You have to specify either 1 individual weight name (which may be ",
           "NULL) or as many as you specified contextual names.")}
    tmp.list <- list()
    for (i in seq(len.mat)) {
      tmp.list[[i]] <- checkIndividualWeights(individual.weight.names[[i]],
                                              contextual.names,
                                              context.variables)}
    return(getHash(contextual.names, tmp.list))
  }
}


## checks dimensions of square matrix
checkWeightMatrix <- function(matrix, dim, name){
  nb.rows <- nrow(matrix)
  if (is.null(matrix)){}
  else if (is(matrix, "matrix")){
    if (nb.rows!=ncol(matrix) || nb.rows!=dim){
      stop(cat("The weight matrix given for", name, "is", nb.rows, "x",
               ncol(matrix), "but should be square", dim, "x", dim, '\n'))}
  } else {
    stop(cat("The weight matrices have to be of type 'matrix' or NULL.",
             "You gave an object of type", class(matrix)))
  }
  return(matrix)
}

## checks whether all all weight matrices are suitable
checkAllWeightMatrices <- function(contextual.weight.matrices,
                                   contextual.names,
                                   nb.area, nb.analyses){
  if (!is(contextual.weight.matrices,"list")) {
    contextual.weight.matrices <- checkWeightMatrix(contextual.weight.matrices,
                                                    nb.area,
                                                    contextual.names[[1]])
    return(getHash(contextual.names, contextual.weight.matrices))
  } else {
    len.mat <- length(contextual.weight.matrices)
    tmp.list <- list()
    if (len.mat != nb.analyses){
      stop("You specified ", len.mat, " weight matrices and ",
           nb.analyses, " contextual names. ",
           "You have to specify either 1 weight (which may be NULL) or as many ",
           "as you specified contextual names")}
    for (i in seq(len.mat)){
      tmp.list[[i]] <- checkWeightMatrix(contextual.weight.matrices[[i]],
                                         nb.area,
                                         contextual.names[[i]])
    }
    length(tmp.list) <- len.mat
    return(getHash(contextual.names, tmp.list))
  }
}

## makes sure that the contextual names are in a list
checkContextualNames <- function(obj, contextual.names, contextual.variables){

  if (!is(contextual.names, "list")) {
    contextual.names <- as.list(contextual.names)}
  for (name in as.list(contextual.names)) {
    if (!name %in% contextual.variables) {
      stop("The variable '", name, "' is not in the context data")}}
  ## obj@contextual.names <- contextual.names

  obj@nb.aggregations <- as.integer(0)
  obj@nb.precise.weightings <- as.integer(0)

  ## make an array of individualised contextual names to avoid name clashes
  counter <- 0
  coder <- function(x){
    counter <<- counter + 1
    return(paste(x,sprintf("%05d",counter), sep=""))
  }
  obj@contextual.names <- lapply(contextual.names, coder)
  obj@nb.analyses <- length(obj@contextual.names)

  ## split the contextual names into names for aggregation and names for
  ## precise contexts
  if (length(obj@contextual.names)>0) {
    for (i in 1:length(obj@contextual.names)) {
      name <- contextual.names[[i]]
      if (name %in% names(GetContextMethod(obj))) {
        obj@nb.aggregations <- obj@nb.aggregations + as.integer(1)
        obj@aggregation.names[[obj@nb.aggregations]] <- obj@contextual.names[[i]]
      } else if (!is.null(obj@precise.data)) {
        if (name %in% names(obj@precise.data)) {
          obj@nb.precise.weightings <- obj@nb.precise.weightings + as.integer(1)
          obj@precise.names[[obj@nb.precise.weightings]] <- obj@contextual.names[[i]]
        }
      }
    }
  }
  return(obj)
}

## makes sure context ids exist in the right dataframes
checkContextId <- function(context.id, individual.level.data.names,
                           contextual.data.names=NULL, precise.data.names=NULL,
                           skip.individual.check=FALSE){
  if (!skip.individual.check && !context.id %in% individual.level.data.names){
    stop("The value given for context.id '", context.id, "' is not in the ",
         "names of individual.level.data")
  }
  if (!is.null(contextual.data.names)){
    if (!context.id %in% contextual.data.names){
      stop("The value given for context.id '", context.id, "' is not in the",
           "names of contextual.data")
    }
  }
  if (!is.null(precise.data.names)){
    if(!context.id %in% precise.data.names){
      stop("The value given for context.id '", context.id, "' is not in the",
           "names of precise.data")
    }
  }
  return(context.id)
}

## makes sure a object is either one of defined types
checkType <- function(obj, classes, name){
  cls <- class(obj)
  if (!cls %in% classes){
    stop(name, " has to be of either of the types {",
         paste(classes, collapse=", "), "}. You gave an object of class ", cls)
  }
  return(obj)
}

## makes sure the aggregation funciton list is suitable
checkAggregationFunctions <- function(aggregation.functions,
                                      nb.analyses,
                                      aggregation.names){
  aggregation.functions <- checkType(aggregation.functions,
                                     c("character", "list"),
                                     "aggregation.functions")
  if (!class(aggregation.functions) == "list") {
    aggregation.functions <- as.list(replicate(nb.analyses, aggregation.functions))
  } else {
    len.mat <- length(aggregation.functions)
    if (len.mat != 1 && len.mat != nb.analyses) {
      stop("You specified ", len.mat, " aggregation functions and ",
           nb.analyses, " contextual names. ",
           "You have to specify either 1 function (which by default is 'mean') or as many ",
           "as you specified contextual names")}
    for (i in seq(len.mat)){
      checkType(aggregation.functions[[i]], c("character"),
                paste("aggregation.functions[[", i, "]]", sep=""))
    }
    
    if (len.mat == 1){
      aggregation.functions <- replicate(nb.analyses, aggregation.functions)
    }
  }
  names(aggregation.functions) <- aggregation.names
  return(aggregation.functions)
}

## makes sure that the object given as formula can be coerced into a formula
checkFormula <- function(formula){
  tryCatch(formula <- as.formula(formula),
           error=function(er) {
             stop(cat('what you gave as formula "', formula,
                      '" could not be coerced into a formula. Please give ',
                      'a character string or a formula\n'))})
  return(formula)
}

## makes sure that the number of resamples is positive and integer
checkNbResamples <- function(nb.resamples) {
  if (!is(nb.resamples, "numeric")) {
    stop("nb.resamples has to be a numeric (integer)  value. You specified ",
         "an object of class ", class(nb.resamples))
  }
  nb.resamples <- as.integer(nb.resamples)
  if (nb.resamples < 1) {
    stop("nb.resamples has to be larger than 1. You specified ",
         nb.resamples)
  }
  return(nb.resamples)
}

## makes sure that the confidence intervals are correctly defined and contain
## the median
checkConfidenceIntervals <- function(confidence.intervals, nb.resamples){
  tryCatch(percentiles <- as.numeric(confidence.intervals),
           error=function(er){
             stop("what you gave as percentiles '", confidence.intervals,
                  "' could not ",
                  "be coerced into numeric. Please give a numeric vector of ",
                  "percentiles")})
  len.confidence.intervals <- length(confidence.intervals)
  len.percentiles = 2*len.confidence.intervals + 1
  percentiles <- rep(0, len.percentiles)
  names <- rep('', len.percentiles)
  ## insure that the median is part of the percentiles
  for (i in seq(len.confidence.intervals)) {
    if (confidence.intervals[i] >= 1 || confidence.intervals[i] <= 0) {
      stop("Confidence interval ", i, " = ", confidence.intervals[i], " ",
           "is outside the interval I = {x|0<x<1}. Make sure all confidence ",
           "intervals are within I")}
    exclude <- 1 - confidence.intervals[i]
    if (exclude*nb.resamples < 50){
      warning("Warning, the confidence interval ", i, " = ",
              confidence.intervals[i], " excludes only ",
              ceiling(exclude*nb.resamples), " individuals. It may be ",
              "unreliable.")
    }
    percentiles[2*i-1] <- .5*exclude
    names[2*i-1] <- paste("CI_", confidence.intervals[i], "_lower", sep='')
    percentiles[2*i] <- 1-percentiles[2*i-1]
    names[2*i] <- paste("CI_", confidence.intervals[i], "_upper", sep='')
  }
  percentiles[len.percentiles] <- .5
  names[len.percentiles] <- "median"
  names(percentiles) <- names
  return(percentiles)
}

## makes sure the kernel function is a function (yeah, I heard it too)
checkKernel <- function(kernel) {
  if (is.null(kernel)) {
    kernel <- function(distance.matrix, h) {
      return(.5^((distance.matrix/h)^2))
    }
  } else if (!is(kernel, "function")) {
    stop("The kernel function has to be of class 'function'. You specified a ",
         "object of class '", class(kernel), "'")
  }

  return(kernel)
}


##makes sure that specified bandwidths are suitable
checkBandwidths <- function(bandwidths){
  if (!is(bandwidths, "numeric")) {
    stop("The bandwidths have to be a vector of numeric. You specified an ",
         "object of class '", class(bandwidths), "'.")
  }
  if (any(bandwidths<0)) {
    stop("only positive bandwidths are allowed. You specified ",
         paste(bandwidths, collapse=", "))
  }
  return(bandwidths)
}
  

## makes sure a sample seed is a valid random seed or a list of samples
checkSeed <- function(sample.seed, nb.resamples){
  
  runif(1) # single call to runif to assure the existence of .Random.seed
  current.seed <- .Random.seed
  if (is(sample.seed, "numeric")){
    if (!length(sample.seed) ==1){
      tryCatch(.Random.seed <- sample.seed, runif(1),
               error=function(er) {
                 stop("You specified a numeric of length ", length(sample.seed), " for ",
                      "the sample seed :", sample.seed, "\n Please specify one of the ",
                      "following instead:\na) a single numeric\nb) a .Random.seed",
                      "\nc) NULL\n",
                      "If you do not understand this error message, you probably want",
                      "to specify either NULL (if you do not to reproduce the same",
                      "samples) or a single digit number (for reproducibility)")})
      return(sample.seed)
    } else {
      set.seed(sample.seed)
      return(.Random.seed)
    }
  } else if (is.null(sample.seed)) {
    return(.Random.seed)
  } else {
    stop("You specified a ", class(sample.seed), " of length ",
         length(sample.seed), " for ",
         "the sample seed :", sample.seed, "\n Please specify one of the ",
         "following instead:\na) a single numeric\nb) a .Random.seed",
         "\nc) NULL\n",
         "If you do not understand this error message, you probably want",
         "to specify either NULL (if you do not to reproduce the same",
         "samples) or a single digit number (for reproducibility)")
  }
}

## if there are contextual weights, apply them and row standardize
performSpatialWeighting <- function(spatial.weight, aggregated.context) {
  if (is.null(spatial.weight)) {
    return(aggregated.context)}
  return (spatial.weight %*% aggregated.context/rowSums(spatial.weight))
}

performAllSpatialWeighting <- function(contextual.names,
                                       context.id,
                                       nb.area,
                                       nb.analyses,
                                       contextual.data,
                                       contextual.weight.matrices){
  return(performAggregation(contextual.names=contextual.names,
                            context.id=context.id,
                            nb.area=nb.area,
                            nb.analyses=nb.analyses,
                            individual.weight.names=NULL,
                            contextual.data=contextual.data,
                            aggregation.function=NULL,
                            contextual.weight.matrices=contextual.weight.matrices))
}

performAggregation <- function(contextual.names,
                               context.id,
                               nb.area,
                               nb.analyses,
                               individual.weight.names,
                               contextual.data,
                               aggregation.functions,
                               contextual.weight.matrices,
                               formula.str=NULL){
  ## prepare a zero-filled named list to keep count
  count.list <- list()
  decoded.names <- lapply(contextual.names, function(x)
                          {return(substr(x, 1, nchar(x)-5))})
  for(name in decoded.names) {
    count.list[[name]] <- 0}
  
  ## prepare the merge dataframe into which the new weighted contextual data
  ## are loaded
  merge.data <- data.frame(1:nb.area)
  names(merge.data) <- context.id
  ## loop through the contextual variables to be analysed

  if (nb.analyses > 0) {
    for (i in 1:nb.analyses) {
      coded.name <- contextual.names[[i]]
      len <- nchar(coded.name)
      name <- substr(coded.name, 1, len-5)

      if (is.null(aggregation.functions)){
        aggregated.context <-
          performSpatialWeighting(contextual.weight.matrices[[coded.name]],
                                  contextual.data[[name]])
      } else {
        if (is.null(individual.weight.names[[coded.name]])) {
          ## if there is only one argument to the aggregation function supplied by
          ## the user, aggregate the data directly
          aggregated.context <-
            matrix(aggregate(contextual.data[[name]],
                             by=list(contextual.data[[context.id]]),
                             FUN=aggregation.functions[[coded.name]])[,2],
                   ncol=1)}
        else {
          ## else, it has to be done explicitely in three steps
          ## a) extract contextual data sorted by area for aggregation
          ## aggr.variable <- contextual.data[[contextual.names[[i]]]]
          ## aggr.weights <-  contextual.data[[individual.weight.names[[i]]]]
          columns.extract <- c(name,
                               individual.weight.names[[coded.name]])
          context.and.weight <- subset(contextual.data, select=columns.extract)
          names(context.and.weight) <- c('x', 'w')

          ## b) split the contextual data by area
          split.context <- split(context.and.weight,
                                 contextual.data[[context.id]])
          
          ## c) aggregate the areas individually
          aggregated.context <- matrix(sapply(X=split.context,
                                              FUN=function(x,y) {do.call(y,x)},
                                              aggregation.functions[[coded.name]]),
                                       ncol=1)
        }
        aggregated.context <-
          performSpatialWeighting(contextual.weight.matrices[[coded.name]],
                                  aggregated.context)
      }

      
      ## store the aggregated context in merge.data for later merge with
      ## the rest, compute the appropriate renaming of the contextual variables
      count.list[[name]] <- count.list[[name]] + 1
      new.name <- paste(name, ".", count.list[[name]], sep='')
      
      merge.data[[new.name]] <- aggregated.context
      ## add the contextual variable to the formula
      if (!is.null(formula.str)) {
        formula.str <- paste(formula.str, " + ", new.name, sep=)
      }
    }
  } else {
    merge.data <- NULL
  }

  if (is.null(formula.str)) {
    return(merge.data)
  } else {
    return(list(merge.data, formula.str))
  }
}

## resamples contextual data at the individual level as used for stratified
## resampling used in descrw and in srawe
resample <- function(area.list){
  sampler <- function(group){
    size <- length(group)
    row.slice <- sample(1:size, size=size, replace=TRUE)
    return(group[row.slice])
  }
  split.contexts = split(1:length(area.list), area.list)
  resampled.groups <- lapply(split.contexts, sampler)
  ## cont.dat$lino=1:nrow(cont.dat)
  ## return(do.call(rbind,
  ##                lapply(split(cont.dat, area.list),
  ##                       function(x) x[sample(1:nrow(x), size=nrow(x), replace=TRUE),])))
  return(do.call(c, resampled.groups))
}

## a stateful iterator which yields either presampled samples or creates them
## on the fly
sampleIterator <- function(sample.seed, area.list = NULL,
                           nb.resamples = NULL){
  if (is.null(area.list)) {
    return(iter(sample.seed))
  } else {
    .Random.seed <<- sample.seed
    ran.seed <- sample.seed
    index <- 1
    nextEl <- function() {
      if (index <= nb.resamples) {
        index <<- index + 1
      } else {
        stop("StopIteration")
      }
      current.seed <- .Random.seed
      .Random.seed <<- ran.seed
      sample <- resample(area.list)
      ran.seed <<- .Random.seed
      .Random.seed <<- current.seed
      return(sample)
    }
    obj <- list(nextElem=nextEl)
    class(obj) <- c("random.sample.iterator", "abstractiter", "iter")
    return(obj)
  }
}
