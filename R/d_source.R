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


## define a method which yields the contextual data, if it exists, or the
## individual level data in the opposite case
GetContextMethod <- function(object) {
  if (is.null(object@ contextual.data)) {
    object@ individual.level.data}
  else {
    object@ contextual.data}}

GetContext <- function(object){}
setGeneric("GetContext")
setMethod("GetContext", signature=c("SpawMLObject"), definition=GetContextMethod)


## provides the SpawMLObject needed for the analysis and 
## performs all consistency checks
MakeSpawMLObject <- function(individual.level.data,
                                         contextual.names,
                                         context.id,
                                         formula,
                                         contextual.data=NULL,
                                         precise.data=NULL,
                                         contextual.weight.matrices=NULL,
                                         individual.weight.names=NULL,
                                         aggregation.functions="mean",
                                         obj=NULL) {
  
  ## create an empty SpawMLObject to fill
  if (is.null(obj)) {
    obj <- new("SpawMLObject")
  }
  
  ## make sure individual.level.data is a data.frame
  obj@individual.level.data <- checkType(individual.level.data, c("data.frame"),
                                         "individual.level.data")

  ## make sure contextual.data is either NULL or a data.frame
  obj@contextual.data <- checkType(contextual.data, c("NULL", "data.frame"),
                                   "contextual.data")
  ## make sure precise.data is either NULL or a data.frame
  obj@precise.data <- checkType(precise.data, c("NULL", "data.frame"),
                                "precise.data")

  ## make sure context.id is a name in contextual.data
  obj@context.id <- checkContextId(context.id, names(obj@individual.level.data),
                                   names(GetContext(obj)),
                                   names(obj@precise.data))

  ## extract number of upper level units
  obj@nb.area <- length(levels(as.factor(GetContext(obj)[[context.id]])))

  ## make sure the contextual names are in the data
  ref.names <- names(GetContext(obj))
  if (!is.null(obj@precise.data)) {
    ref.names <- c(ref.names, names(obj@precise.data))
  }
  obj <- checkContextualNames(obj, contextual.names, ref.names)
  
  ## check that every name appears only in either aggregation.names or in
  ## precise.names
  for (name in obj@aggregation.names) {
    if (name %in% obj@precise.names) {
      stop("The contextual name '", name, "' is present in both the precise ",
           "data and the contextual data for aggregation, but can appear only ",
           "in one of them")
    }
  }
  
  
  ## make sure the weight matrix list is suitable
  obj@contextual.weight.matrices <-
    checkAllWeightMatrices(contextual.weight.matrices, obj@contextual.names,
                           obj@nb.area, obj@nb.analyses)

  ## make sure the individual.weight.names are suitable
  obj@individual.weight.names <-
    checkAllIndividualWeights(individual.weight.names, obj@aggregation.names,
                              ref.names,
                              obj@nb.aggregations)

  ## make sure the aggregation function list is suitable
  obj@aggregation.functions <- checkAggregationFunctions(aggregation.functions,
                                                         obj@nb.aggregations,
                                                         obj@aggregation.names)
  ## make sure the formula isn't a string anymore
  obj@formula <- checkFormula(formula)

  
  return(obj)
}

summary.SpawMLOutput <- function(object,...) {
  print(summary(object@lme))
  message("\nStandardised fixed effects:")
  print(object@beta)}

setGeneric("summary")
setMethod("summary",
          signature=c("SpawMLOutput"),
          definition=summary.SpawMLOutput)

setMethod("fixef",
          signature=c("SpawMLOutput"),
          definition=function(object, ...) {lme4::fixef(object@lme, ...)})
setMethod("ranef",
          signature=c("SpawMLOutput"),
          definition=function(object, ...) {lme4::ranef(object@lme, ...)})
setMethod("VarCorr",
          signature=c("SpawMLOutput"),
          definition=function(x, ...) {lme4::VarCorr(x@lme, ...)})
setMethod("AIC",
          signature=c("SpawMLOutput"),
          definition=function(object, ..., k=2) {lme4::AIC(object@lme, ..., k)})
setMethod("BIC",
          signature=c("SpawMLOutput"),
          definition=function(object, ..., k=2) {lme4::BIC(object@lme, ..., k)})

printSmlObject <- function(x, is.print=TRUE) {
  print(x@lme)
  cat("\nStandardised fixed effects:\n")
  print(x@beta)
}

setGeneric("print")
setMethod("print", signature="SpawMLOutput",
          definition=printSmlObject)
setMethod("show", signature="SpawMLOutput",
          definition=function(object) {printSmlObject(object, is.print=FALSE)})

PerformSpawML <- function(obj) {
  ## create an output object 
  output.obj <- new("SpawMLOutput")

  ## chop up the formula in order to manipulate it
  formula.str<- paste(as.character(obj@ formula[2]), '~',
                      as.character(obj@ formula[3]))

  ## perform aggregation of contextual data for aggregation
  l <- performAggregation(obj@aggregation.names,
                          obj@context.id,
                          obj@nb.area,
                          obj@nb.aggregations,
                          obj@individual.weight.names,
                          GetContext(obj),
                          obj@aggregation.functions,
                          obj@contextual.weight.matrices,
                          formula.str)
  
  aggretation.data <- l[[1]]
  
  formula.str <- l[[2]]

  ## perform weighting of exact data (same function, but different behaviour
  l <- performAggregation(contextual.names=obj@precise.names,
                          context.id=obj@context.id,
                          nb.area=obj@nb.area,
                          nb.analyses=obj@nb.precise.weightings,
                          individual.weight.names=NULL,
                          contextual.data=obj@precise.data,
                          aggregation.functions=NULL,
                          contextual.weight.matrices=obj@contextual.weight.matrices,
                          formula.str=formula.str)
  precise.data <- l[[1]]
  formula.str <- l[[2]]

  if (!is.null(aggretation.data) && !is.null(precise.data)) {
    merge.data <- merge(aggretation.data,
                        precise.data,
                        by=obj@context.id)
  } else if (is.null(aggretation.data) && !is.null(precise.data)) {
    merge.data <- precise.data
  } else if (!is.null(aggretation.data) && is.null(precise.data)) {
    merge.data <- aggretation.data
  } else {
    merge.data <- NULL
  }

  if (!is.null(merge.data)) {
    merged.data <- merge(obj@ individual.level.data,
                         merge.data, by=obj@ context.id)
  } else {
    merged.data <- obj@individual.level.data
  }
  fixed.effect.formula = as.formula(formula.str)

  output.obj@lme <- lmer(formula=fixed.effect.formula,
                         data=merged.data)

  ## compute the standardised coefficients for contextual data
  coefficient.names <- names(output.obj@lme@fixef)[-1]
  nb.coefficients <- length(coefficient.names)
  output.obj@beta <- numeric(nb.coefficients)
  names(output.obj@beta) <- coefficient.names
  cont.names <- names(merge.data)
  sd.outcome <- sd(obj@individual.level.data[[as.character(obj@formula[2])]])
  for (name in coefficient.names) {
    if (name %in% cont.names) {
      sd.explanatory <- sd(as.numeric(merge.data[[name]]))
    } else {
      sd.explanatory <- sd(as.numeric(obj@individual.level.data[[name]]))
    }
    output.obj@beta[[name]] <- output.obj@lme@fixef[[name]]*sd.explanatory/sd.outcome
  }
  return(output.obj)
}

SpawML <- function(individual.level.data,
                               contextual.names,
                               context.id,
                               formula,
                               contextual.data=NULL,
                               precise.data=NULL,
                               contextual.weight.matrices=NULL,
                               individual.weight.names=NULL,
                               aggregation.functions="mean") {
  ## build the SpawMLobject
  obj <- MakeSpawMLObject(individual.level.data,
                                      contextual.names,
                                      context.id,
                                      formula,
                                      contextual.data,
                                      precise.data,
                                      contextual.weight.matrices,
                                      individual.weight.names,
                                      aggregation.functions)
  model <- PerformSpawML(obj)

}

