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


## provides a ResampleAggregateSpawMLObject needed for the analysis and 
## performs all consistency checks
MakeResampleAggregateSpawMLObject <-function(individual.level.data,
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

                                                     contextual.sample.seed=NULL) {

  ## create an empty ResampleAggregateSpawMLObject to fill, using the
  ## parent class SpawMLObject
  obj <-
    MakeResampleExactSpawMLObject(
      individual.level.data=individual.level.data,
      contextual.names=contextual.names,
      context.id=context.id,
      formula=formula,
      contextual.data=contextual.data,
      precise.data=precise.data,
      contextual.weight.matrices=contextual.weight.matrices,
      individual.weight.names=individual.weight.names,
      aggregation.functions=aggregation.functions,
      confidence.intervals,
      nb.resamples=nb.resamples,
      individual.sample.seed=individual.sample.seed,

      new("ResampleAggregateSpawMLObject"))
  
  ## make sure the sample seed is ok
  obj@contextual.sample.seed <- checkSeed(contextual.sample.seed)

  ## internal use of the package creator, neglect this
  obj@same.survey=FALSE
  return(obj)
}


ResampleAggregateSpawML <-function(
           individual.level.data,
           contextual.names,
           context.id,
           formula,
           contextual.data=NULL,
           precise.data=NULL,
           contextual.weight.matrices=NULL,
           individual.weight.names=NULL,
           aggregation.functions="mean",

           confidence.intervals=c(.95),
           nb.resamples=1000,
           individual.sample.seed=NULL,
           contextual.sample.seed=NULL) {
  obj <-
    MakeResampleAggregateSpawMLObject(individual.level.data,
                                              contextual.names,
                                              context.id,
                                              formula,
                                              contextual.data,
                                              precise.data,
                                              contextual.weight.matrices,
                                              individual.weight.names,
                                              aggregation.functions,

                                              confidence.intervals,
                                              nb.resamples,
                                              individual.sample.seed,
                                              contextual.sample.seed)
  model <- PerformResampledMultilevel(obj, srawe=TRUE)
}
