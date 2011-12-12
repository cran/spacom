\name{traces_ind}
     \docType{data}
     \alias{traces_ind}
     \title{Dataset traces_ind}
     \description{
       Individual level data set, with individual level predictor and outcome variables. The traces_ind dataframe includes responses of 2254 respondents on 12 variables. 
     }
     \usage{traces_ind}
     \format{This dataframe contains the following columns:

       \code{area} - context ID variable, ranges from 1 to 80

       \code{label} - labels of areas

       \code{entity} - political entities of former Yugoslavia in 2006 (1 - Slovenia, 2 - Croatia, 3 - Bosnia, 4 - Serbia, 5 - Vojvodina, 6 - Montenegro, 7 - Kosovo, 8 - Macedonia)

       \code{victim_d} - personal experience of war victimization (0-no, 1-yes)

       \code{comb_d} - participation in combats (0-no, 1-yes)

       \code{male} - gender (0-woman, 1-man)

       \code{cg_acc} - composite results on 5 items of the Collective guilt acceptance scale, ranges from 1 to 7

       \code{cg_ass} - composite results on 5 items of the Collective guilt assignment scale, ranges from 1 to 7

       \code{age_1990} - age of respondents, in years

       \code{elem_school} - elementary school as a highest level of education (0-no, 1-yes)

       \code{high_school} - high school as a highest level of education (0-no, 1-yes)

       \code{higher_edu} - achieved higher level education (0-no, 1-yes)
}
     \source{TRACES. Full  TRACES datasets can be accessed through the Data and Research Information 
Services from the Swiss Foundation for Research in the Social Sciences (www.unil.ch/daris).
}
     \references{
Spini, D., Elcheroth, G., & Fasel, R. (2011). TRACES: Methodological and technical report. \emph{LIVES Working Papers, 4}.  

     }
     \examples{data(traces_ind)
}
