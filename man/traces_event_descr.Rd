\name{traces_event}
     \docType{data}
     \alias{traces_event}
     \title{Dataset traces_event}
     \description{
       A dataset with micro-level data to generate contextual indicator of risk of war victimization by aggregation.
The dataset consists of 243110 cases and three variables.
The cases are individuals by three-month time periods.

     }
     \usage{traces_event}
     \format{The dataframe contains the following columns:

       \code{area} - context ID variable, ranges from 1 to 80

       \code{w_all} - experience of a war event (0 - no, 1 - yes)

       \code{weight} - design weights

}
     \source{TRACES dataset. Full  TRACES datasets can be accessed through the Data and Research Information
Services from the Swiss Foundation for Research in the Social Sciences (www.unil.ch/daris).
}
     \references{
Elcheroth, G., Penic, S., Fasel, R., Giudici, F., Glaeser, S., Joye, D., Le Goff, J.-M., Morselli, D., & Spini, D. (2012). Spatially weighted context data: a new approach for modelling the impact of collective experiences.\emph{LIVES Working Papers, 19}.

Spini, D., Elcheroth, G., & Fasel, R. (2011). TRACES: Methodological and technical report. \emph{LIVES Working Papers, 4}.

     }
     \examples{data(traces_event)
}
