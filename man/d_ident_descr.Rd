\name{d_ident}
     \docType{data}
     \alias{d_ident}
     \title{A distance matrix for TRACES areas}
     \description{
      A matrix of dimension \eqn{80\times 80}{80 by 80}, containing distances \eqn{d_{ij}}{d_ij} between centers of TRACES areas \eqn{i} and \eqn{j} defined as a lack of
common identification.

     }
     \usage{d_ident}
     \source{Elcheroth, G., Penic, S., Fasel, R., Giudici, F., Glaeser, S., Joye, D., Le Goff, J.-M., Morselli, D., & Spini, D. (2012). Spatially weighted context data: a new approach for modelling the impact of collective experiences.\emph{LIVES Working Papers, 19}. 
}
     \examples{data(d_ident)
w.id <- WeightMatrix(distance.matrix=d_ident, bandwidth=2)
}
