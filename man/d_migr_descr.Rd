\name{d_migr}
     \docType{data}
     \alias{d_migr}
     \title{A distance matrix for TRACES areas}
     \description{
     A matrix of dimension \eqn{80\times 80}{80 by 80}, containing distances \eqn{d_{ij}}{d_ij} between centers of TRACES areas \eqn{i} and \eqn{j} defined as a lack of
contact between two populations.

     }
     \usage{d_migr}
     \source{Elcheroth, G., Penic, S., Fasel, R., Giudici, F., Glaeser, S., Joye, D., Le Goff, J.-M., Morselli, D., & Spini, D. (forthcoming). Spatially weighted context data: a new approach for modelling the impact of collective experiences.\emph{LIVES Working Papers, 4}. 
}
     \examples{data(d_migr)
w_migr <- WeightMatrix(distance.matrix=d_migr, bandwidth=25)
}
