\name{d_geo}
     \docType{data}
     \alias{d_geo}
     \title{A distance matrix for TRACES areas}
     \description{

A matrix of dimension \eqn{80\times 80}{80 by 80}, containing geographical distances (in km) \eqn{d_{ij}}{d_ij} between centers of TRACES areas \eqn{i} and \eqn{j}.
     }
     \usage{d_geo}
     \source{Elcheroth, G., Penic, S., Fasel, R., Giudici, F., Glaeser, S., Joye, D., Le Goff, J.-M., Morselli, D., & Spini, D. (2012). Spatially weighted context data: a new approach for modelling the impact of collective experiences.\emph{LIVES Working Papers, 19}. 
}
     \examples{data(d_geo)
geow.50 <- WeightMatrix(distance.matrix=d_geo, bandwidth=50)
}
