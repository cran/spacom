\name{homog_census}
     \docType{data}
     \alias{homog_census}
     \title{Dataset homog_census}
     \description{
A dataset with an indicator of ethnic homogeneity in TRACES areas of former Yugoslavia, calculated on the basis of national censuses and official estimates from 2000 to 2006. It should be used as a precise contextual indicator in multilevel analyses. 

     }
     \usage{homog_census}
     \format{The dataframe contains the following columns:
       
       \code{area} - context ID variable, ranges from 1 to 80
       
       \code{Homog_00} - proportion of largest ethnic group in the area

}
     \source{Our own calculation based on national censuses and official estimates from 2000 to 2006.
}
  
     \examples{data(homog_census)
}
