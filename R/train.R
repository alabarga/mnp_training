#-----------------------------------------------------------------------------------
# parallelized version of random Forest
#                                                                     
# Martin Sill
# m.sill@dkfz.de                                                                  
# 
# 2018-03-14 UTC
#------------------------------------------------------------------------------------                   
#options(max.print = 1000)
#options(stringsAsFactors = FALSE)
#options(scipen = 999)

library(randomForest)
library(parallel)

rfp <- function(xx, ..., ntree = ntree, mc = mc, seed = NULL)
{
  if(!is.null(seed)) set.seed(seed, "L'Ecuyer")
  rfwrap <- function(ntree, xx, ...) randomForest::randomForest(x=xx,ntree=ntree,norm.votes=FALSE, ...)
  rfpar <- mclapply(rep(ceiling(ntree/mc),mc),mc.cores=mc, rfwrap, xx=xx, ...)
  do.call(randomForest::combine, rfpar)
}