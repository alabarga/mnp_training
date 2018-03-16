#--------------------------------------------------------------------
# mutliclass brier
#
#
# Martin Sill
# m.sill@dkfz.de                                                                  
# 
# 2018-03-14 UTC
#--------------------------------------------------------------------        

brier <- function(scores,y){
  ot <- matrix(0,nrow=nrow(scores),ncol=ncol(scores))
  arr.ind <- cbind(1:nrow(scores),match(y,colnames(scores)))
  ot[arr.ind] <- 1
  sum((scores - ot)^2)/nrow(scores)
}
