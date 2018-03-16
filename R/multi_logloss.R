#--------------------------------------------------------------------
# multiclass logloss
#
#
# Martin Sill
# m.sill@dkfz.de                                                                  
# 
# 2018-03-14 UTC
#--------------------------------------------------------------------        
mlogloss <- function(scores,y){
  N <- nrow(scores)
  y_true <- matrix(0,nrow=nrow(scores),ncol=ncol(scores))
  arr.ind <- cbind(1:nrow(scores),match(y,colnames(scores)))
  y_true[arr.ind] <- 1
  eps <- 1e-15
  scores <- pmax(pmin(scores, 1 - eps), eps)
  (-1 / N) * sum(y_true * log(scores))
}
