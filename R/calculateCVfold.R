#-----------------------------------------------------------------------------------
# script to train the classifier in each CV fold, including batch adjustment and 
# feature selection.
#
# Note, in this example we reduced the number of features to 20k probes by sd filtering before applying the random forest for 
# feature selection. To perform feature selection as described in the paper remove line 25,
# this will increase the computation time significantly.
# 
# Martin Sill
# m.sill@dkfz.de                                                                  
# 
# 2018-03-14 UTC
#------------------------------------------------------------------------------------                   
#options(max.print = 1000)
#options(stringsAsFactors = FALSE)
#options(scipen = 999)

calcultateCVfold <- function(Mset,y,batch,fold,p,cores,ntrees){
  
  
  message("adjusting for batch effects ...",Sys.time())
  badj <- batchadjust(Mset,batch,fold)
  
  # sd pre filtering to 20k probes, to speed up the example
  badj$betas.train <- badj$betas.train[,order(-apply(badj$betas.train,2,sd))[1:20000]]
  
  message("performing variable selection ...",Sys.time())
  message("cores: ",cores)
  message("ntrees: ",ntrees)  
  message("n: ",nrow(badj$betas.train))
  message("p: ",ncol(badj$betas.train))  
  
  rf.varsel <- rfp(badj$betas.train,
                   y=y[fold$train],
                   mc=cores,
                   ntree=ntrees,
                   sampsize=rep(min(table(y[fold$train])),length(table(y[fold$train]))),
                   importance=TRUE)
  
  # get permutation variable importance
  imp.meandecrease <- rf.varsel$importance[,dim(rf.varsel$importance)[2]-1]
  
  # reduce data matrix
  or <- order(imp.meandecrease,decreasing=T)
  
  message("training classifier ...",Sys.time())
  message("cores: ",cores)
  message("ntrees: ",ntrees)  
  message("n: ",nrow(badj$betas.train))
  message("p: ",p)  
  
  rf <- rfp(badj$betas.train[,or[1:p]],y[fold$train],
            sampsize=rep(min(table(y[fold$train])),length(table(y[fold$train])))
            ,mc=cores,ntree=ntrees,importance=TRUE) 
  
  message("predicting test set ...",Sys.time())
  
  rf.scores <- predict(rf,badj$betas.test[,match(rownames(rf$importance),
                                                 colnames(badj$betas.test))],type="prob")
  
  err <- sum(colnames(rf.scores)[apply(rf.scores,1,which.max)]!=y[fold$test])/length(fold$test)
  message("misclassification error: ",err)
  
  return(rf.scores)
}
