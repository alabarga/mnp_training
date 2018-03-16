#-----------------------------------------------------------------------------------
# training Random Forest classifier
#                                                                     
# Martin Sill
# m.sill@dkfz.de                                                                  
# 
# Note, in this example we reduce the number of features to 20k probes by sd filtering before applying the random forest for 
# feature selection. To perform feature selection as described in the paper remove line 34,
# this will increase the computation time significantly.
#
# 2018-03-14 UTC
#------------------------------------------------------------------------------------                   
#options(max.print = 1000)
#options(stringsAsFactors = FALSE)
#options(scipen = 999)
rm(list=ls())

library(randomForest)
library(parallel)

ntrees <- 500  # 10000 in the paper, here 500 to speed up the example
cores <- 4
seed <- 180314
p <- 10000   

message("loading preprocessed data ...",Sys.time())
load("./results/betas.ba.RData")

message("performing variable selection ...",Sys.time())
source("./R/train.R")
y <- as.factor(anno$`methylation class:ch1`)

# sd pre filtering to 20k probes, to speed up the example
betas <- betas[,order(-apply(betas,2,sd))[1:20000]]

set.seed(seed,kind ="L'Ecuyer-CMRG") 
message("seed: ",seed)
message("cores: ",cores)
message("ntrees: ",ntrees)  
message("n: ",nrow(betas))
message("p: ",ncol(betas))  

rf.varsel <- rfp(betas,
                 y,
                 mc=cores,
                 ntree=ntrees,
                 sampsize=rep(min(table(y)),length(table(y))),
                 importance=TRUE)

# get permutation variable importance
imp.meandecrease <- rf.varsel$importance[,dim(rf.varsel$importance)[2]-1]

# save selection forest
save(rf.varsel,file="./results/varsel.RData")
rm(rf.varsel)

# reduce data matrix
or <- order(imp.meandecrease,decreasing=T)
betasy <- betas[,or[1:p]]

gc()

message("finished ...",Sys.time())

message("training classifier ...",Sys.time())

message("single core")
message("ntrees: ",ntrees)  
message("n: ",nrow(betasy))
message("p: ",ncol(betasy))

rf.pred <- randomForest(betasy,
                        y,
                        #mc=cores,
                        ntree=ntrees,
                        #strata=y,
                        #mtry=sqrt(ncol(betas)),
                        sampsize=rep(min(table(y)),length(table(y))),
                        proximity=TRUE,
                        oob.prox=TRUE,
                        importance=TRUE,
                        keep.inbag=TRUE,
                        do.trace=FALSE,
                        seed=seed
)

message("finished ...",Sys.time())

save(rf.pred,file="./results/rf.pred.RData")

message("finished ...",Sys.time())
