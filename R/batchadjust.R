#-----------------------------------------------------------------------------------
# script to include the batchadjustment into the cross validation
#                                                                     
# Martin Sill
# m.sill@dkfz.de                                                                  
# 
# 2018-03-14 UTC
#------------------------------------------------------------------------------------                   
#options(max.print = 1000)
#options(stringsAsFactors = FALSE)
#options(scipen = 999)


batchadjust <- function(Mset,batchall,fold){
  methy <- getMeth(Mset[,fold$train])
  unmethy <- getUnmeth(Mset[,fold$train])
  batch <- batchall[fold$train]
  methy.ba <- 2^removeBatchEffect(log2(methy +1), batch)
  unmethy.ba <- 2^removeBatchEffect(log2(unmethy +1), batch)
  # recalculate betas, illumina like
  betas <- methy.ba / (methy.ba +unmethy.ba +100)
  betas <- as.data.frame(t(betas))
  # extract coefficients that can be applied to diagnostic samples
  s.frozen <- min(which(batch == "Frozen"))
  s.ffpe <- min(which(batch == "FFPE"))
  methy.coef <- unmethy.coef <- list()
  methy.coef[["Frozen"]] <- log2(methy.ba[, s.frozen]) - log2(methy[, s.frozen] +1)
  methy.coef[["FFPE"]] <- log2(methy.ba[, s.ffpe]) - log2(methy[, s.ffpe] +1)
  unmethy.coef[["Frozen"]] <- log2(unmethy.ba[, s.frozen]) - log2(unmethy[, s.frozen] +1)
  unmethy.coef[["FFPE"]] <- log2(unmethy.ba[, s.ffpe]) - log2(unmethy[, s.ffpe] +1)
  methy <- getMeth(Mset[,fold$test])
  unmethy <- getUnmeth(Mset[,fold$test])
  batch <- batchall[fold$test]
  # perform batch adjustment
  methy.b <- log2(methy +1) + matrix(unlist(methy.coef[match(batch,names(methy.coef))]),ncol=length(batch))
  unmethy.b <- log2(unmethy +1) + matrix(unlist(unmethy.coef[match(batch,names(unmethy.coef))]),ncol=length(batch))
  methy.b[methy.b < 0] <- 0
  unmethy.b[unmethy.b < 0] <- 0
  methy.ba <- 2^methy.b
  unmethy.ba <- 2^unmethy.b
  # illumina-like beta values
  betas.test <- methy.ba / (methy.ba +unmethy.ba +100)
  betas.test <- as.data.frame(t(betas.test))
  return(list(betas.train=betas,betas.test=betas.test))
}
