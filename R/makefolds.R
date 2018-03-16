makefolds <- function(y,cv.fold=3){
  n <- length(y)
  nlvl <- table(y)
  idx <- numeric(n)
  folds <- list()
  for (i in 1:length(nlvl)) {
    idx[which(y == levels(y)[i])] <- sample(rep(1:cv.fold,length = nlvl[i]))
  }
  for (i in 1:cv.fold){
    folds[[i]] <- list(train = which(idx!=i),
                       test =  which(idx==i)) 
  }  
  return(folds)
}
makenestedfolds <- function(y,cv.fold=3){
  nfolds <- list()
  folds <- makefolds(y,cv.fold)
  names(folds) <- paste0("outer",1:length(folds))
  for(k in 1:length(folds)){
    inner = makefolds(y[folds[[k]]$train],cv.fold)
    names(inner) <- paste0("inner",1:length(folds))
    for(i in 1:length(inner)){
      inner[[i]]$train <- folds[[k]]$train[inner[[i]]$train]
      inner[[i]]$test <- folds[[k]]$train[inner[[i]]$test]
    }
    nfolds[[k]] <- list(folds[k],inner) 
  }
  names(nfolds) <- paste0("outer",1:length(nfolds))
  return(nfolds)
}

# test if any outer test data in inner folds
#test <- logical()
#for(i in 1:length(nfolds)) test[i] <- any(nfolds[[i]][[1]][[1]]$test %in% unlist(nfolds[[i]][[2]]))
#any(test)
