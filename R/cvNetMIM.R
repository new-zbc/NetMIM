#### predict response of test data
predictY <- function(result, n.burnin, max_iters, data){
  Y = data$Y
  C = as.matrix(data$C)
  E = as.matrix(data$E)
  M = as.matrix(data$M)
  map = data$map
  n = dim(C)[1]
  p = dim(E)[2]
  output = matrix(NA, nrow = n, ncol=(max_iters - n.burnin))
  for(j in 1:(max_iters - n.burnin)){
    E_hat = matrix(0, n, p)
    
    for(k in 1:p){
      E_hat[ ,k] = as.matrix(M[, which(map == k)]) %*% result$Omega[j, which(map == k)]
    }
    
    output[,j] = C %*% result$beta_c[j, ] + cbind((E - E_hat), E_hat) %*% result$beta[j, ]
  }
  
  est = rowMeans(output)
  return(est)
}

cvNIBAMM <- function(data, full_data_index, prior, n.burnin, max_iters, n.fold = 5){
  source("R/NIBAMM.R")
  library(caret)
  folds = createFolds(full_data_index, k=n.fold)
  mse1 = rep(0, n.fold)
  mse2 = rep(0, n.fold)
  for( i in 1:n.fold){
    folds[[i]] = full_data_index[folds[[i]]]
    fold_test = list()
    fold_test$Y = data$Y[folds[[i]]]
    fold_test$C = data$C[folds[[i]], ]
    fold_test$E = data$E[folds[[i]], ]
    fold_test$M = data$M[folds[[i]], ]
    fold_test$map = map
    
    fold_train = list()
    fold_train$Y = data$Y[-folds[[i]]]
    fold_train$C = data$C[-folds[[i]], ]
    fold_train$E = data$E[-folds[[i]], ]
    fold_train$M = data$M[-folds[[i]], ]
    fold_train$map = map
    
    
    result1 = NIBAMM(fold_train, prior = prior, max_iters = max_iters)
   
    hatY = predictY(result1, n.burnin, max_iters, data = fold_test)
    mse1[i] = mean((hatY - fold_test$Y)^2)
    
    
    #### complete case
    fold_train_cc = list()
    sampleid = setdiff(1:length(data$Y), full_data_index)
    sampleid = union(sampleid, folds[[i]])
    fold_train$Y = data$Y[-sampleid]
    fold_train$C = data$C[-sampleid, ]
    fold_train$E = data$E[-sampleid, ]
    fold_train$M = data$M[-sampleid, ]
    fold_train$map = map
    
    result2 = NIBAMM(fold_train, prior = prior, max_iters = max_iters)
    
    hatY = predictY(result2, n.burnin, max_iters, data = fold_test)
    mse2[i] = mean((hatY - fold_test$Y)^2)
  }
  
  if(mean(mse1) < mean(mse2)){
    return(NIBAMM(data, prior = prior, max_iters = max_iters))
  }
  else{
    dataNew = list()
    sampleid = setdiff(1:length(data$Y), full_data_index)
    dataNew$Y = data$Y[-sampleid]
    dataNew$C = data$C[-sampleid, ]
    dataNew$E = data$E[-sampleid, ]
    dataNew$M = data$M[-sampleid, ]
    dataNew$map = map
    
    return(NIBAMM(dataNew, prior = prior, max_iters = max_iters))
  }
}