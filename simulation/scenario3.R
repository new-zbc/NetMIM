load("simulation/data.Rdata")
library(doParallel)
library(foreach)
library(MASS)
n.cluster = 15

#### Assessment function for variable selection #######
selection.perf <- function(gamma, gamma_truth, threshold = 0.5)
{
  library(ROCR)
  p = length(gamma)
  subset.selection = rep(0, p)
  subset.selection[gamma >= threshold] = 1
  if(sum(subset.selection) >= p)
  {
    if(sum(gamma) == p){
      return(c(sensitivity = 1, 
               Specificity = 0, 
               MCC = NA, 
               AUC = NA))
    }
    else{
      pred.auc = prediction(gamma, gamma_truth)
      perf <- performance(pred.auc,'auc')
      auc = perf@y.values[[1]]
      return(c(sensitivity = 1, 
               Specificity = 0, 
               MCC = NA, 
               AUC = auc))
    }
    
  }
  
  if(sum(subset.selection) == 0)
  {
    if(sum(gamma) == 0){
      return(c(sensitivity = 0, 
               Specificity = 1, 
               MCC = NA, 
               AUC = NA))
    }
    else{
      pred.auc = prediction(gamma, gamma_truth)
      perf <- performance(pred.auc,'auc')
      auc = perf@y.values[[1]]
      return(c(sensitivity = 0, 
               Specificity = 1, 
               MCC = NA, 
               AUC = auc))
    }
  }
  
  result = table(subset.selection, gamma_truth)
  sens = result[2,2]/sum(result[,2])
  spec = result[1,1]/sum(result[,1])
  mcc = (result[2,2]*result[1,1] - result[1,2]*result[2,1]) / sqrt(sum(result[1,]) * sum(result[2,]))
  mcc = mcc /  sqrt(sum(result[,1]) * sum(result[,2]))
  pred.auc = prediction(gamma, gamma_truth)
  perf <- performance(pred.auc,'auc')
  auc = perf@y.values[[1]]
  return(c(sensitivity = sens, 
           Specificity = spec, 
           MCC = mcc, 
           AUC = auc))
}

#### predict response of test data
predictY <- function(result, n.burnin, max_iters, data){
  Y = data$Y
  C = data$C
  E = data$E
  M = data$M
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


#### generate data
simData <- function(n, gamma, Z, map, Methyl, sd=1, seed = 0){
  
  library(MASS)
  p = as.integer(length(gamma)/2)
  pm = length(Z)
  set.seed(seed)
  sampleID = sample(n)
  
  M = Methyl[sampleID, sample(pm)]
  
  Omega = rep(0, pm)
  for(j in 1:pm){
    if(Z[j] == 1){Omega[j] = sample(c(-1,1), size = 1) * runif(1, 0.5, 1)}
  }
  
  E_hat = matrix(0, nrow = n, ncol = p)
  for(j in 1:p){
    idx = which(map == j)
    E_hat[, j] = M[ ,idx] %*% as.matrix(Omega[idx])
  }
  
  E_bar = mvrnorm(n, mu=rep(0, p), Sigma = 1*diag(1, p) + 0*matrix(1, p, p))
  E = E_hat + E_bar
  
  
  C = mvrnorm(n=n, mu=rep(0, 3), Sigma = diag(1, 3))
  C = cbind(rep(1, n), C)
  beta_c = rep(1,4)
  
  beta = rep(0, 2*p)
  beta[which(gamma == 1)] = sample(c(-1,1), size = sum(gamma), replace = TRUE) * runif(sum(gamma), min = 1, max = 1.5)
  
  Y = as.vector(C %*% beta_c + cbind(E_bar, E_hat) %*% beta + rnorm(n, mean = 0, sd=sd))
  
  testSample  = sample(n, size = 50)
  trainSample = setdiff(1:n, testSample)
  
  data = list()
  data$Y = Y[trainSample]
  data$C = C[trainSample, ]
  data$E = E[trainSample, ]
  data$M = M[trainSample, ]
  data$map = map
  
  dataTest = list()
  dataTest$Y = Y[testSample]
  dataTest$C = C[testSample, ]
  dataTest$E = E[testSample, ]
  dataTest$M = M[testSample, ]
  dataTest$map = map
  
  return(list(train=data, test=dataTest))
}



result_wd = "simulation/scenario3"

if(!file.exists(result_wd)){dir.create(result_wd)}

res.sum = data.frame()

for(p in c(100, 200, 500)){
  for(ratio in c(0.1, 0.2, 0.3, 0.4, 0.5)){
    
    n = 100
    p = p
    q = 3
    pm = p*q
    set.seed(0)
    
    
    # generate gamma
    gamma_ground_truth = c(1:10, p+(1:5), p+(11:15))
    gamma = rep(0, 2*p)
    gamma[gamma_ground_truth] = 1
    ## generate map
    map = rep(0, pm)
    for(j in 1:p){
      map[q*j] = j
      #map[q*j -1] = j
      for(k in 1:(q-1)){
        map[q*j - k] = sample(1:p,size = 1)
      }
    }
    map = sort(map)
    
    Z = rep(0, pm)
    ## generate Z
    prob = 0.6
    for(j in 1:p){
      if(j < 20){
        Z[q*j] = 1
        for(k in 1:(q-1)){
          Z[q*j - k] = sample(c(0,1),size = 1, prob = c(prob, 1-prob))
        }
      }
      else{
        Z[q*j] = sample(c(0,1),size = 1, prob = c(prob, 1-prob))
        for(k in 1:(q-1)){
          Z[q*j - k] = sample(c(0,1),size = 1, prob = c(prob, 1-prob))
        }
      }
    }
    
    ### graph
    graph = diag(1, p, p)
    for(i in 1:15){
      for(j in 1:(i-1)){
        graph[i,j] = graph[j, i] = sample(c(0,1), size = 1, prob = c(0.9, 0.1))
      }
    }
    
    # ###### our method #######
    # if(p < 500){
    #   prior = list(tau_c = 0.001,
    #                tau = 1,
    #                tau_g = rep(1,p),
    #                delta1 = 0.001,
    #                delta2 = 0.001,
    #                a = -3,
    #                b = 0.5,
    #                e = 0.5,
    #                f = 0.5,
    #                graph = graph)
    # }else{
    #   prior = list(tau_c = 0.001,
    #                tau = 1,
    #                tau_g = rep(1,p),
    #                delta1 = 0.001,
    #                delta2 = 0.001,
    #                a = -4,
    #                b = 0.5,
    #                e = 0.5,
    #                f = 0.5,
    #                graph = graph)
    # }
    # 
    # 
    # cl = makeCluster(n.cluster)
    # registerDoParallel(cl)
    # mydata1 <- foreach(i=1:(2*n.cluster), .combine=rbind, .noexport = c('<Functions that were implemented in C++>')) %dopar%
    #   {
    #     source("graphIBAG/iBAG_v2.R")
    #     
    #     n.burnin = 2000
    #     max_iters = 10000
    #     data_sim = simData(n+50, gamma, Z, map, MethylSample1, sd=1, seed = i^2)
    #     
    #     sampleid = sample(n, size = n*ratio)
    #     
    #     data_sim$train$E[sampleid, ] <- NA
    #     
    #     result = GraphIBAG(data_sim$train, prior = prior, max_iters = max_iters)
    #     varSelection = colMeans(result$gamma[(5*n.burnin):dim(result$gamma)[1], ])
    #     perf = selection.perf(varSelection, gamma)
    #     hatY = predictY(result, n.burnin, max_iters, data = data_sim$test)
    #     mse = mean((hatY - data_sim$test$Y)^2)
    #     c(perf, mse, varSelection)
    #     
    #   }
    # stopCluster(cl)
    # 
    # filename = paste0(result_wd, "/", "graphIBAG", "_ratio_", ratio, "_dim_", p, ".txt")
    # write.table(mydata1, file = filename, quote = FALSE, row.names = FALSE)
    # 
    # vec1 = colMeans(mydata1, na.rm = T)
    # vec2 = apply(mydata1, 2, sd, na.rm = T)
    # 
    # 
    # print(vec1[1:5])
    # print(vec2[1:5])
    # 
    # res.sum = rbind(res.sum, vec1)
    # res.sum = rbind(res.sum, vec2)
    # 
    # 
    # ###### our method with b = 0 #######
    # if(p < 500){
    #   prior = list(tau_c = 0.001,
    #                tau = 1,
    #                tau_g = rep(1,p),
    #                delta1 = 0.001,
    #                delta2 = 0.001,
    #                a = -3,
    #                b = 0,
    #                e = 0.5,
    #                f = 0.5,
    #                graph = graph)
    # }else{
    #   prior = list(tau_c = 0.001,
    #                tau = 1,
    #                tau_g = rep(1,p),
    #                delta1 = 0.001,
    #                delta2 = 0.001,
    #                a = -4,
    #                b = 0,
    #                e = 0.5,
    #                f = 0.5,
    #                graph = graph)
    # }
    # 
    # 
    # cl = makeCluster(n.cluster)
    # registerDoParallel(cl)
    # mydata2 <- foreach(i=1:(2*n.cluster), .combine=rbind, .noexport = c('<Functions that were implemented in C++>')) %dopar%
    #   {
    #     source("graphIBAG/iBAG_v2.R")
    #     
    #     n.burnin = 2000
    #     max_iters = 10000
    #     data_sim = simData(n+50, gamma, Z, map, MethylSample1, sd=1, seed = i^2)
    #     
    #     sampleid = sample(n, size = n*ratio)
    #     
    #     data_sim$train$E[sampleid, ] <- NA
    #     
    #     result = GraphIBAG(data_sim$train, prior = prior, max_iters = max_iters)
    #     varSelection = colMeans(result$gamma[(5*n.burnin):dim(result$gamma)[1], ])
    #     perf = selection.perf(varSelection, gamma)
    #     hatY = predictY(result, n.burnin, max_iters, data = data_sim$test)
    #     mse = mean((hatY - data_sim$test$Y)^2)
    #     c(perf, mse, varSelection)
    #     
    #   }
    # stopCluster(cl)
    # 
    # filename = paste0(result_wd, "/", "graphIBAG_0", "_ratio_", ratio, "_dim_", p,".txt")
    # write.table(mydata2, file = filename, quote = FALSE, row.names = FALSE)
    # 
    # vec1 = colMeans(mydata2, na.rm = T)
    # vec2 = apply(mydata2, 2, sd, na.rm = T)
    # 
    # 
    # print(vec1[1:5])
    # print(vec2[1:5])
    # 
    # res.sum = rbind(res.sum, vec1)
    # res.sum = rbind(res.sum, vec2)
    
    
    #### FBM method #######
    cl = makeCluster(n.cluster)
    registerDoParallel(cl)
    mydata3 <- foreach(i=1:(2*n.cluster), .combine=rbind, .noexport = c('<Functions that were implemented in C++>')) %dopar%
      {
        library(Rcpp)
        library(RcppArmadillo)
        library(RcppGSL)
        library(scales)
        library(MASS)
        Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
        sourceCpp("FBM/FBM.cpp")
        n.burnin = 2000
        max_iters = 10000
        
        dataList = simData(n+50, gamma, Z, map, MethylSample1, sd=1, seed = i^2)
        
        data_sim = dataList$train
        
        sampleid = sample(n, size = n*ratio)
        
        data <- list()
        data$Y <-  data_sim$Y
        
        data$geneObs <- data_sim$E[-sampleid, ]
        row.names(data$geneObs) <- setdiff(1:length(data$Y), sampleid)
        
        
        data$methyObs <- data_sim$M
        data$methyObs <- as.matrix(data$methyObs)
        row.names(data$methyObs) <- 1:dim(data$methyObs)[1]
        
        data$mapMethy <-  map
        
        data$C <- data_sim$C
        row.names(data$C) = 1:dim(data_sim$C)[1]
        
        data$missGeneIdx <- sampleid
        data$obsGeneIdx <- setdiff(1:length(data$Y), data$missGeneIdx)
        
        data$missMethyIdx <- integer(0)
        data$obsMethyIdx <- setdiff(1:length(data$Y), data$missMethyIdx)
        
        result = FBMcpp(data, max_iters, 15217, "Full")
        varSelection = colMeans(cbind(result$IMbar[n.burnin:max_iters, ], result$IM[n.burnin:max_iters, ]))
        perf = selection.perf(varSelection, gamma)
        
        parEst = list(beta = cbind(result$gammaMbar, result$gammaM), Omega = result$omega, beta_c = result$gammaC)
        
        hatY = predictY(parEst, n.burnin, max_iters, data = dataList$test)
        mse = mean((hatY - dataList$test$Y)^2)
        c(perf, mse, varSelection)
        
        
      }
    stopCluster(cl)
    
    filename = paste0(result_wd, "/", "FBM", "_ratio_", ratio, "_dim_", p, ".txt")
    write.table(mydata3, file = filename, quote = FALSE, row.names = FALSE)
    
    vec1 = colMeans(mydata3, na.rm = T)
    vec2 = apply(mydata3, 2, sd, na.rm = T)
    
    
    print(vec1[1:5])
    print(vec2[1:5])
    
    res.sum = rbind(res.sum, vec1)
    res.sum = rbind(res.sum, vec2)
    
    
    ############## complete case #####################
    ###### our method remove missing data #######
    # if(p < 500){
    #   prior = list(tau_c = 0.001,
    #                tau = 1,
    #                tau_g = rep(1,p),
    #                delta1 = 0.001,
    #                delta2 = 0.001,
    #                a = -3,
    #                b = 0.5,
    #                e = 0.5,
    #                f = 0.5,
    #                graph = graph)
    # }else{
    #   prior = list(tau_c = 0.001,
    #                tau = 1,
    #                tau_g = rep(1,p),
    #                delta1 = 0.001,
    #                delta2 = 0.001,
    #                a = -4,
    #                b = 1,
    #                e = 0.5,
    #                f = 0.5,
    #                graph = graph)
    # }
    # 
    # 
    # cl = makeCluster(n.cluster)
    # registerDoParallel(cl)
    # mydata1 <- foreach(i=1:(2*n.cluster), .combine=rbind, .noexport = c('<Functions that were implemented in C++>')) %dopar%
    #   {
    #     source("graphIBAG/iBAG_v2.R")
    #     
    #     n.burnin = 2000
    #     max_iters = 10000
    #     data_sim = simData(n+50, gamma, Z, map, MethylSample1, sd=1, seed = i^2)
    #     
    #     sampleid = sample(n, size = n*ratio)
    #     
    #     data_sim$train$Y <- data_sim$train$Y[-sampleid]
    #     data_sim$train$C <- data_sim$train$C[-sampleid, ]
    #     data_sim$train$E <- data_sim$train$E[-sampleid, ]
    #     data_sim$train$M <- data_sim$train$M[-sampleid, ]
    #     
    #     result = GraphIBAG(data_sim$train, prior = prior, max_iters = max_iters)
    #     
    #     varSelection = colMeans(result$gamma[(5*n.burnin):dim(result$gamma)[1], ])
    #     perf = selection.perf(varSelection, gamma)
    #     hatY = predictY(result, n.burnin, max_iters, data = data_sim$test)
    #     mse = mean((hatY - data_sim$test$Y)^2)
    #     c(perf, mse, varSelection)
    #     
    #   }
    # stopCluster(cl)
    # 
    # filename = paste0(result_wd, "/", "graphIBAG_cc", "_ratio_", ratio, "_dim_", p, ".txt")
    # write.table(mydata1, file = filename, quote = FALSE, row.names = FALSE)
    # 
    # vec1 = colMeans(mydata1, na.rm = T)
    # vec2 = apply(mydata1, 2, sd, na.rm = T)
    # 
    # 
    # print(vec1[1:5])
    # print(vec2[1:5])
    # 
    # res.sum = rbind(res.sum, vec1)
    # res.sum = rbind(res.sum, vec2)
    # 
    # 
    # 
    # ###### our method with b = 0 remove missing data #######
    # if(p < 500){
    #   prior = list(tau_c = 0.001,
    #                tau = 1,
    #                tau_g = rep(1,p),
    #                delta1 = 0.001,
    #                delta2 = 0.001,
    #                a = -3,
    #                b = 0,
    #                e = 0.5,
    #                f = 0.5,
    #                graph = graph)
    # }else{
    #   prior = list(tau_c = 0.001,
    #                tau = 1,
    #                tau_g = rep(1,p),
    #                delta1 = 0.001,
    #                delta2 = 0.001,
    #                a = -4,
    #                b = 0,
    #                e = 0.5,
    #                f = 0.5,
    #                graph = graph)
    # }
    # 
    # 
    # cl = makeCluster(n.cluster)
    # registerDoParallel(cl)
    # mydata2 <- foreach(i=1:(2*n.cluster), .combine=rbind, .noexport = c('<Functions that were implemented in C++>')) %dopar%
    #   {
    #     source("graphIBAG/iBAG_v2.R")
    #     
    #     n.burnin = 2000
    #     max_iters = 10000
    #     data_sim = simData(n+50, gamma, Z, map, MethylSample1, sd=1, seed = i^2)
    #     
    #     sampleid = sample(n, size = n*ratio)
    #     
    #     data_sim$train$Y <- data_sim$train$Y[-sampleid]
    #     data_sim$train$C <- data_sim$train$C[-sampleid, ]
    #     data_sim$train$E <- data_sim$train$E[-sampleid, ]
    #     data_sim$train$M <- data_sim$train$M[-sampleid, ]
    #     
    #     result = GraphIBAG(data_sim$train, prior = prior, max_iters = max_iters)
    #     
    #     varSelection = colMeans(result$gamma[(5*n.burnin):dim(result$gamma)[1], ])
    #     perf = selection.perf(varSelection, gamma)
    #     hatY = predictY(result, n.burnin, max_iters, data = data_sim$test)
    #     mse = mean((hatY - data_sim$test$Y)^2)
    #     c(perf, mse, varSelection)
    #     
    #   }
    # stopCluster(cl)
    # 
    # filename = paste0(result_wd, "/", "graphIBAG_0_cc", "_ratio_", ratio, "_dim_", p,".txt")
    # write.table(mydata2, file = filename, quote = FALSE, row.names = FALSE)
    # 
    # vec1 = colMeans(mydata2, na.rm = T)
    # vec2 = apply(mydata2, 2, sd, na.rm = T)
    # 
    # 
    # print(vec1[1:5])
    # print(vec2[1:5])
    # 
    # res.sum = rbind(res.sum, vec1)
    # res.sum = rbind(res.sum, vec2)
    
    
    #### FBM method #######
    cl = makeCluster(n.cluster)
    registerDoParallel(cl)
    mydata3 <- foreach(i=1:(2*n.cluster), .combine=rbind, .noexport = c('<Functions that were implemented in C++>')) %dopar%
      {
        library(Rcpp)
        library(RcppArmadillo)
        library(RcppGSL)
        library(scales)
        library(MASS)
        Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
        sourceCpp("FBM/FBM.cpp")
        n.burnin = 2000
        max_iters = 10000
        
        dataList = simData(n+50, gamma, Z, map, MethylSample1, sd=1, seed = i^2)
        
        data_sim = dataList$train
        
        sampleid = sample(n, size = n*ratio)
        
        data <- list()
        data$Y <-  data_sim$Y[-sampleid]
        
        data$geneObs <- data_sim$E[-sampleid, ]
        row.names(data$geneObs) <- 1:dim(data$geneObs)[1]
        
        
        data$methyObs <- data_sim$M[-sampleid, ]
        data$methyObs <- as.matrix(data$methyObs)
        row.names(data$methyObs) <- 1:dim(data$methyObs)[1]
        
        data$mapMethy <-  map
        
        data$C <- data_sim$C[-sampleid, ]
        row.names(data$C) = 1:dim(data$C)[1]
        
        data$missGeneIdx <- integer(0)
        data$obsGeneIdx <- setdiff(1:length(data$Y), data$missGeneIdx)
        
        data$missMethyIdx <- integer(0)
        data$obsMethyIdx <- setdiff(1:length(data$Y), data$missMethyIdx)
        
        result = FBMcpp(data, max_iters, 15217, "Full")
        
        varSelection = colMeans(cbind(result$IMbar[n.burnin:max_iters, ], result$IM[n.burnin:max_iters, ]))
        perf = selection.perf(varSelection, gamma)
        
        parEst = list(beta = cbind(result$gammaMbar, result$gammaM), Omega = result$omega, beta_c = result$gammaC)
        
        hatY = predictY(parEst, n.burnin, max_iters, data = dataList$test)
        mse = mean((hatY - dataList$test$Y)^2)
        c(perf, mse, varSelection)
        
        
      }
    stopCluster(cl)
    
    filename = paste0(result_wd, "/", "FBM_cc", "_ratio_", ratio, "_dim_", p, ".txt")
    write.table(mydata3, file = filename, quote = FALSE, row.names = FALSE)
    
    vec1 = colMeans(mydata3, na.rm = T)
    vec2 = apply(mydata3, 2, sd, na.rm = T)
    
    
    print(vec1[1:5])
    print(vec2[1:5])
    
    res.sum = rbind(res.sum, vec1)
    res.sum = rbind(res.sum, vec2)
    
  }
}

colnames(res.sum)[1:5] = c("sens", "spec", "MCC", "AUC", "MSE")
filename = paste0(result_wd, "/", "summary_p_ratio.csv")
write.csv(res.sum, file = filename, quote = FALSE, row.names = FALSE)
