library(doParallel)
library(foreach)
n.fold = 5
load("real_case_2/data/TrainDataSet.RData")
C = cbind(rep(1, dim(survivalTrain)[1]),survivalTrain$age_at_index, survivalTrain$gender)

death = as.numeric(!is.na(survivalTrain$days_to_death))
Y = survivalTrain$days_to_last_follow_up

Y[!is.na(survivalTrain$days_to_death)] = survivalTrain$days_to_death[!is.na(survivalTrain$days_to_death)]

Y = log(Y)

data = list()
data$Y = cbind(Y, death)
data$C = C
data$E = ExpDataTrain
data$M = MethylDataTrain
data$map = probesMap

index1 = which(death == 0)
index2 = which(death == 1)

library(caret)
set.seed(0)
folds1 = createFolds(index1, k=n.fold)
folds2 = createFolds(index2, k=n.fold)

folds = list()
for(k in 1:n.fold){
  result = union(index1[folds1[[k]]], index2[folds2[[k]]])
  folds = append(folds, list(result))
}

names(folds) = paste0("Fold", 1:n.fold)

data_train = list()
data_val = list()
for(i in 1:n.fold){
  data1 = list()
  data1$Y = data$Y[-folds[[i]], ]
  data1$C = data$C[-folds[[i]], ]
  data1$E = data$E[-folds[[i]], ]
  data1$M = data$M[-folds[[i]], ]
  data1$map = data$map
  data_train = append(data_train, list(data1))
  
  data2 = list()
  data2$Y = data$Y[folds[[i]], ]
  data2$C = data$C[folds[[i]], ]
  data2$E = data$E[folds[[i]], ]
  data2$M = data$M[folds[[i]], ]
  data2$map = data$map
  
  data_val = append(data_val, list(data2))
}


prior = list(tau_c = 0.001, 
             tau = 1, 
             tau_g = rep(1,dim(inputeGraph)[1]), 
             delta1 = 0.001, 
             delta2 = 0.001, 
             a = -4, 
             b = 0.5, 
             e = 0.2, 
             f = 0.8, 
             graph = inputeGraph)


n.cluster = n.fold

result_wd = "real_case_2/complete_case/cv_multipleGamma40_b_0.5_a_4"

if(!file.exists(result_wd)){dir.create(result_wd)}


cl = makeCluster(n.cluster)
registerDoParallel(cl)
t1 = Sys.time()
mydata1 <- foreach(i=1:n.cluster, .noexport = c('<Functions that were implemented in C++>')) %dopar%
{
  source("graphIBAG/iBAG_survival_v3.R")
  result1 = GraphIBAG_survival(data=data_train[[i]], prior = prior, max_iters = 10000, seed = 1)
  filename = paste0(result_wd, "/", "fold", i, ".RData")
  save(result1, file = filename)
  result1
}
stopCluster(cl)
t2 = Sys.time()
print(t2-t1)

filename = paste0(result_wd, "/", "cross_validation.RData")
save(mydata1, file = filename)


#### calculate the c index in validation data set
predictionSurv <- function(result, n.burnin, max_iters, data){
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


c_index = rep(0, n.fold)
library(survival)
library(Hmisc)

for(i in 1:n.fold){
  test_pred = predictionSurv(mydata1[[i]], 5000, 10000, data_val[[i]])
  c_index[i] = rcorr.cens(test_pred, Surv(data_val[[i]]$Y[,1], data_val[[i]]$Y[,2]))[1]
}

print(c_index)
filename = paste0(result_wd, "/", "c_index.RData")
save(c_index, file = filename)



