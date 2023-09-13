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

load("real_case_2/missing_case/multipleGamma30_b_0.5_a_4/four_chain.RData")

###### calculate the correlation between different chains
gamma1 = colMeans(mydata1[[1]]$gamma[50000:400000, ])
gamma2 = colMeans(mydata1[[2]]$gamma[50000:400000, ])
gamma3 = colMeans(mydata1[[3]]$gamma[50000:400000, ])
gamma4 = colMeans(mydata1[[4]]$gamma[50000:400000, ])

cor(gamma1, gamma2)
cor(gamma1, gamma3)
cor(gamma1, gamma4)
cor(gamma2, gamma3)
cor(gamma2, gamma4)
cor(gamma3, gamma4)


load("real_case_2/data/TrainDataSetMissing.RData")

C = cbind(rep(1, dim(survivalTrainMissing)[1]),survivalTrainMissing$age_at_index, survivalTrainMissing$gender)

death = as.numeric(!is.na(survivalTrainMissing$days_to_death))
Y = survivalTrainMissing$days_to_last_follow_up

Y[!is.na(survivalTrainMissing$days_to_death)] = survivalTrainMissing$days_to_death[!is.na(survivalTrainMissing$days_to_death)]

Y = log(Y)

dataTrain = list()
dataTrain$Y = cbind(Y, death)
dataTrain$C = C
dataTrain$E = ExpDataTrainMissing
dataTrain$M = MethylDataTrainMissing
dataTrain$map = probesMap

index = which(survivalTrainMissing$missing_index == 0)
dataTrain$Y = dataTrain$Y[index, ]
dataTrain$C = dataTrain$C[index, ]
dataTrain$E = dataTrain$E[index, ]
dataTrain$M = dataTrain$M[index, ]

library(doParallel)
library(foreach)
cl = makeCluster(2)
registerDoParallel(cl)

train_pred <- foreach(i=1:4, .combine = "cbind") %dopar%
{
  predictionSurv(mydata1[[i]], 5000, 10000, dataTrain)
}
stopCluster(cl)


library(survival)
library(Hmisc)
c_index_train = rcorr.cens(rowMeans(train_pred), Surv(dataTrain$Y[, 1], dataTrain$Y[, 2]))



load("real_case_2/data/TestDataSet.RData")

C = cbind(rep(1, dim(survivalTest)[1]), survivalTest$age_at_index, survivalTest$gender)

death = as.numeric(!is.na(survivalTest$days_to_death))
Y = survivalTest$days_to_last_follow_up

Y[!is.na(survivalTest$days_to_death)] = survivalTest$days_to_death[!is.na(survivalTest$days_to_death)]

Y = log(Y)

dataTest = list()
dataTest$Y = cbind(Y, death)
dataTest$C = C
dataTest$E = ExpDataTest
dataTest$M = MethylDataTest
dataTest$map = probesMap


library(doParallel)
library(foreach)
cl = makeCluster(2)
registerDoParallel(cl)

test_pred <- foreach(i=1:4, .combine = "cbind") %dopar%
{
  predictionSurv(mydata1[[i]], 5000, 10000, dataTest)
}
stopCluster(cl)


library(survival)
library(Hmisc)



c_index_test = rcorr.cens(rowMeans(test_pred), Surv(dataTest$Y[,1], dataTest$Y[,2]))

filename = "real_case_2/missing_case/multipleGamma30_b_0.5_a_4/predict.RData"
save(train_pred, c_index_train, test_pred, c_index_test, file = filename)