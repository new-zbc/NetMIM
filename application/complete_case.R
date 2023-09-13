
library(doParallel)
library(foreach)
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


seeds = c(1, 12, 123, 1234)
n.cluster = length(seeds)

result_wd = "real_case_2/complete_case/multipleGamma40_b_0.5_a_4"

if(!file.exists(result_wd)){dir.create(result_wd)}


cl = makeCluster(n.cluster)
registerDoParallel(cl)
t1 = Sys.time()
mydata1 <- foreach(i=1:n.cluster, .noexport = c('<Functions that were implemented in C++>')) %dopar%
{
  source("graphIBAG/iBAG_survival_v3.R")
  result1 = GraphIBAG_survival(data=data, prior = prior, max_iters = 10000, seed = seeds[i])
  filename = paste0(result_wd, "/", "chain", i, ".RData")
  save(result1, file = filename)
  result1
}
stopCluster(cl)
t2 = Sys.time()
print(t2-t1)

filename = paste0(result_wd, "/", "four_chain.RData")
save(mydata1, file = filename)
