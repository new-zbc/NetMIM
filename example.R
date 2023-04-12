
clinical_data = read.table(file = "example_data/clinical_covariates.txt")

##first variable is the response
Y = clinical_data[, 1]
C = clinical_data[, -1]

E = read.table(file = "example_data/gene_expression.txt", row.names = 1)
M = read.table(file = "example_data/DNA_methylation.txt", row.names = 1)
map = read.table(file = "example_data/map.txt", header = TRUE)
gene_network = read.table(file = "example_data/gene_network.txt", header = FALSE)
gene_network = as.matrix(gene_network)

## construct data set
data = list()
data$Y = Y
data$C = C
data$E = E
data$M = M
data$map = map

p = dim(E)[2]
n = dim(E)[1]
### construct prior list
prior = list(tau_c = 0.001,
             tau = 1,
             tau_g = rep(1,p),
             delta1 = 0.001,
             delta2 = 0.001,
             a = -3,
             b = 1,
             e = 0.5,
             f = 0.5,
             graph = gene_network)

source("R/NetMIM.R")

result = NIBAMM(data, prior = prior, max_iters = 10000)
colMeans(result$gamma)


### cross validation
source("R/cvNetMIM.R")
full_data_index = intersect(which(!is.na(E[, 1])), which(!is.na(M[, 1])))

result = cvNIBAMM(data, full_data_index, prior, n.burnin = 2000, max_iters = 8000)
colMeans(result$gamma)


#### for survival response
### we simulate a censor indicator
censor_indicator = c(rep(1, 70), rep(0, n-70))
Y = cbind(Y, censor_indicator)
data = list()
data$Y = Y
data$C = C
data$E = E
data$M = M
data$map = map

source("R/survivalNetMIM.R")
result = survivalNIBAMM(data, prior = prior, max_iters = 10000)
colMeans(result$gamma)
