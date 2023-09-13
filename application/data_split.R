##### split the data into training data and test data
ratio = 0.3
load("real_case_2/data/real_application_data_800.RData")
p = dim(inputeGraph)[1]
set.seed(0)
survival = survival[-which(row.names(survival) %in% c("100", "147", "446", "50", "53", "265", "5", "48", "49", "26", "202", "205", "206", "221", "483")), ]
mean(is.na(survival$days_to_death))
mean(survival$missing_index == 1)
gender = rep(1, dim(survival)[1])
gender[which(survival$gender == "male")] = 0
survival$gender = gender
survival$age_at_index = (survival$age_at_index - mean(survival$age_at_index))/sd(survival$age_at_index)

#### instructed missing data imputation in methylation 

for(j in 1:dim(MethylData)[2]){
  MethylData[is.na(MethylData[,j]),j] = mean(MethylData[, j], na.rm = T)
}

MethylData = log2(MethylData) / (1 - log2(MethylData))

#### gene expression normalization
for(j in 1:dim(ExpData)[2]){
  ExpData[,j] = (ExpData[,j]- mean(ExpData[, j])) /sd(ExpData[, j])
}


###### obtain test data set
survival_cc = survival[survival$missing_index == 0, ]
mean(is.na(survival_cc$days_to_death))

testSample_death = sample(which(!is.na(survival_cc$days_to_death)), size = length(which(!is.na(survival_cc$days_to_death))) * ratio)
testSample_censored = sample(which(is.na(survival_cc$days_to_death)), size = length(which(is.na(survival_cc$days_to_death))) * ratio)

testSample = sort(union(testSample_death, testSample_censored))

survivalTest = survival_cc[testSample, ]
ExpDataTest = ExpData[row.names(ExpData) %in% row.names(survivalTest), ]
MethylDataTest = MethylData[row.names(MethylData) %in% row.names(survivalTest), ]
### check the sample
row.names(survivalTest) == row.names(ExpDataTest)
row.names(survivalTest) == row.names(MethylDataTest)

save(survivalTest, ExpDataTest, MethylDataTest, probesMap, file = "real_case_2/data/TestDataSet.RData")


#### obtain training data without missing data
trainSample = setdiff(1:dim(survival_cc)[1], testSample)
survivalTrain = survival_cc[trainSample, ]
ExpDataTrain = ExpData[row.names(ExpData) %in% row.names(survivalTrain), ]
MethylDataTrain = MethylData[row.names(MethylData) %in% row.names(survivalTrain), ]
### check the sample
row.names(survivalTrain) == row.names(ExpDataTrain)
row.names(survivalTrain) == row.names(MethylDataTrain)

save(survivalTrain, ExpDataTrain, MethylDataTrain, probesMap, inputeGraph, file = "real_case_2/data/TrainDataSet.RData")


#### obtain training data with missing data
testSampleNames = row.names(survival_cc[testSample, ])
trainSampleNames = setdiff(row.names(survival), testSampleNames)

survivalTrainMissing = survival[row.names(survival) %in% trainSampleNames, ]
ExpDataTrainMissing = ExpData[row.names(ExpData) %in% trainSampleNames, ]
MethylDataTrainMissing  = MethylData[row.names(MethylData) %in% trainSampleNames, ]

##### NA replace the missing case

# dataExp = matrix(NA, nrow = dim(survivalTrainMissing)[1], ncol = dim(ExpDataTrainMissing)[2])
# dataExp[which(survivalTrainMissing$missing_index!=2), ] = ExpDataTrainMissing
dataMet = matrix(NA, nrow = dim(survivalTrainMissing)[1], ncol = dim(MethylDataTrainMissing)[2])
row.names(dataMet) = row.names(survivalTrainMissing)
dataMet[which(survivalTrainMissing$missing_index!=1), ] = MethylDataTrainMissing
colnames(dataMet) = colnames(MethylDataTrainMissing)
### check sample
row.names(survivalTrainMissing[which(survivalTrainMissing$missing_index!=1),]) == row.names(dataMet[which(survivalTrainMissing$missing_index!=1),])
sum(is.na(dataMet[which(survivalTrainMissing$missing_index!=1),]))
row.names(dataMet[which(!is.na(dataMet[,1])), ]) == row.names(MethylDataTrainMissing)


dataExp = matrix(NA, nrow = dim(survivalTrainMissing)[1], ncol = dim(ExpDataTrainMissing)[2])
row.names(dataExp) = row.names(survivalTrainMissing)
dataExp[which(survivalTrainMissing$missing_index !=2), ] = ExpDataTrainMissing
colnames(dataExp) = colnames(ExpDataTrainMissing)

### check sample
row.names(survivalTrainMissing[which(survivalTrainMissing$missing_index!=2),]) == row.names(dataExp[which(survivalTrainMissing$missing_index!=2),])
sum(is.na(dataExp[which(survivalTrainMissing$missing_index!=2),]))
row.names(dataExp[which(!is.na(dataExp[,1])), ]) == row.names(ExpDataTrainMissing)


ExpDataTrainMissing = dataExp
MethylDataTrainMissing = dataMet

save(survivalTrainMissing, ExpDataTrainMissing, MethylDataTrainMissing, probesMap,inputeGraph, file = "real_case_2/data/TrainDataSetMissing.RData")


