#step4:modelling################################################################

rm(list = ls())
library(caret)
library(randomForest)
library(ggplot2)
load(file = 'antimicrobial_toxin/rf_cluster/optvar51_cluster.Rdata')

set.seed(51)
dat <- read.csv(file = 'antimicrobial_toxin/train_new.csv',header = T)
dat2 <- na.omit(dat)
dat3 <- dat2[,-1:-4]
x <- dat3[,optvar]
logBBB <- dat2$toxin_clusters
load(file ="antimicrobial_toxin/rf_cluster/inTrain_cluster.Rdata")
trainx <- x[inTrain, ]
testx <- x[-inTrain, ]
trainy <- logBBB[inTrain]
testy <- logBBB[-inTrain]

customRF <- list(type = "Regression", library = "randomForest", loop = NULL,
                 parameters = data.frame(parameter = c("mtry", "ntree","nodesize"), 
                                         class = rep("numeric", 3), 
                                         label = c("mtry", "ntree","nodesize")),
                 
                 grid = function(x, y, len = NULL, search = "grid") {
                   if(search == "grid") {
                     out <- expand.grid(mtry = caret::var_seq(p = ncol(x),
                                                              classification = is.factor(y),
                                                              len = len),
                                        ntree = c(500,700,900,1000,1500),
                                        nodesize = c(4,5,6,7,8,9))
                   } else {
                     out <- data.frame(mtry = unique(sample(1:ncol(x), size = len, replace = TRUE)),
                                       ntree = unique(sample(c(500,700,900,1000,1500), 
                                                             size = len, replace = TRUE)),
                                       nodesize = unique(sample(c(4,5,6,7,8,9), 
                                                                size = len, replace = TRUE)))
                   }
                 },
                 fit = function(x, y, wts, param, lev, last, weights, classProbs, ...) {
                   randomForest(x, y, mtry = param$mtry, ntree=param$ntree, nodesize = param$nodesize, ...)
                 },
                 predict = function(modelFit, newdata, preProc = NULL, submodels = NULL)
                   predict(modelFit, newdata),
                 prob = function(modelFit, newdata, preProc = NULL, submodels = NULL)
                   predict(modelFit, newdata, type = "prob"),
                 sort = function(x) x[order(x[,1]),],
                 levels = function(x) x$classes
)

load(file = 'antimicrobial_toxin/rf_cluster/rf_gridsearch_cluster.Rdata')
testset <- testx
rsq <- function(x, y) summary(lm(y~x))$r.squared

modeldata <- data.frame(matrix(nrow = 100,ncol = 3))
control_grid <- trainControl(method = "cv", number = 10, search = "grid")

for (i in 1:100) {
  set.seed(i)
  rfmodel <- train(trainx,trainy, method=customRF,
                   tuneGrid=rf_gridsearch$bestTune, trControl=control_grid, 
                   metric = "Rsquared",impotance=T)
  modeldata[i,1] <- rfmodel[["results"]][["Rsquared"]]
  modeldata[i,2] <- rfmodel[["results"]][["RMSE"]]
  predata <- predict(rfmodel,newdata = testset)
  rsq <- function(x, y) summary(lm(y~x))$r.squared
  modeldata[i,3] <- rsq(testy, predata)
}

rfmodel1 <- rf_gridsearch$finalModel
rfmodel1
predictRF_1 <- predict(rfmodel1,newdata = testset)
rsq(testy, predictRF_1)

set.seed(69)
control_grid <- trainControl(method = "cv", number = 10, search = "grid")
rfmodel2 <- train(trainx,trainy, method=customRF,
                  tuneGrid=rf_gridsearch$bestTune, trControl=control_grid, 
                  metric = "Rsquared",impotance=T)
rfmodel2
predictRF_2 <- predict(rfmodel2,newdata = testset)
rsq(testy, predictRF_2)

set.seed(9)
control_grid <- trainControl(method = "cv", number = 10, search = "grid")
rfmodel3 <- train(trainx,trainy, method=customRF,
                  tuneGrid=rf_gridsearch$bestTune, trControl=control_grid, 
                  metric = "Rsquared",impotance=T)
rfmodel3
predictRF_3 <- predict(rfmodel3,newdata = testset)
rsq(testy, predictRF_3)

set.seed(4)
control_grid <- trainControl(method = "cv", number = 10, search = "grid")
rfmodel4 <- train(trainx,trainy, method=customRF,
                  tuneGrid=rf_gridsearch$bestTune, trControl=control_grid, 
                  metric = "Rsquared",impotance=T)
rfmodel4
predictRF_4 <- predict(rfmodel4,newdata = testset)
rsq(testy, predictRF_4)

set.seed(42)
control_grid <- trainControl(method = "cv", number = 10, search = "grid")
rfmodel5 <- train(trainx,trainy, method=customRF,
                  tuneGrid=rf_gridsearch$bestTune, trControl=control_grid, 
                  metric = "Rsquared",impotance=T)
rfmodel5
predictRF_5 <- predict(rfmodel5,newdata = testset)
rsq(testy, predictRF_5)

set.seed(80)
control_grid <- trainControl(method = "cv", number = 10, search = "grid")
rfmodel6 <- train(trainx,trainy, method=customRF,
                  tuneGrid=rf_gridsearch$bestTune, trControl=control_grid, 
                  metric = "Rsquared",impotance=T)
rfmodel6
predictRF_6 <- predict(rfmodel6,newdata = testset)
rsq(testy, predictRF_6)

set.seed(26)
control_grid <- trainControl(method = "cv", number = 10, search = "grid")
rfmodel7 <- train(trainx,trainy, method=customRF,
                  tuneGrid=rf_gridsearch$bestTune, trControl=control_grid, 
                  metric = "Rsquared",impotance=T)
rfmodel7
predictRF_7 <- predict(rfmodel7,newdata = testset)
rsq(testy, predictRF_7)

set.seed(7)
control_grid <- trainControl(method = "cv", number = 10, search = "grid")
rfmodel8 <- train(trainx,trainy, method=customRF,
                  tuneGrid=rf_gridsearch$bestTune, trControl=control_grid, 
                  metric = "Rsquared",impotance=T)
rfmodel8
predictRF_8 <- predict(rfmodel8,newdata = testset)
rsq(testy, predictRF_8)

set.seed(51)
control_grid <- trainControl(method = "cv", number = 10, search = "grid")
rfmodel9 <- train(trainx,trainy, method=customRF,
                  tuneGrid=rf_gridsearch$bestTune, trControl=control_grid, 
                  metric = "Rsquared",impotance=T)
rfmodel9
predictRF_9 <- predict(rfmodel9,newdata = testset)
rsq(testy, predictRF_9)

set.seed(94)
control_grid <- trainControl(method = "cv", number = 10, search = "grid")
rfmodel10 <- train(trainx,trainy, method=customRF,
                   tuneGrid=rf_gridsearch$bestTune, trControl=control_grid, 
                   metric = "Rsquared",impotance=T)
rfmodel10
predictRF_10 <- predict(rfmodel10,newdata = testset)
rsq(testy, predictRF_10)

mapset <- read.csv(file = 'antimicrobial_toxin/rf_cluster/map_cluster.csv')
predictRF1 <- predict(rfmodel1,newdata = mapset)
predictRF2 <- predict(rfmodel2,newdata = mapset)
predictRF3 <- predict(rfmodel3,newdata = mapset)
predictRF4 <- predict(rfmodel4,newdata = mapset)
predictRF5 <- predict(rfmodel5,newdata = mapset)
predictRF6 <- predict(rfmodel6,newdata = mapset)
predictRF7 <- predict(rfmodel7,newdata = mapset)
predictRF8 <- predict(rfmodel8,newdata = mapset)
predictRF9 <- predict(rfmodel9,newdata = mapset)
predictRF10 <- predict(rfmodel10,newdata = mapset)

predictsum <- data.frame(predictRF1,predictRF2,predictRF3,predictRF4,predictRF5,
                         predictRF6,predictRF7,predictRF8,predictRF9,predictRF10,
                         mapset$longitude,mapset$latitude)

for (i in 1:nrow(predictsum)){
  predictsum$prediction_sd[i] <- sd(as.numeric(predictsum[i,1:10]))
  predictsum$prediction_mean[i] <- mean(as.numeric(predictsum[i,1:10]))
}
predictsum$cv <- predictsum$prediction_sd / predictsum$prediction_mean
write.csv(predictsum,file = 'antimicrobial_toxin/rf_cluster/final_cluster.csv')

var6 <- varImp(rf_gridsearch, scale=FALSE)$importance  
var6$pRF <- var6$Overall / sum(var6$Overall) * 100
var6$rfround <- round(var6$pRF,2)
write.csv(var6,file = 'antimicrobial_toxin/rf_cluster/importance_cluster.csv')
