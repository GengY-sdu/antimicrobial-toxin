###############
rm(list = ls())
library(caret)
library(ggplot2)

set.seed(51)
dat <- read.csv(file = 'train_new.csv',header = T)
dat2 <- na.omit(dat)
dat3 <- dat2[,-1:-4]
x <- dat3[,-nearZeroVar(dat3)]   
logBBB <- dat2$toxin_xxx
inTrain <- createDataPartition(logBBB, p = 0.8, list = FALSE)[,1]  
save(inTrain,file = 'inTrain.Rdata')
trainx <- x[inTrain, ]
testx <- x[-inTrain, ]
trainy <- logBBB[inTrain]
testy <- logBBB[-inTrain]

set.seed(51)
control<-rfeControl(functions = rfFuncs,method = "cv",number = 10,rerank = T) 
rfFuncs <- rfe(trainx, trainy,sizes = c(1:97),rfeControl = control,metric = "Rsquared")
save(rfFuncs,file = 'rfFuncs.Rdata')
variable <- rfFuncs$results$Variables
rmse <- rfFuncs$results$RMSE
mae <- rfFuncs$results$MAE
rsq <- rfFuncs$results$Rsquared
sum <- cbind(variable,rmse,mae,rsq)
optvar <- rfFuncs$optVariables
save(optvar,file = 'optvar51.Rdata')

plot <- data.frame(rfFuncs$results$Rsquared,rfFuncs$results$Variables)
ggplot(data = plot, aes(x = rfFuncs.results.Variables, y = rfFuncs.results.Rsquared)) + 
  geom_point(size = 4,color="red",alpha=0.5) + 
  geom_line(size = 1,color='black',alpha=0.5) + 
  labs(x = "Number of variables", y = expression("10-fold Cross-Validation"~R^2)) +
  theme_classic() +
  theme(axis.text=element_text(size=12,color = 'black'),axis.title=element_text(size=14,color = 'black'),
        axis.line.x=element_line(linetype=1,color="black",size=0),
        axis.line.y=element_line(linetype=1,color="black",size=0),
        axis.ticks.x=element_line(color="black",size=1,lineend = 2),
        axis.ticks.y=element_line(color="black",size=1,lineend = 2),
        panel.border = element_rect(linetype=1,color = "black", size = 1.5, fill = NA))

################
rm(list = ls())
library(caret)
library(randomForest)
library(ggplot2)
load(file = 'optvar51.Rdata')

set.seed(51)
dat <- read.csv(file = 'train_new.csv',header = T)
dat2 <- na.omit(dat)
dat3 <- dat2[,-1:-4]
x <- dat3[,optvar]
logBBB <- dat2$toxin_xxx
inTrain <- createDataPartition(logBBB, p = 0.8, list = FALSE)[,1]
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

set.seed(51)
control_grid <- trainControl(method="cv", number=10, search="grid",savePredictions = 'final')
tunegrid <- expand.grid(mtry=c(1,5,10,15,20,25,30,35,40,45,50,55,60),
                        ntree = c(500,800,1000,1200,1500,2000,2500),
                        nodesize=c(4:16))
rf_gridsearch <- train(trainx,trainy, method=customRF,
                       tuneGrid=tunegrid, trControl=control_grid, 
                       metric = "Rsquared", importance=T)

print(rf_gridsearch)
save(rf_gridsearch,file = 'rf_gridsearch.Rdata')

plot(rf_gridsearch,ylab=expression("10-fold Cross-Validation"~R^2),
     col=c('#F26F66','#8572AD','#95D1C0','#9EC3DB','#CC9C5A','#CECCE6',
           '#F47D38','#C0C0C0','#F1B719','#97CD79','#5DB7E8','#CCAF47','#BFDED6'))

plot <- data.frame(rf_gridsearch$pred$pred,rf_gridsearch$pred$obs)
ggplot(data = plot, aes(x = rf_gridsearch.pred.obs, y = rf_gridsearch.pred.pred)) + 
  geom_point(size = 4,fill='red',color="black",alpha=0.4,shape=21) +
  labs(x = "Observed antimicrobial toxin xxx", y = "Predicted antimicrobial toxin xxx") +
  theme_classic() +
  theme(axis.text=element_text(size=12,color = 'black'),axis.title=element_text(size=14,color = 'black'),
        axis.line.x=element_line(linetype=1,color="black",size=0),
        axis.line.y=element_line(linetype=1,color="black",size=0),
        axis.ticks.x=element_line(color="black",size=1,lineend = 2),
        axis.ticks.y=element_line(color="black",size=1,lineend = 2),
        panel.border = element_rect(linetype=1,color = "black", size = 1.5, fill = NA)) +
  geom_smooth(formula = y ~ x, method = "lm",alpha=0.3,linetype = 0 ) +
  stat_smooth (formula = y ~ x,method = 'lm',geom="line",alpha=0.5, color="#9EC3DB",size=1.5)

############################

testset <- testx
predictRF <- predict(rf_gridsearch,newdata = testset)
predictRF1 <- predict(rf_gridsearch$finalModel,newdata = testset)

R2(testy, predictRF)
rsq <- function(x, y) summary(lm(y~x))$r.squared
rsq(testy, predictRF)

plot <- data.frame(predictRF,testy)
ggplot(data = plot, aes(x = testy, y = predictRF)) + 
  geom_point(size = 4,fill='red',color="black",alpha=0.4,shape=21) +
  labs(x = "Observed antimicrobial toxin xxx", y = "Predicted antimicrobial toxin xxx") +
  theme_classic() +
  theme(axis.text=element_text(size=12,color = 'black'),axis.title=element_text(size=14,color = 'black'),
        axis.line.x=element_line(linetype=1,color="black",size=0),
        axis.line.y=element_line(linetype=1,color="black",size=0),
        axis.ticks.x=element_line(color="black",size=1,lineend = 2),
        axis.ticks.y=element_line(color="black",size=1,lineend = 2),
        panel.border = element_rect(linetype=1,color = "black", size = 1.5, fill = NA)) +
  geom_smooth(formula = y ~ x, method = "lm",alpha=0.3,linetype = 0 ) +
  stat_smooth (formula = y ~ x,method = 'lm',geom="line",alpha=0.5, color="#9EC3DB",size=1.5)

############################

set.seed(x)
control_grid <- trainControl(method = "cv", number = 10, search = "grid")
rfmodelx <- train(trainx,trainy, method=customRF,
                 tuneGrid=rf_gridsearch$bestTune, trControl=control_grid, 
                 metric = "Rsquared",impotance=T)
rfmodelx
predictRF_x <- predict(rfmodelx,newdata = testset)
rsq(testy, predictRF_x)

mapset <- read.csv(file = 'map.csv')
predictRF1 <- predict(rfmodelx,newdata = mapset)
predictRF2 <- predict(rfmodel2,newdata = mapset)
predictRF3 <- predict(rfmodel3,newdata = mapset)
predictRF4 <- predict(rfmodel4,newdata = mapset)
predictRF5 <- predict(rfmodel5,newdata = mapset)
predictRF6 <- predict(rfmodel6,newdata = mapset)
predictRF7 <- predict(rfmodel7,newdata = mapset)
predictRF8 <- predict(rfmodel8,newdata = mapset)
predictRF9 <- predict(rfmodel9,newdata = mapset)
predictRF10 <- predict(rfmodel10,newdata = mapset)

##########################################################
predictsum <- data.frame(predictRF1,predictRF2,predictRF3,predictRF4,predictRF5,
                         predictRF6,predictRF7,predictRF8,predictRF9,predictRF10,
                         mapset$longitude,mapset$latitude)

for (i in 1:nrow(predictsum)){
  predictsum$prediction_sd[i] <- sd(as.numeric(predictsum[i,1:10]))
  predictsum$prediction_mean[i] <- mean(as.numeric(predictsum[i,1:10]))
}
predictsum$cv <- predictsum$prediction_sd / predictsum$prediction_mean
write.csv(predictsum,file = 'final_xxx.csv')

########################################################
var6 <- varImp(rfmodel6, scale=FALSE)$importance  
var6$pRF <- var6$Overall / sum(var6$Overall) * 100
var6$rfround <- round(var6$pRF,2)

write.csv(var6,file = 'importance.csv')
