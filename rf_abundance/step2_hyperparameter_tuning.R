#step2:hyperparameter tuning####################################################

rm(list = ls())
library(caret)
library(randomForest)
library(ggplot2)
load(file = 'antimicrobial_toxin/rf_abundance/optvar51_abundance.Rdata')

set.seed(51)
dat <- read.csv(file = 'antimicrobial_toxin/train_new.csv',header = T)
dat2 <- na.omit(dat)
dat3 <- dat2[,-1:-4]
x <- dat3[,optvar]
logBBB <- dat2$toxin_abundance
load(file ="antimicrobial_toxin/rf_abundance/inTrain_abundance.Rdata")
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
save(rf_gridsearch,file = 'antimicrobial_toxin/rf_abundance/rf_gridsearch_abundance.Rdata')

plot(rf_gridsearch,ylab=expression("10-fold Cross-Validation"~R^2),
     col=c('#F26F66','#8572AD','#95D1C0','#9EC3DB','#CC9C5A','#CECCE6',
           '#F47D38','#C0C0C0','#F1B719','#97CD79','#5DB7E8','#CCAF47','#BFDED6'))

plot <- data.frame(rf_gridsearch$pred$pred,rf_gridsearch$pred$obs)
ggplot(data = plot, aes(x = rf_gridsearch.pred.obs, y = rf_gridsearch.pred.pred)) + 
  geom_point(size = 4,fill='red',color="black",alpha=0.4,shape=21) +
  labs(x = "Observed antimicrobial toxin abundance", y = "Predicted antimicrobial toxin abundance") +
  theme_classic() +
  theme(axis.text=element_text(size=12,color = 'black'),axis.title=element_text(size=14,color = 'black'),
        axis.line.x=element_line(linetype=1,color="black",size=0),
        axis.line.y=element_line(linetype=1,color="black",size=0),
        axis.ticks.x=element_line(color="black",size=1,lineend = 2),
        axis.ticks.y=element_line(color="black",size=1,lineend = 2),
        panel.border = element_rect(linetype=1,color = "black", size = 1.5, fill = NA)) +
  geom_smooth(formula = y ~ x, method = "lm",alpha=0.3,linetype = 0 ) +
  stat_smooth (formula = y ~ x,method = 'lm',geom="line",alpha=0.5, color="#9EC3DB",size=1.5)
