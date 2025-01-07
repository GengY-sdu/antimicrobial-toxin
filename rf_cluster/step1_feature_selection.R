#step1:feature selection########################################################

rm(list = ls())
library(caret)
library(ggplot2)

set.seed(51)
dat <- read.csv(file = 'antimicrobial_toxin/train_new.csv',header = T)
dat2 <- na.omit(dat)
dat3 <- dat2[,-1:-4]
x <- dat3[,-nearZeroVar(dat3)]
logBBB <- dat2$toxin_clusters
inTrain <- createDataPartition(logBBB, p = 0.8, list = FALSE)[,1]
save(inTrain,file = 'antimicrobial_toxin/rf_cluster/inTrain_cluster.Rdata')
trainx <- x[inTrain, ]
testx <- x[-inTrain, ]
trainy <- logBBB[inTrain]
testy <- logBBB[-inTrain]

set.seed(51)
control<-rfeControl(functions = rfFuncs,method = "cv",number = 10,rerank = T)  
rfFuncs <- rfe(trainx, trainy,sizes = c(1:97),rfeControl = control,metric = "Rsquared")
save(rfFuncs,file = 'antimicrobial_toxin/rf_cluster/rfFuncs_cluster.Rdata')
variable <- rfFuncs$results$Variables
rmse <- rfFuncs$results$RMSE
mae <- rfFuncs$results$MAE
rsq <- rfFuncs$results$Rsquared
sum <- cbind(variable,rmse,mae,rsq)
optvar <- rfFuncs$optVariables
save(optvar,file = 'antimicrobial_toxin/rf_cluster/optvar51_cluster.Rdata')

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
