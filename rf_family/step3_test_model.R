#step3:test model###############################################################

rm(list = ls())
library(caret)
library(randomForest)
library(ggplot2)
load(file = 'antimicrobial_toxin/rf_family/optvar51_family.Rdata')

set.seed(51)
dat <- read.csv(file = 'antimicrobial_toxin/train_new.csv',header = T)
dat2 <- na.omit(dat)
dat3 <- dat2[,-1:-4]
x <- dat3[,optvar]
logBBB <- dat2$toxin_families
load(file ="antimicrobial_toxin/rf_family/inTrain_family.Rdata")
trainx <- x[inTrain, ]
testx <- x[-inTrain, ]
trainy <- logBBB[inTrain]
testy <- logBBB[-inTrain]

load(file = 'antimicrobial_toxin/rf_family/rf_gridsearch_family.Rdata')

testset <- testx
predictRF <- predict(rf_gridsearch,newdata = testset)

R2(testy, predictRF)
rsq <- function(x, y) summary(lm(y~x))$r.squared
rsq(testy, predictRF)

plot <- data.frame(predictRF,testy)
ggplot(data = plot, aes(x = testy, y = predictRF)) + 
  geom_point(size = 4,fill='red',color="black",alpha=0.4,shape=21) +
  labs(x = "Observed antimicrobial toxin families", y = "Predicted antimicrobial toxin families") +
  theme_classic() +
  theme(axis.text=element_text(size=12,color = 'black'),axis.title=element_text(size=14,color = 'black'),
        axis.line.x=element_line(linetype=1,color="black",size=0),
        axis.line.y=element_line(linetype=1,color="black",size=0),
        axis.ticks.x=element_line(color="black",size=1,lineend = 2),
        axis.ticks.y=element_line(color="black",size=1,lineend = 2),
        panel.border = element_rect(linetype=1,color = "black", size = 1.5, fill = NA)) +
  geom_smooth(formula = y ~ x, method = "lm",alpha=0.3,linetype = 0 ) +
  stat_smooth (formula = y ~ x,method = 'lm',geom="line",alpha=0.5, color="#9EC3DB",size=1.5)
