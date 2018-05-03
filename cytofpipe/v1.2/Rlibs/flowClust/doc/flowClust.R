### R code from vignette source 'flowClust.Rnw'

###################################################
### code chunk number 1: stage0
###################################################
library(flowCore)
library(flowClust)


###################################################
### code chunk number 2: stage1
###################################################
data(rituximab)
summary(rituximab)
res1 <- flowClust(rituximab, varNames=c("FSC.H", "SSC.H"), K=1, B=100)


###################################################
### code chunk number 3: stage2
###################################################
rituximab2 <- rituximab[rituximab %in% res1,]
res2 <- flowClust(rituximab2, varNames=c("FL1.H", "FL3.H"), K=1:6, B=100)


###################################################
### code chunk number 4: stage2Result
###################################################
criterion(res2, "BIC")
summary(res2[[3]])


###################################################
### code chunk number 5: stage2ChangeRule1
###################################################
ruleOutliers(res2[[3]]) <- list(level=0.95)
summary(res2[[3]])


###################################################
### code chunk number 6: stage2ChangeRule2
###################################################
ruleOutliers(res2[[3]]) <- list(z.cutoff=0.6)
summary(res2[[3]])


###################################################
### code chunk number 7: stage2Alternative
###################################################
flowClust(rituximab2, varNames=c("FL1.H", "FL3.H"), K=2, B=100, min=c(0,0), max=c(400,800))


###################################################
### code chunk number 8: stage2Scatter
###################################################
plot(res2[[3]], data=rituximab2, level=0.8, z.cutoff=0)


###################################################
### code chunk number 9: stage2Contour
###################################################
res2.den <- density(res2[[3]], data=rituximab2)
plot(res2.den)


###################################################
### code chunk number 10: stage2Image
###################################################
plot(res2.den, type="image")


###################################################
### code chunk number 11: stage2Hist
###################################################
hist(res2[[3]], data=rituximab2, subset="FL1.H")


###################################################
### code chunk number 12: stage2Hist
###################################################
hist(res2[[3]], data=rituximab2, subset=1)


###################################################
### code chunk number 13: stage2f
###################################################
s2filter <- tmixFilter("s2filter", c("FL1.H", "FL3.H"), K=3, B=100)
res2f <- filter(rituximab2, s2filter)


###################################################
### code chunk number 14: stage2fSubsetting
###################################################
Subset(rituximab2, res2f)
split(rituximab2, res2f, population=list(sc1=1:2, sc2=3))


###################################################
### code chunk number 15: stage2fRectGate
###################################################
rectGate <- rectangleGate(filterId="rectRegion", "FL1.H"=c(0, 400), "FL3.H"=c(0, 800))
MBCfilter <- tmixFilter("MBCfilter", c("FL1.H", "FL3.H"), K=2, B=100)
filter(rituximab2, MBCfilter %subset% rectGate)


###################################################
### code chunk number 16: priors
###################################################
prior<-flowClust2Prior(res2[[2]],kappa=1,Nt=5000)
prior2<-prior
prior2$Mu0[1,]<-rep(box(200,prior2$lambda),2)
prior2$Lambda0<-prior2$Lambda0/2
pfit2<-flowClust(rituximab2,varNames=c("FL1.H","FL3.H"),K=2,prior=prior2,usePrior="yes")

par(mfrow=c(1,2))
plot(res2[[2]],data=rituximab2)
plot(pfit2,data=rituximab2)


