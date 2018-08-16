### R code from vignette source 'ncdfFlow.Rnw'

###################################################
### code chunk number 1: loadPackage
###################################################
library(ncdfFlow)


###################################################
### code chunk number 2: ncdfFlowSet
###################################################
path<-system.file("extdata","compdata","data",package="flowCore")
files<-list.files(path,full.names=TRUE)[1:3]
nc1 <- read.ncdfFlowSet(files=files)
nc1


###################################################
### code chunk number 3: flowSet
###################################################
fs1  <- read.flowSet(files=files)


###################################################
### code chunk number 4: ncdfFlowSet
###################################################
unlink(nc1)
rm(nc1)


###################################################
### code chunk number 5: ncdfFlowSet
###################################################
nc1  <- read.ncdfFlowSet(files=files, isWriteSlice= FALSE)
nc1[[1]]


###################################################
### code chunk number 6: ncdfFlowSet
###################################################
targetSampleName<-sampleNames(fs1)[1]
nc1[[targetSampleName]] <- fs1[[1]]
nc1[[1]]
nc1[[2]]


###################################################
### code chunk number 7: ncdfFlowSet
###################################################
nc2 <- clone.ncdfFlowSet(nc1, isEmpty = TRUE)
nc2[[1]]
nc2[[sampleNames(fs1)[1]]] <- fs1[[1]]
nc2[[1]]


###################################################
### code chunk number 8: ncdfFlowSet
###################################################
unlink(nc2)
rm(nc2)

unlink(nc1)
rm(nc1)


###################################################
### code chunk number 9: ncdfFlowSet
###################################################
data(GvHD)
GvHD <- GvHD[pData(GvHD)$Patient %in% 6:7][1:4]
nc1<-ncdfFlowSet(GvHD)


###################################################
### code chunk number 10: ncdfFlowSet
###################################################
fs1<-as.flowSet(nc1,top=2)


###################################################
### code chunk number 11: metaData
###################################################
phenoData(nc1)
pData(nc1)
varLabels(nc1)
varMetadata(nc1)
sampleNames(nc1)
keyword(nc1,"FILENAME")
identifier(nc1)
colnames(nc1)
colnames(nc1,prefix="s6a01")
length(nc1)
getIndices(nc1,"s6a01")


###################################################
### code chunk number 12: extraction
###################################################

nm<-sampleNames(nc1)[1]
expr1<-paste("nc1$'",nm,"'",sep="")
eval(parse(text=expr1))
nc1[[nm]]

nm<-sampleNames(nc1)[c(1,3)]
nc2<-nc1[nm]
summary(nc2)


###################################################
### code chunk number 13: fsApply
###################################################
fsApply(nc1,range)
fsApply(nc1, each_col, median)


###################################################
### code chunk number 14: Transformation and compensation
###################################################
cfile <- system.file("extdata","compdata","compmatrix", package="flowCore")
comp.mat <- read.table(cfile, header=TRUE, skip=2, check.names = FALSE)
comp <- compensation(comp.mat)

#compensation
summary(nc1)[[1]]
nc2<-compensate(nc1, comp)
summary(nc2)[[1]]
unlink(nc2)
rm(nc2)

#transformation
asinhTrans <- arcsinhTransform(transformationId="ln-transformation", a=1, b=1, c=1)
nc2 <- transform(nc1,`FL1-H`=asinhTrans(`FL1-H`))
summary(nc1)[[1]]
summary(nc2)[[1]]
unlink(nc2)
rm(nc2)


###################################################
### code chunk number 15: Gating
###################################################
rectGate <- rectangleGate(filterId="nonDebris","FSC-H"=c(200,Inf))
fr <- filter(nc1,rectGate)
summary(fr)

rg2 <- rectangleGate(filterId="nonDebris","FSC-H"=c(300,Inf))
rg3 <- rectangleGate(filterId="nonDebris","FSC-H"=c(400,Inf))
flist <- list(rectGate, rg2, rg3)
names(flist) <- sampleNames(nc1[1:3])
fr3 <- filter(nc1[1:3], flist)
summary(fr3[[3]])


###################################################
### code chunk number 16: Subsetting
###################################################
nc2 <- Subset(nc1,rectGate)
summary(nc2[[1]])
nc2 <- Subset(nc1, fr)
summary(nc2[[1]])
rm(nc2)


morphGate <- norm2Filter("FSC-H", "SSC-H", filterId = "MorphologyGate",scale = 2)
smaller <- Subset(nc1[c(1,3)], morphGate,c("FSC-H", "SSC-H"))
smaller[[1]]
nc1[[1]]
rm(smaller)


###################################################
### code chunk number 17: split
###################################################

##splitting by a gate
qGate <- quadGate(filterId="qg", "FSC-H"=200, "SSC-H"=400)
fr<-filter(nc1,qGate)
ncList<-split(nc1,fr)
ncList
nc1[[1]]
ncList[[2]][[1]]
ncList[[1]][[1]]


###################################################
### code chunk number 18: split
###################################################
ncList_new<-split(nc1,fr,isNew=TRUE)


