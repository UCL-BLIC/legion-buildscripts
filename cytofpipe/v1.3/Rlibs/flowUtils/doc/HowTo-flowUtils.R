### R code from vignette source 'HowTo-flowUtils.Rnw'

###################################################
### code chunk number 1: loadPackage
###################################################
library("flowUtils")
library("flowCore")
options(width=60)


###################################################
### code chunk number 2: ReadGatingML1
###################################################
gateFile <- system.file("extdata", "GatingML2.0_Example1.xml", 
    package = "flowUtils")
flowEnv  <- new.env()
read.gatingML(gateFile, flowEnv)
for (x in ls(flowEnv)) 
    if (is(flowEnv[[x]], "filter"))
        cat(paste("Gate", x, "of class", class(flowEnv[[x]]), "\n"))


###################################################
### code chunk number 3: ReadGatingML2
###################################################
flowEnv1.5  <- new.env()
g1.5Example <- system.file("extdata/Gating-MLFiles", "01Rectangular.xml",
  package="gatingMLData")
read.gatingML(g1.5Example, flowEnv1.5)
ls(flowEnv1.5)

flowEnv2.0  <- new.env()
g2.0Example <- system.file("extdata/Gml2/Gating-MLFiles",
  "gates3.xml", package = "gatingMLData")
read.gatingML(g2.0Example, flowEnv2.0)
ls(flowEnv2.0)


###################################################
### code chunk number 4: ReadGatingML3
###################################################
flowEnv2.0[['myRectangleGate4LogicleArcSinHFCSCompensated']]


###################################################
### code chunk number 5: ReadGatingML4
###################################################
str(flowEnv2.0[['myRectangleGate4LogicleArcSinHFCSCompensated']])


###################################################
### code chunk number 6: ReadGatingML5
###################################################
str(flowEnv2.0[['myLogicle.FCS.PE-A']])


###################################################
### code chunk number 7: ReadGatingML6
###################################################
str(flowEnv2.0[['Tr_Arcsinh.FCS.APC-Cy7-A']])


###################################################
### code chunk number 8: ReadGatingML7
###################################################
str(flowEnv2.0[['myPolygonGateWithCustomSpillover']])
flowEnv2.0[['MySpill']]


###################################################
### code chunk number 9: ApplyGatingML1
###################################################
fcsFile <- system.file("extdata/Gml2/FCSFiles", "data1.fcs", 
    package = "gatingMLData")
myFrame <- read.FCS(fcsFile, transformation="linearize-with-PnG-scaling")
for (x in ls(flowEnv)) if (is(flowEnv[[x]], "filter")) {
    result <- filter(myFrame, flowEnv[[x]])
    print(summary(result))
}


###################################################
### code chunk number 10: WriteGatingML1
###################################################
flowEnv <- new.env()
flowEnv[['myGate']] <- rectangleGate(filterId="myGate", 
    list("FSC-H"=c(150, 300), "SSC-H"=c(200, 600)))
outputFile <- tempfile(fileext=".gating-ml2.xml")
write.gatingML(flowEnv, outputFile)


###################################################
### code chunk number 11: WriteGatingML2
###################################################
flowEnv <- new.env()
covM <- matrix(c(62.5, 37.5, 37.5, 62.5), nrow = 2, byrow=TRUE)
colnames(covM) <- c("FL1-H", "FL2-H")
compPars <- list(
  compensatedParameter(parameters="FL1-H", spillRefId="SpillFromFCS", 
    transformationId=paste("FL1-H", "_compensated_according_to_FCS", sep=""), 
    searchEnv=flowEnv),
  compensatedParameter(parameters="FL2-H", spillRefId="SpillFromFCS", 
    transformationId=paste("FL2-H", "_compensated_according_to_FCS", sep=""), 
    searchEnv=flowEnv)
)
myEl <- ellipsoidGate(mean=c(12, 16), distance=1, .gate=covM, filterId="myEl")
myEl@parameters <- new("parameters", compPars)
flowEnv[['myEl']] <- myEl
write.gatingML(flowEnv)


###################################################
### code chunk number 12: WriteGatingML3
###################################################
spillM <- matrix(c(1, 0.03, 0.07, 1), nrow = 2, byrow=TRUE)
colnames(spillM) <- c("FL1-H", "FL2-H")
rownames(spillM) <- c("Comp-FL1-H", "Comp-FL2-H")
pars <- new("parameters", list("FL1-H", "FL2-H"))
myComp <- compensation(spillover=spillM, compensationId='myComp', pars)
flowEnv[['myComp']] <- myComp
compPars <- list(
  compensatedParameter(parameters="FL1-H", spillRefId="myComp", 
    transformationId="Comp-FL1-H", searchEnv=flowEnv),
  compensatedParameter(parameters="FL2-H", spillRefId="myComp", 
    transformationId="Comp-FL2-H", searchEnv=flowEnv)
)
myEl@parameters <- new("parameters", compPars)
flowEnv[['myEl']] <- myEl
write.gatingML(flowEnv)


###################################################
### code chunk number 13: WriteGatingML4
###################################################
flowEnv <- new.env()
specM <- matrix(c(0.78, 0.13, 0.22, 0.05, 0.57, 0.89), nrow = 2, byrow=TRUE)
colnames(specM) <- c("FL1-H", "FL2-H", "FL3-H")
rownames(specM) <- c("Deconvoluted-P1", "Deconvoluted-P2")
pars <- new("parameters", list("FL1-H", "FL2-H", "FL3-H"))
mySpecM <- compensation(spillover=specM, compensationId='specM', pars)
flowEnv[['mySpecM']] <- mySpecM
compPars <- list(
  compensatedParameter(parameters="FL1-H", spillRefId="mySpecM", 
    transformationId="Deconvoluted-P1", searchEnv=flowEnv),
  compensatedParameter(parameters="FL2-H", spillRefId="mySpecM", 
    transformationId="Deconvoluted-P2", searchEnv=flowEnv)
)
myEl@parameters <- new("parameters", compPars)
flowEnv[['myEl']] <- myEl
write.gatingML(flowEnv)


###################################################
### code chunk number 14: WriteGatingML5
###################################################
flowEnv <- new.env()
myTrQuad <- quadGate(filterId = "myTrQuad", "APC-A" = 0.5, "APC-Cy7-A" = 0.6)
trArcSinH1 <- asinhtGml2(parameters = "APC-A", T = 1000, M = 4.5, A = 0, 
  transformationId="trArcSinH1")
trLogicle1 <- logicletGml2(parameters = "APC-Cy7-A", T = 1000, W = 0.5, 
  M = 4.5, A = 0, transformationId="trLogicle1")
flowEnv[['trArcSinH1']] <- trArcSinH1
flowEnv[['trLogicle1']] <- trLogicle1
trPars <- list(
  transformReference("trArcSinH1", flowEnv),
  transformReference("trLogicle1", flowEnv)
)
myTrQuad@parameters <- new("parameters", trPars)
flowEnv[['myTrQuad']] <- myTrQuad
write.gatingML(flowEnv)


###################################################
### code chunk number 15: WriteGatingML5Bound (eval = FALSE)
###################################################
## trArcSinH1 <- asinhtGml2(parameters = "APC-A", T = 1000, M = 4.5, A = 0, 
##   transformationId="trArcSinH1", boundMin = 0.02, boundMax = 0.96)
## trLogicle1 <- logicletGml2(parameters = "APC-Cy7-A", T = 1000, W = 0.5, 
##   M = 4.5, A = 0, transformationId="trLogicle1", boundMin = -0.04)


###################################################
### code chunk number 16: WriteGatingML6
###################################################
rm(list=ls(flowEnv), envir=flowEnv)
trArcSinH1@parameters <- compensatedParameter(parameters="APC-A", 
  spillRefId="SpillFromFCS", searchEnv=flowEnv,
  transformationId= "APC-A_compensated_according_to_FCS")
trLogicle1@parameters <- compensatedParameter(parameters="APC-Cy7-A", 
  spillRefId="SpillFromFCS", searchEnv=flowEnv,
  transformationId="APC-Cy7-A_compensated_according_to_FCS")
trPars <- list(trArcSinH1,trLogicle1)
myTrQuad@parameters <- new("parameters", trPars)
flowEnv[['myTrQuad']] <- myTrQuad
write.gatingML(flowEnv)


###################################################
### code chunk number 17: WriteGatingML7
###################################################
flowEnv <- new.env()
rat1 <- ratio("FSC-A", "SSC-A", transformationId = "rat1")
myRectGate <- rectangleGate(filterId="myRectGate", "rat1"=c(0.8, 1.4))
myRectGate@parameters <- new("parameters", list(rat1))
flowEnv[['myRectGate']] <- myRectGate
write.gatingML(flowEnv)


###################################################
### code chunk number 18: WriteGatingML8
###################################################
flowEnv <- new.env()
myASinH <- asinht("FL3-W", a = 1.5828, b = 0.0965, transformationId = "myASinH")
gate1 <- rectangleGate(filterId="gate1", "myASinH"=c(0.3, 0.7))
gate1@parameters <- new("parameters", list(myASinH))
flowEnv[['gate1']] <- gate1
write.gatingML(flowEnv)


###################################################
### code chunk number 19: WriteGatingML9
###################################################
flowEnv <- new.env()
# Creation of a simplified spillover matrix
spillM <- matrix(c(1, 0, 0.03, 0, 0, 1, 0, 0.07, 0.1, 0, 1, 0, 0, 0.05, 0, 1), 
    nrow = 4, byrow=TRUE)
colnames(spillM) <- c("FL1-A", "FL1-W", "FL2-A", "FL2-W")
rownames(spillM) <- c("cFL1-A", "cFL1-W", "cFL2-A", "cFL2-W")
pars <- new("parameters", list("FL1-A", "FL1-W", "FL2-A", "FL2-W"))
myComp <- compensation(spillover=spillM, compensationId='myComp', pars)
flowEnv[['myComp']] <- myComp
myComp


###################################################
### code chunk number 20: WriteGatingML10
###################################################
# First dimension is a log(cFL1-A / cFL1-W)  
myRatio <- ratio("FL1-A", "FL1-W", transformationId = "myRatio")
myRatio@numerator <- compensatedParameter(parameters="FL1-A", 
    spillRefId="myComp", transformationId="cFL1-A", searchEnv=flowEnv)
myRatio@denominator <- compensatedParameter(parameters="FL1-W", 
    spillRefId="myComp", transformationId="cFL1-W", searchEnv=flowEnv)
myLog <- logtGml2(myRatio, T = 1, M = 1, transformationId="myLog")
# Second dimension is a Hyperlog(cFL2-A)
secPar <- compensatedParameter(parameters="FL2-A", spillRefId="myComp", 
    transformationId="cFL2-A", searchEnv=flowEnv)
myHLog <- hyperlogtGml2(secPar, T=262144, M=4.5, W=0.5, A=0, "myHLog")
# A Polygon gate in the two defined dimensions
vertices <- matrix(c(0.9, 0.5, 1.2, 0.6, 1.1, 0.8), nrow=3, ncol=2, byrow=TRUE)
myGate <- polygonGate(filterId="myGate", .gate=vertices, 
    new("parameters", list(myLog, myHLog)))
flowEnv[['myGate']] <- myGate
# Finally, write the Gating-ML output
write.gatingML(flowEnv)


###################################################
### code chunk number 21: WriteGatingML11
###################################################
flowEnv <- new.env()
logicle1 <- logicletGml2(parameters="FL1-H", T=10000, M=4.5, A=0, W=.5, "logicle1")
logicle2 <- logicletGml2(parameters="FL2-H", T=10000, M=4.5, A=0, W=.5, "logicle2")
lin1 <- lintGml2(parameters = "FL1-H", T = 10000, A = 0, "lin1")
lin2 <- lintGml2(parameters = "FL2-H", T = 10000, A = 0, "lin2")
rectG <- rectangleGate(filterId="rectG", "logicle1"=c(.1, .6), "lin2"=c(.2, .6))
rectG@parameters <- new("parameters", list(logicle1, lin2))
rangeG1 <- rectangleGate(filterId="rangeG1", "logicle2"=c(0.1, 0.5))
rangeG1@parameters <- new("parameters", list(logicle2))
rangeG2 <- rectangleGate(filterId="rangeG2", "lin1"=c(0.6, 0.9))
rangeG2@parameters <- new("parameters", list(lin1))
flowEnv[['rectG']] <- rectG
flowEnv[['rangeG1']] <- rangeG1
flowEnv[['rangeG2']] <- rangeG2
write.gatingML(flowEnv)


###################################################
### code chunk number 22: WriteGatingMLUnsupportedCase1
###################################################
logicle1 <- logicletGml2(parameters = "FL1-H", T = 1000, M = 4.5, A = 0, 
  W = 0.5, transformationId="logicle1")
logicle2 <- logicletGml2(parameters = "logicle1", T = 1000, M = 4.5, A = 0, 
  W = 0.5, transformationId="logicle2")
logicle2@parameters <- logicle1
myRect <- rectangleGate(filterId="myRect", list("logicle2"=c(0, .6)))
myRect@parameters <- new("parameters", list(logicle2))
flowEnv[['myRect']] <- myRect
x <- tryCatch(write.gatingML(flowEnv), error = function(e) { e })
x$message


###################################################
### code chunk number 23: WriteGatingMLUnsupportedCase2
###################################################
flowEnv <- new.env()
tSS <- splitscale(parameters = "FL1-H", r = 1024, maxValue = 10000, 
    transitionChannel = 256, transformationId = "tSS")
myRect <- rectangleGate(filterId="myRect", list("tss"=c(100, 700)))
myRect@parameters <- new("parameters", list(tSS))
flowEnv[['myRect']] <- myRect
x <- tryCatch(write.gatingML(flowEnv), error = function(e) { e })
x$message


###################################################
### code chunk number 24: WriteGatingMLUnsupportedCase3
###################################################
flowEnv <- new.env()
# Gating-ML 1.5 example 5.3.4.c
a <- matrix(c(-1, 0, 0, 0, -1, 0, 0, 0, -1, 1, 0, 0, 0, 0, 1), ncol=3)
b <- c(100, 50, 0, 250, 300)
myPolytope = polytopeGate(filterId='myPolytope', .gate=a, b=b, 
  list("FSC-H", "SSC-H", "FL1-H"))
flowEnv[['myPolytope']] <- myPolytope
x <- tryCatch(write.gatingML(flowEnv), error = function(e) { e })
x$message


###################################################
### code chunk number 25: WriteGatingMLUnsupportedCase4
###################################################
flowEnv <- new.env()
myNorm2Filter <- norm2Filter("FSC-H", "SSC-H", filterId="myNorm2Filter")
flowEnv[['myNorm2Filter']] <- myNorm2Filter
x <- tryCatch(write.gatingML(flowEnv), error = function(e) { e })
x$message


###################################################
### code chunk number 26: TestGatingMLCompliance1 (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite("gatingMLData")


###################################################
### code chunk number 27: TestGatingMLCompliance2 (eval = FALSE)
###################################################
## testGatingMLCompliance("ComplianceReport_v1.5.html", version=1.5)


###################################################
### code chunk number 28: TestGatingMLCompliance3 (eval = FALSE)
###################################################
## testGatingMLCompliance("ComplianceReport_v2.0.html", version=2.0)


###################################################
### code chunk number 29: PrecisionIssuesDemo
###################################################
fcsFile  <- system.file("extdata/examples", "166889.fcs", 
    package = "gatingMLData")
gateFile <- system.file("extdata/examples", "GatingML2.0_Export_166889.xml", 
    package = "gatingMLData")
myFrame  <- read.FCS(fcsFile, transformation="linearize-with-PnG-scaling")
flowEnv  <- new.env()
read.gatingML(gateFile, flowEnv)
for (x in ls(flowEnv)) if (is(flowEnv[[x]], "filter")) {
    result <- filter(myFrame, flowEnv[[x]])
    print(summary(result))
}


###################################################
### code chunk number 30: PrecisionIssuesDemo2
###################################################
class(flowEnv[['Gate_1_UEUtQSBTU0MtQSBFMQ.._UEUtQQ.._U1NDLUE.']])
str(flowEnv[['Gate_1_UEUtQSBTU0MtQSBFMQ.._UEUtQQ.._U1NDLUE.']])
class(flowEnv[['GateSet_1_UEUtQQ.._U1NDLUE.']])
str(flowEnv[['GateSet_1_UEUtQQ.._U1NDLUE.']])


###################################################
### code chunk number 31: PrecisionIssuesDemo3
###################################################
cat(readChar(gateFile, file.info(gateFile)$size))


