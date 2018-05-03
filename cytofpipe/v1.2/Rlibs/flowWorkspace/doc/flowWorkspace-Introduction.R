## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, results = "markup", message = FALSE)

## ----findxml-------------------------------------------------------------
library(flowWorkspace)
path <- system.file("extdata",package="flowWorkspaceData");
wsfile <- list.files(path, pattern="A2004Analysis.xml", full = TRUE)

## ----openws,results='markup'---------------------------------------------
ws <- openWorkspace(wsfile)
ws

## ----getsamples-ws-------------------------------------------------------
getSamples(ws)

## ----getgroups-----------------------------------------------------------
getSampleGroups(ws)

## ------------------------------------------------------------------------
sn <- "a2004_O1T2pb05i_A1_A01.fcs"
getKeywords(ws, sn)[1:5]

## ----parsews,message=FALSE-----------------------------------------------
gs <- parseWorkspace(ws,name = 1); #import the first group
#Lots of output here suppressed for the vignette.
gs

## ----sampleNames---------------------------------------------------------
sampleNames(gs)

## ----parseGatingML, eval=FALSE-------------------------------------------
#  library(CytoML)
#  xmlfile <- system.file("extdata/cytotrol_tcell_cytobank.xml", package = "CytoML")
#  fcsFiles <- list.files(pattern = "CytoTrol", system.file("extdata", package = "flowWorkspaceData"), full = T)
#  gs1 <- parse.gatingML(xmlfile, fcsFiles)

## ----subset--------------------------------------------------------------
gs[1]

## ----plotTree------------------------------------------------------------
plot(gs)

## ----getNodes-path-1-----------------------------------------------------
getNodes(gs, path = 1)

## ----getNodes-path-full--------------------------------------------------
getNodes(gs, path = "full")

## ----getNodes-path-auto--------------------------------------------------
nodelist <- getNodes(gs, path = "auto")
nodelist

## ----getGate-------------------------------------------------------------
node <- nodelist[3]
g <- getGate(gs, node)
g

## ----getStats------------------------------------------------------------
getPopStats(gs)[1:10,]

## ----plotGate-nodeName---------------------------------------------------
plotGate(gs, "pDC")

## ----annotate------------------------------------------------------------
d <- data.frame(sample=factor(c("sample 1", "sample 2")),treatment=factor(c("sample","control")) )
pd <- pData(gs)
pd <- cbind(pd,d)
pData(gs) <- pd
pData(gs)

## ------------------------------------------------------------------------
subset(gs, treatment == "control")

## ------------------------------------------------------------------------
fs <- getData(gs)
class(fs)
nrow(fs[[1]])

## ----getData-gh----------------------------------------------------------
fs <- getData(gs, node)
nrow(fs[[1]])

## ----gh------------------------------------------------------------------
gh <- gs[[1]]
gh

## ------------------------------------------------------------------------
head(getPopStats(gh))

## ----plotPopCV-----------------------------------------------------------
plotPopCV(gh)

## ------------------------------------------------------------------------
plotGate(gh)

## ----getInd--------------------------------------------------------------
table(getIndices(gh,node))

## ----getCMAT-------------------------------------------------------------
C <- getCompensationMatrices(gh);
C

## ----getTrans,results='markup'-------------------------------------------
T <- getTransformations(gh)
names(T)
T[[1]]

## ----create gs-----------------------------------------------------------
data(GvHD)
#select raw flow data
fs <- GvHD[1:2]

## ----GatingSet constructor-----------------------------------------------
gs <- GatingSet(fs)

## ----compensate----------------------------------------------------------
cfile <- system.file("extdata","compdata","compmatrix", package="flowCore")
comp.mat <- read.table(cfile, header=TRUE, skip=2, check.names = FALSE)
## create a compensation object 
comp <- compensation(comp.mat)
#compensate GatingSet
gs <- compensate(gs, comp)

## ----eval=FALSE----------------------------------------------------------
#  gs <- compensate(gs, comp.list)

## ----user-transformation-------------------------------------------------
require(scales)
trans.func <- asinh
inv.func <- sinh
trans.obj <- trans_new("myAsinh", trans.func, inv.func)

## ----transform-build-in--------------------------------------------------
trans.obj <- asinhtGml2_trans()
trans.obj

## ----transformerList-----------------------------------------------------
chnls <- colnames(fs)[3:6] 
transList <- transformerList(chnls, trans.obj)

## ----estimateLogicle-----------------------------------------------------
estimateLogicle(gs[[1]], chnls)

## ----transform-gs--------------------------------------------------------
gs <- transform(gs, transList)
getNodes(gs) 

## ----add-rectGate--------------------------------------------------------
rg <- rectangleGate("FSC-H"=c(200,400), "SSC-H"=c(250, 400), filterId="rectangle")
nodeID <- add(gs, rg)
nodeID
getNodes(gs)  

## ----add-quadGate--------------------------------------------------------
qg <- quadGate("FL1-H"= 0.2, "FL2-H"= 0.4)
nodeIDs <- add(gs,qg,parent="rectangle")
nodeIDs 
getNodes(gs)

## ----add-boolGate--------------------------------------------------------
bg <- booleanFilter(`CD15 FITC-CD45 PE+|CD15 FITC+CD45 PE-`)
bg
nodeID2 <- add(gs,bg,parent="rectangle")
nodeID2
getNodes(gs)

## ----plot-gh,eval=FALSE--------------------------------------------------
#  plot(gs, bool=TRUE)

## ----recompute-----------------------------------------------------------
recompute(gs)

## ----plotGate-rect-------------------------------------------------------
plotGate(gs,"rectangle") #plot one Gate

## ----plotGate-multiple---------------------------------------------------
plotGate(gs,getChildren(gs[[1]], "rectangle")) 

## ----plotGate-gh-bool,eval=FALSE-----------------------------------------
#  plotGate(gs[[1]], bool=TRUE)

## ----rm------------------------------------------------------------------
Rm('rectangle', gs)
getNodes(gs)

## ----archive,eval=FALSE--------------------------------------------------
#  tmp <- tempdir()
#  save_gs(gs,path = file.path(tmp,"my_gs"))
#  gs <- load_gs(file.path(tmp,"my_gs"))

## ----clone,eval=FALSE----------------------------------------------------
#  gs_cloned <- clone(gs)

