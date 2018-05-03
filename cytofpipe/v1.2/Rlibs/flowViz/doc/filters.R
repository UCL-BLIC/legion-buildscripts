## ----setup, include=FALSE------------------------------------------------
library(knitr)
opts_chunk$set(out.extra='style="display:block; margin: auto"', fig.align="center", message = FALSE, warning = FALSE, knitr.table.format = "markdown")
assign("depthtrigger", 60, data.table:::.global)

## ----loadPackage, echo=FALSE,results='hide'------------------------------
library(flowViz)

## ----loadData------------------------------------------------------------
library(flowViz)
data(GvHD)
head(pData(GvHD))

## ----subsetLattice, eval=FALSE-------------------------------------------
#  xyplot(`SSC-H` ~ `FSC-H` | Visit, data=GvHD,
#         subset=Patient=="6")

## ----subsetFlowCore------------------------------------------------------
GvHD <- GvHD[pData(GvHD)$Patient==6]

## ----transform-----------------------------------------------------------
tf <- transformList(from=colnames(GvHD)[3:7], tfun=asinh)
GvHD <- tf %on% GvHD

## ----simpleRG, fig=TRUE, include=FALSE, echo=TRUE, prefix=FALSE, height=4.5, width=8----
rgate <- rectangleGate("FSC-H"=c(0, 400),"SSC-H"=c(-50, 300))
print(xyplot(`SSC-H` ~ `FSC-H` | Visit, data=GvHD, filter=rgate))

## ----threeDRG------------------------------------------------------------
rgate2 <- rgate * rectangleGate("FL1-H"=c(2, 4))
xyplot(`SSC-H` ~ `FSC-H` | Visit, data=GvHD, filter=rgate2)

## ----nonSmoothScatter----------------------------------------------------
xyplot(`SSC-H` ~ `FSC-H` | Visit, data=GvHD, filter=rgate, smooth=FALSE)

## ----hexbinScatter-------------------------------------------------------
xyplot(`SSC-H` ~ `FSC-H` | Visit, data=GvHD, filter=rgate, smooth=FALSE,xbin=128)

## ----plotnonSmoothScatter, fig=TRUE, include=FALSE, echo=FALSE, results='hide', prefix=FALSE, height=4.5, width=8----
print(plot(trellis.last.object()))

## ----nonSmoothScatterOutline---------------------------------------------
xyplot(`SSC-H` ~ `FSC-H` | Visit, data=GvHD, filter=rgate,
smooth=FALSE, stat=TRUE,pos=0.5,abs=TRUE)

## ----multi-filters-------------------------------------------------------
recGate1<-rectangleGate("FL4-H"=c(3.3,5.3),"FL2-H"=c(2.5,5))
recGate2<-rectangleGate("FL4-H"=c(4,6),"FL2-H"=c(6,8))
filters1<-filters(list(recGate1,recGate2))
xyplot(`FL2-H` ~ `FL4-H`
	   ,data=GvHD[[1]]
       ,filter=filters1
       ,stat=TRUE
       ,margin=FALSE
       )	

## ----plotnonSmoothScatterOutline, fig=TRUE, include=FALSE, echo=FALSE, results='hide', prefix=FALSE, height=4.5, width=8----
print(plot(trellis.last.object()))

## ----nonSmoothScatterOutline2--------------------------------------------
xyplot(`SSC-H` ~ `FSC-H` | Visit, data=GvHD, filter=rgate2,
       smooth=FALSE, outline=TRUE)

## ----plotnonSmoothScatter2, fig=TRUE, include=FALSE, echo=FALSE, results='hide', prefix=FALSE, height=4.5, width=8----
print(plot(trellis.last.object()))

## ----norm2Filter---------------------------------------------------------
n2Filter <- norm2Filter("SSC-H", "FSC-H", scale=2, filterId="Lymphocytes")
n2Filter.results <- filter(GvHD, n2Filter)
xyplot(`SSC-H` ~ `FSC-H` | Visit, data=GvHD, filter=n2Filter.results)

## ----plotnorm2Filter, fig=TRUE, include=FALSE, echo=FALSE, results='hide', prefix=FALSE, height=4.5, width=8----
print(plot(trellis.last.object()))

## ----curv2Filter---------------------------------------------------------
library(flowStats)
c2f <- curv2Filter("FSC-H", "FL4-H", bwFac=1.8)
c2f.results <- filter(GvHD, c2f)
xyplot(`FL4-H` ~ `FSC-H` | Visit, data=GvHD, filter=c2f.results)

## ----plotcurv2Filter, fig=TRUE, include=FALSE, echo=FALSE, results='hide', prefix=FALSE, height=4.5, width=8----
print(plot(trellis.last.object()))

## ----densityplots--------------------------------------------------------
densityplot(~ `FSC-H`, GvHD, filter=curv1Filter("FSC-H"))

## ----plotdensityplots, fig=TRUE, include=FALSE, echo=FALSE, results='hide', prefix=FALSE, height=4.5, width=8----
print(plot(trellis.last.object()))

## ----multDensityplots----------------------------------------------------
densityplot(~ ., GvHD, channels=c("FSC-H", "SSC-H", "FL1-H"), 
filter=list(curv1Filter("FSC-H"), NULL, rgate2))

## ----plotmultDensityplots, fig=TRUE, include=FALSE, echo=FALSE, results='hide', prefix=FALSE, height=4, width=8----
print(plot(trellis.last.object()))

## ----Pars1---------------------------------------------------------------
xyplot(`SSC-H` ~ `FSC-H` | Visit, data=GvHD, filter=rgate,
       par.settings=list(gate=list(fill="black", alpha=0.2)))

## ----plotPars1, fig=TRUE, include=FALSE, echo=FALSE, results='hide', prefix=FALSE, height=4.5, width=8----
print(plot(trellis.last.object()))

## ----plotPars2, eval=FALSE-----------------------------------------------
#  flowViz.par.set(gate=list(fill="black", alpha=0.2))
#  xyplot(`SSC-H` ~ `FSC-H` | Visit, data=GvHD, filter=rgate)

## ----Pars3---------------------------------------------------------------
xyplot(`SSC-H` ~ `FSC-H` | Visit, data=GvHD, filter=rgate,
       smooth=FALSE, 
       par.settings=list(gate=list(col="orange", alpha=0.04,
                                   pch=20, cex=0.7),
                         flow.symbol=list(alpha=0.04, pch=20, cex=0.7)))

## ----plotPars3, fig=TRUE, include=FALSE, echo=FALSE, results='hide', prefix=FALSE, height=4.2, width=8----
print(plot(trellis.last.object()))

## ----Pars4---------------------------------------------------------------
xyplot(`FL4-H` ~ `FSC-H` | Visit, data=GvHD, filter=c2f.results,
par.settings=list(gate=list(fill=rainbow(10), alpha=0.5, col="transparent")))

## ----plotPars4, fig=TRUE, include=FALSE, echo=FALSE, results='hide', prefix=FALSE, height=4.5, width=8----
print(plot(trellis.last.object()))

## ----norm2FilterNames----------------------------------------------------
xyplot(`SSC-H` ~ `FSC-H` | Visit, data=GvHD, filter=n2Filter.results, 
       names=TRUE, 
       par.settings=list(gate=list(fill="black", alpha=0.2, 
                                   col="transparent"),
                         gate.text=list(col="darkred", alpha=0.7, cex=0.6)))

## ----plotnorm2FilterNames, fig=TRUE, include=FALSE, echo=FALSE, results='hide', prefix=FALSE, height=4.5, width=8----
print(plot(trellis.last.object()))

## ----parCoord------------------------------------------------------------
parallel(~ . | Visit, GvHD, filter=n2Filter.results, alpha = 0.01)

## ----plotparCoord, fig=TRUE, include=FALSE, echo=FALSE, results='hide', prefix=FALSE, height=6, width=8----
print(plot(trellis.last.object()))

