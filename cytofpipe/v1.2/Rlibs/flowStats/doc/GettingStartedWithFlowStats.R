### R code from vignette source 'GettingStartedWithFlowStats.Rnw'

###################################################
### code chunk number 1: loadGvHD
###################################################
library(flowStats)
data(ITN)


###################################################
### code chunk number 2: transform
###################################################
wf <- workFlow(ITN)
tl <- transformList(colnames(ITN)[3:7], asinh, transformationId="asinh")
add(wf, tl)


###################################################
### code chunk number 3: lymphGate
###################################################
lg <- lymphGate(Data(wf[["asinh"]]), channels=c("SSC", "CD3"),
                preselection="CD4", filterId="TCells", eval=FALSE,
                scale=2.5)
add(wf, lg$n2gate, parent="asinh")


###################################################
### code chunk number 4: lymphGatePlot
###################################################
library(flowViz)
print(xyplot(SSC ~ CD3| PatientID, wf[["TCells+"]],
      par.settings=list(gate=list(col="red", 
         fill="red", alpha=0.3))))


###################################################
### code chunk number 5: variation
###################################################
pars <- colnames(Data(wf[["base view"]]))[c(3,4,5,7)]
print(densityplot(PatientID~., Data(wf[["TCells+"]]), channels=pars, groups=GroupID,
                  scales=list(y=list(draw=F)), filter=lapply(pars, curv1Filter),
                  layout=c(4,1)))


###################################################
### code chunk number 6: norm
###################################################
norm <- normalization(normFun=function(x, parameters, ...)
                      warpSet(x, parameters, ...),
                      parameters=pars,
                      arguments=list(grouping="GroupID", monwrd=TRUE),
                      normalizationId="Warping")
add(wf, norm, parent="TCells+")


###################################################
### code chunk number 7: normPlot
###################################################
print(densityplot(PatientID~., Data(wf[["Warping"]]), channels=pars, groups=GroupID,
                  scales=list(y=list(draw=F)), filter=lapply(pars, curv1Filter),
                  layout=c(4,1)))


###################################################
### code chunk number 8: quadGate
###################################################
qgate <- quadrantGate(Data(wf[["Warping"]]), stains=c("CD4", "CD8"),
                      filterId="CD4CD8", sd=3)
add(wf, qgate, parent="Warping")


###################################################
### code chunk number 9: quadGatePlot
###################################################
print(xyplot(CD8 ~ CD4 | PatientID, wf[["CD4+CD8+"]],
             par.settings=list(gate=list(fill="transparent", 
                               col="red"))))


###################################################
### code chunk number 10: rangeGate
###################################################
CD69rg <- rangeGate(Data(wf[["Warping"]]), stain="CD69",
                    alpha=0.75, filterId="CD4+CD8-CD69", sd=2.5)
add(wf, CD69rg, parent="CD4+CD8-")


###################################################
### code chunk number 11: rangeGatePlot
###################################################
print(densityplot(PatientID ~ CD69, Data(wf[["CD4+CD8-"]]), main = "CD4+",
            groups=GroupID, refline=CD69rg@min))


###################################################
### code chunk number 12: createData
###################################################
dat <- Data(wf[["Warping"]])



###################################################
### code chunk number 13: rawData
###################################################
print( xyplot(CD8 ~ CD4 , dat, main= "Experimental data set"))


###################################################
### code chunk number 14: createControlData
###################################################
datComb <- as(dat,"flowFrame")
subCount <- nrow(exprs(datComb))/length(dat)
	sf <- sampleFilter(filterId="mySampleFilter", size=subCount)
	fres <- filter(datComb, sf)
	ctrlData <- Subset(datComb, fres)
	ctrlData <- ctrlData[,-ncol(ctrlData)] ##remove the  column name "original"
	


###################################################
### code chunk number 15: BinControlData
###################################################
minRow=subCount*0.05
refBins<-proBin(ctrlData,minRow,channels=c("CD4","CD8"))
	


###################################################
### code chunk number 16: controlBinsPlot
###################################################
plotBins(refBins,ctrlData,channels=c("CD4","CD8"),title="Control Data")



###################################################
### code chunk number 17: binSampleData
###################################################
sampBins <- fsApply(dat,function(x){
		   binByRef(refBins,x)
		   })


###################################################
### code chunk number 18: pearsonStat
###################################################
pearsonStat <- lapply(sampBins,function(x){
		      calcPearsonChi(refBins,x)
                     })


###################################################
### code chunk number 19: Roderers PBin metric
###################################################
sCount <- fsApply(dat,nrow)
pBStat <-lapply(seq_along(sampBins),function(x){
		calcPBChiSquare(refBins,sampBins[[x]],subCount,sCount[x])
		})


###################################################
### code chunk number 20: plotBinsresiduals
###################################################
par(mfrow=c(4,4),mar=c(1.5,1.5,1.5,1.5))

plotBins(refBins,ctrlData,channels=c("CD4","CD8"),title="Control Data")
patNames <-sampleNames(dat)
tm<-lapply(seq_len(length(dat)),function(x){
		plotBins(refBins,dat[[x]],channels=c("CD4","CD8"),
			title=patNames[x],
			residuals=pearsonStat[[x]]$residuals[2,],
			shadeFactor=0.7)
		
		}
      )


###################################################
### code chunk number 21: chiSqstatisticvalues
###################################################
library(xtable)
chi_Square_Statistic <- unlist(lapply(pearsonStat,function(x){
		x$statistic
		}))

pBin_Statistic <-unlist(lapply(pBStat,function(x){
                x$pbStat
						                        })) 
	
frame <- data.frame(chi_Square_Statistic, pBin_Statistic)
rownames(frame) <- patNames
 	


###################################################
### code chunk number 22: GettingStartedWithFlowStats.Rnw:336-337
###################################################
print(xtable(frame))


###################################################
### code chunk number 23: GettingStartedWithFlowStats.Rnw:343-344
###################################################
toLatex(sessionInfo())


