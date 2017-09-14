## @knitr parameters

jobid <- as.character(Sys.getenv("JOB_ID"))

input <- paste0(jobid, ".txt")
lines <- readLines(input, n = 5)

inputfiles <- lines[1]
outputdir <- lines[2]
markersFile <- lines[3]
configFile <- lines[4]
template <- lines[5]

## @knitr libraries

library(cytofkit) 
library(flowCore)
library(ini)
library(hash)
library(openCyto)
library(mvtnorm)

##require(ggplot2)
##require(plyr)
require(reshape2)
require(VGAM)
require(colourpicker)
require(gplots)


#---------------------------------------------------------------------------------------------
#- Functions
#---------------------------------------------------------------------------------------------

## @knitr functions_opencyto


#- A custom gating function for a DNA/DNA gate on CyTOF data.
#- Finds the intersection between a quantile of a multivariate normal fit
#- of a population and a boundary along y = -x+b 
#- author: jfreling@fhcrc.org
boundry <-  function(xs) {
    # find the boundry events that are above a quantile and below a line 

    cxs <- scale(xs) # scale data so that it can be compaired to the results from qnorm
    f <- qnorm(0.95) # set a boundry level
    pd <- dmvnorm(c(f, f))[1] # and find the p(x) for that level

    pxs <- dmvnorm(x=cxs)  
    idxs <- (pxs > pd) # find those points who are above the boundy level

    idxs2 <- ((-1*cxs[,1]) + 1.96) > cxs[,2] # find points that are below the line with y=-1*x+b 
    pos_xs <- xs[idxs&idxs2,] # intersection of points below line and above threshold level

    hpts <- chull(pos_xs) # find the boundry points of the intersection of cells
    return(pos_xs[hpts,])
}

.dnaGate <- function(fr, pp_res, channels = NA, filterId="", ...){
   xs <- exprs(fr[,channels])
  pnts <- boundry(xs)
  return(polygonGate(.gate=pnts, filterId=filterId))
}

registerPlugins(fun=.dnaGate,methodName='dnaGate', dep='mvtnorm','gating')


## @knitr functions_cytofkit

#- Function to plot all level plots for all markers (https://github.com/JinmiaoChenLab/cytofkit/blob/master/inst/shiny/global.R)
cytof_wrap_colorPlot <- function(data, xlab, ylab, markers, scaleMarker = FALSE,
                             colorPalette = c("bluered", "spectral1", "spectral2", "heat"), 
                             pointSize=1, 
                             removeOutlier = TRUE){
     
     remove_outliers <- function(x, na.rm = TRUE, ...) {
         qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
         H <- 1.5 * IQR(x, na.rm = na.rm)
         y <- x
         y[x < (qnt[1] - H)] <- qnt[1] - H
         y[x > (qnt[2] + H)] <- qnt[2] + H
         y
     }
     
     data <- as.data.frame(data)
     title <- "Marker Expression Level Plot"
     data <- data[,c(xlab, ylab, markers)]
     
     if(removeOutlier){
         for(m in markers){
             data[[m]] <- remove_outliers(data[ ,m])
         }
     }
     
     if(scaleMarker){
         data[ ,markers] <- scale(data[ ,markers], center = TRUE, scale = TRUE)
         ev <- "ScaledExpression"
         data <- melt(data, id.vars = c(xlab, ylab), 
                      measure.vars = markers,
                      variable.name = "markers", 
                      value.name = ev)
     }else{
         ev <- "Expression"
         data <- melt(data, id.vars = c(xlab, ylab), 
                      measure.vars = markers,
                      variable.name = "markers", 
                      value.name = ev)
     }
     
 
     colorPalette <- match.arg(colorPalette)
     switch(colorPalette,
            bluered = {
                myPalette <- colorRampPalette(c("blue", "white", "red"))
            },
            spectral1 = {
                myPalette <- colorRampPalette(c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4",
                                                "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61",
                                                "#F46D43", "#D53E4F", "#9E0142"))
            },
            spectral2 = {
                myPalette <- colorRampPalette(rev(c("#7F0000","red","#FF7F00","yellow","white", 
                                                    "cyan", "#007FFF", "blue","#00007F")))
            },
            heat = {
                myPalette <- colorRampPalette(heat.colors(50))
            }
     )
     zlength <- nrow(data)
     grid_row_num <- round(sqrt(length(markers)))
     gp <- ggplot(data, aes_string(x = xlab, y = ylab, colour = ev)) + 
         facet_wrap(~markers, nrow = grid_row_num, scales = "fixed") +
         scale_colour_gradientn(name = ev, colours = myPalette(zlength)) +
         geom_point(size = pointSize) + theme_bw() + coord_fixed() +
         theme(legend.position = "right") + xlab(xlab) + ylab(ylab) + ggtitle(title) +
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
         theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
     
     return(gp)
}

#- A function to normalize expression values to a 0-1 range
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}


#-----------------
#- Get input data
#-----------------

## @knitr fcs

files <- list.files(inputfiles,pattern='.fcs$', full=TRUE)
files_short <- list.files(inputfiles,pattern='.fcs$', full=F)
usermarkers <- as.character(read.table(markersFile, header = FALSE)[,1])


## @knitr fcs1

fcs1<-read.FCS(files[1])

#— Ideally the user should upload the marker names as provided in the “Description” column of Flowjo. However, if the user uploads the shorter version (i.e., CD38 instead of 141Pr_CD38), the below will do the trick. This is irrelevant for flow data, as the description field given by Flowjo for Flow is always in the short form
#- Also, the colnames in the expression data are on the “Name<Desc>” form, so I am making another hash to be able to substitute these colnames for the markers given by the user (i.e., description or the short version of the description)
 
allMarkerNames<-pData(parameters(fcs1))$name
allMarkerDesc<-pData(parameters(fcs1))$desc
allMarkerNameAndDesc<-paste0(allMarkerNames, "<", allMarkerDesc, ">")
 
UserName2Desc <-hash()
Desc2UserName <-hash()
for(i in 1:length(allMarkerDesc)){
	UserName2Desc[[ allMarkerDesc[i] ]] <- allMarkerDesc[i]
	Desc2UserName[[ allMarkerDesc[i] ]] <- allMarkerDesc[i]
}
if (sum(has.key( usermarkers, UserName2Desc )) == 0) {
 	clear(UserName2Desc)
	clear(Desc2UserName)
 	for(i in 1:length(allMarkerDesc)){
 		id <- gsub( "^[^_]+_", "", allMarkerDesc[i])
 		UserName2Desc[[ id ]] <- allMarkerDesc[i]
 		Desc2UserName[[ allMarkerNames[i] ]] <- id
 	}
}

Desc2NameDesc <-hash()
for(i in 1:length(allMarkerDesc)){
	Desc2NameDesc[[ allMarkerDesc[i] ]] <- allMarkerNameAndDesc[i]
}
 
NameDesc2UserName<-hash()
for(i in 1:length(allMarkerNameAndDesc)){
	UserName=values(Desc2UserName, keys=allMarkerDesc[i])
	if(is.na(UserName)){
	 	NameDesc2UserName[[ allMarkerNameAndDesc[i] ]] <- allMarkerNames[i]
	}else{
	 	NameDesc2UserName[[ allMarkerNameAndDesc[i] ]] <- UserName
	}
}  

markersDesc <-vector()
for(i in 1:length(usermarkers)){
	markersDesc[i]<-values(UserName2Desc,  keys=usermarkers[i])
 }
markersNameDesc <-vector()
for(i in 1:length(markersDesc)){
	markersNameDesc[i]<-values(Desc2NameDesc, keys=markersDesc[i])
 }
markersUserName <-vector()
for(i in 1:length(markersNameDesc)){
	markersUserName[i] <-values(NameDesc2UserName, keys=markersNameDesc[i])
 }
 

#------------------------------------------------------------------
#- Parse config file
#------------------------------------------------------------------

## @knitr parseConfig

projectName = "cytofpipe"

dimReductionMethod="tsne"
clusterMethods<-vector()
visualizationMethods<-vector()
visualizationMethods<-c(visualizationMethods,"tsne")

config<-read.ini(configFile)

autogating=config$cytofpipe$GATING
transformMethod = config$cytofpipe$TRANSFORM
mergeMethod = config$cytofpipe$MERGE
fixedNum = 10000
flowsom_num = 15
perplexity = 30
theta = 0.5
max_iter = 1000
	
if(length(config$cytofpipe$MERGE)==1){tolower(config$cytofpipe$MERGE);if(config$cytofpipe$MERGE == "fixed" || config$cytofpipe$MERGE == "ceil"){fixedNum=config$cytofpipe$DOWNSAMPLE}}
if(length(config$cytofpipe$PERPLEXITY)==1){perplexity=config$cytofpipe$PERPLEXITY}
if(length(config$cytofpipe$THETA)==1){theta=config$cytofpipe$THETA}
if(length(config$cytofpipe$MAX_ITER)==1){max_iter=config$cytofpipe$MAX_ITER}

if(length(config$cytofpipe$PHENOGRAPH)==1){tolower(config$cytofpipe$PHENOGRAPH);if(config$cytofpipe$PHENOGRAPH == "yes"){clusterMethods<-c(clusterMethods,"Rphenograph")}}
if(length(config$cytofpipe$CLUSTERX)==1){tolower(config$cytofpipe$CLUSTERX);if(config$cytofpipe$CLUSTERX == "yes"){clusterMethods<-c(clusterMethods,"ClusterX")}}
if(length(config$cytofpipe$DENSVM)==1){tolower(config$cytofpipe$DENSVM);if(config$cytofpipe$DENSVM == "yes"){clusterMethods<-c(clusterMethods,"DensVM")}}
if(length(config$cytofpipe$FLOWSOM)==1){tolower(config$cytofpipe$FLOWSOM);if(config$cytofpipe$FLOWSOM == "yes"){clusterMethods<-c(clusterMethods,"FlowSOM");flowsom_num=config$cytofpipe$FLOWSOM_K}}
if(length(clusterMethods) == 0){clusterMethods<-c(clusterMethods,"NULL")}

if(length(config$cytofpipe$PCA)==1){tolower(config$cytofpipe$PCA);if(config$cytofpipe$PCA == "yes"){visualizationMethods<-c(visualizationMethods,"pca")}}
if(length(config$cytofpipe$ISOMAP)==1){tolower(config$cytofpipe$ISOMAP);if(config$cytofpipe$ISOMAP == "yes"){visualizationMethods<-c(visualizationMethods,"isomap")}}


#------------------------------------------------------------------
#- Do automatic gating
#------------------------------------------------------------------

## @knitr gating

if(autogating == 'yes'){

	gt<-gatingTemplate(template)

	#------------------------------------------------------------------------------------------------
	#- gates are based on arcSinh transformed data, so raw files need to be transformed before gating
	#------------------------------------------------------------------------------------------------
	
	arcsinh <- arcsinhTransform("arcsinh transformation")	
	dataTransform <- transform(read.flowSet(files), 
			"arcsinh_Ce142Di"= arcsinh(Ce142Di),
			"arcsinh_Ce140Di"= arcsinh(Ce140Di),
			"arcsinh_Ir191Di"= arcsinh(Ir191Di),
			"arcsinh_Ir193Di"= arcsinh(Ir193Di),
			"arcsinh_Y89Di"= arcsinh(Y89Di),
			"arcsinh_Pt195Di"= arcsinh(Pt195Di)
	)

	gs <- GatingSet(dataTransform)
	gating(x = gt, y = gs)
	fs_live <- getData(gs,"Live")

	pdf(paste0(outputdir,"/gating_scheme.pdf"))
	plot(gs)
	dev.off()

	write.flowSet(fs_live, paste0(outputdir, "/gating_fs_live"))

	rm(files)
	rm(files_short)

	files<-list.files(paste0(outputdir, "/gating_fs_live"), patter=".fcs", full=T)
	files_short<-list.files(paste0(outputdir, "/gating_fs_live"), patter=".fcs", full=F)

	for(i in 1:length(files_short)){
		pdf(paste0(outputdir,"/gating_",files_short[i],".pdf"))
		plotGate(gs[[i]])
		dev.off()
	}

}


#------------------------------------------------------------------
#- Run cytofkit wraper
#------------------------------------------------------------------

## @knitr cytofkit

#analysis_results <- cytofkit(fcsFiles = files, 
#		markers = parameters2, 
#		projectName = projectName,
#		transformMethod = transform, 
#		mergeMethod = merge,
#		fixedNum = as.numeric(num),
#		dimReductionMethod = "tsne",
#		clusterMethods = cluster,
#		visualizationMethods = visualization,
#		FlowSOM_k = as.numeric(flowsom_num),
#		progressionMethod = "NULL",
#		resultDir = outputdir,
#		saveResults = TRUE, 
#		saveObject = TRUE)

#- cytof_exprsMerge calls cytof_exprsExtract, which excludes Time and Event channels from the expression matrix, and excludes FSC/SSC from transformation
exprs_data <- cytof_exprsMerge(fcsFiles = files, comp = FALSE, verbose = FALSE, 
                                   transformMethod = transformMethod, 
                                   mergeMethod = mergeMethod, fixedNum = as.numeric(fixedNum))

#- change the colnames here so that the plots show the markers as uploaded by the user
for(i in 1:length(colnames(exprs_data))){
	colnames(exprs_data)[i]<-values(NameDesc2UserName, keys=colnames(exprs_data)[i])

 }

## dimension reduced data, a list
alldimReductionMethods <- unique(c(visualizationMethods, dimReductionMethod))
allDimReducedList <- lapply(alldimReductionMethods, 
                             cytof_dimReduction, data = exprs_data, 
			     markers = markersUserName,
			     perplexity = as.numeric(perplexity),
			     theta = as.numeric(theta),
			     max_iter = as.numeric(max_iter))
names(allDimReducedList) <- alldimReductionMethods

## cluster results, a list
cluster_res <- lapply(clusterMethods, cytof_cluster, 
                          ydata = allDimReducedList[[dimReductionMethod]], 
                          xdata = exprs_data,
                          FlowSOM_k = as.numeric(flowsom_num))
names(cluster_res) <- clusterMethods


## wrap the results
analysis_results <- list(expressionData = exprs_data,
                             dimReductionMethod = dimReductionMethod,
                             visualizationMethods = alldimReductionMethods,
                             dimReducedRes = allDimReducedList,
                             clusterRes = cluster_res, 
                             projectName = projectName,
                             rawFCSdir = inputfiles,
                             resultDir = outputdir,
			     dimRedMarkers = markersUserName)
        
## save the results
cytof_writeResults(analysis_results = analysis_results,
                       saveToRData = TRUE,
                       saveToFCS = TRUE,
                       saveToFiles = TRUE)


#------------------------------------------------------------------
#- Get scaled and norm01 heatmaps for median and percentage
#-	and level Plots
#------------------------------------------------------------------

## @knitr scaledHeatmaps

exprs <- as.data.frame(analysis_results$expressionData)
clusterData <- analysis_results$clusterRes
dimRed<-as.data.frame(analysis_results$dimReducedRes)

ifMultiFCS <- length(unique(sub("_[0-9]*$", "", row.names(exprs)))) > 1

## Level plots
data_all<-cbind(exprs, dimRed)
for(i in 1:length(visualizationMethods)){
	vis<-visualizationMethods[i]
	
	pdf(paste0(outputdir,"/",projectName, "_", vis, "_level_plot.pdf"))
	gp<-cytof_wrap_colorPlot(data=data_all, xlab=paste0(vis, ".",vis,"_1"), ylab=paste0(vis, ".",vis,"_2"), markers=colnames(exprs), colorPalette = c("spectral1"), pointSize=0.1)
	print(gp)
	dev.off()

}

if(!is.null(clusterData) && length(clusterData) > 0){

	## Heatmaps
	for(j in 1:length(clusterData)){
		methodj <- names(clusterData)[j]
		dataj <- clusterData[[j]]
		if(!is.null(dataj)){
                    
			exprs_cluster_sample <- data.frame(exprs, cluster = dataj, check.names = FALSE)
		
			## cluster median 
			cluster_median <- cytof_clusterStat(data= exprs_cluster_sample, cluster = "cluster", statMethod = "median")

			pdf(paste0(outputdir,"/",projectName, "_",methodj, "_cluster_median_heatmap_scaled.pdf"))
			cytof_heatmap(cluster_median, scaleMethod="column", paste(projectName, methodj, "\ncluster median (scaled)", sep = " "))
			dev.off()

			cluster_median_norm01<-as.data.frame( apply(cluster_median, 2, range01))
			pdf(paste0(outputdir,"/",projectName, "_",methodj, "_cluster_median_heatmap_norm01.pdf"))
			cytof_heatmap(cluster_median_norm01, paste(projectName, methodj, "\ncluster median (norm01)", sep = " "))
			dev.off()

			## cluster percentage
			if (ifMultiFCS) {
				cluster_percentage <- cytof_clusterStat(data= exprs_cluster_sample, cluster = "cluster", statMethod = "percentage")
				pdf(paste0(outputdir,"/",projectName, "_",methodj, "_cluster_percentage_heatmap_scaled.pdf"))
				cytof_heatmap(cluster_percentage,scaleMethod="column", paste(projectName, methodj, "cluster\ncell percentage (scaled)", sep = " "))
				dev.off()
			}

		}
	}
}


sessionInfo()
