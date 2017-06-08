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

require(cytofkit)
require(ggplot2)
require(reshape2)
require(plyr)
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
parameters <- as.character(read.table(markersFile, header = FALSE)[,1])


## @knitr fcs1

fcs1<-read.FCS(files[1])
markernames<-pData(parameters(fcs1))$name
markerdesc<-pData(parameters(fcs1))$desc
h<-hash()
for(i in 1:length(markerdesc)){
        h[[ markerdesc[i] ]] <- markernames[i]
}
if (sum(has.key( parameters, h )) == 0) {
        clear(h)
        for(i in 1:length(markerdesc)){
                id <- gsub( "^[^_]+_", "", markerdesc[i])
                h[[ id ]] <- markernames[i]
        }
}

parameters2<-values(h[parameters])


#------------------------------------------------------------------
#- Parse config file
#------------------------------------------------------------------

## @knitr parseConfig

projectName = "cytofpipe"

cluster<-vector()
visualization<-vector()
visualization<-c(visualization,"tsne")

config<-read.ini(configFile)

autogating=config$cytofpipe$GATING
transform = config$cytofpipe$TRANSFORM
merge = config$cytofpipe$MERGE
num = 10000
flowsom_num = 15
	
if(length(config$cytofpipe$MERGE)==1){tolower(config$cytofpipe$MERGE);if(config$cytofpipe$MERGE == "fixed" || config$cytofpipe$MERGE == "ceil"){num=config$cytofpipe$DOWNSAMPLE}}

if(length(config$cytofpipe$PHENOGRAPH)==1){tolower(config$cytofpipe$PHENOGRAPH);if(config$cytofpipe$PHENOGRAPH == "yes"){cluster<-c(cluster,"Rphenograph")}}
if(length(config$cytofpipe$CLUSTERX)==1){tolower(config$cytofpipe$CLUSTERX);if(config$cytofpipe$CLUSTERX == "yes"){cluster<-c(cluster,"ClusterX")}}
if(length(config$cytofpipe$DENSVM)==1){tolower(config$cytofpipe$DENSVM);if(config$cytofpipe$DENSVM == "yes"){cluster<-c(cluster,"DensVM")}}
if(length(config$cytofpipe$FLOWSOM)==1){tolower(config$cytofpipe$FLOWSOM);if(config$cytofpipe$FLOWSOM == "yes"){cluster<-c(cluster,"FlowSOM");flowsom_num=config$cytofpipe$FLOWSOM_K}}
if(length(cluster) == 0){cluster<-c(cluster,"NULL")}
	
if(length(config$cytofpipe$PCA)==1){tolower(config$cytofpipe$PCA);if(config$cytofpipe$PCA == "yes"){visualization<-c(visualization,"pca")}}
if(length(config$cytofpipe$ISOMAP)==1){tolower(config$cytofpipe$ISOMAP);if(config$cytofpipe$ISOMAP == "yes"){visualization<-c(visualization,"isomap")}}


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

analysis_results <- cytofkit(fcsFiles = files, 
		markers = parameters2, 
		projectName = projectName,
		transformMethod = transform, 
		mergeMethod = merge,
		fixedNum = as.numeric(num),
		dimReductionMethod = "tsne",
		clusterMethods = cluster,
		visualizationMethods = visualization,
		FlowSOM_k = as.numeric(flowsom_num),
		progressionMethod = "NULL",
		resultDir = outputdir,
		saveResults = TRUE, 
		saveObject = TRUE)




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
for(i in 1:length(visualization)){
	vis<-visualization[i]
	
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


