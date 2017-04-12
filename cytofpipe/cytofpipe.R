## @knitr parameters

jobid <- as.character(Sys.getenv("JOB_ID"))
#jobid<-"/home/regmond/Scratch/cytof/cytofkit/2818957"

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


#---------------------------------------------------------------------------------------------
#- Functions
#---------------------------------------------------------------------------------------------

## @knitr functions


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
	if(grepl("_",markerdesc[i])){
		h[[ strsplit(markerdesc, "[_]")[[i]][2] ]] <- markernames[i]
	}else{
		h[[ markerdesc[i] ]] <- markernames[i]
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

config<-read.ini(configFile)

autogating=config$cytofpipe$GATING
transform = config$cytofpipe$TRANSFORM
merge = config$cytofpipe$MERGE
	
if(length(config$cytofpipe$MERGE)==1){tolower(config$cytofpipe$MERGE);if(config$cytofpipe$MERGE == "fixed" || config$cytofpipe$MERGE == "ceil"){num=config$cytofpipe$DOWNSAMPLE}}

if(length(config$cytofpipe$PHENOGRAPH)==1){tolower(config$cytofpipe$PHENOGRAPH);if(config$cytofpipe$PHENOGRAPH == "yes"){cluster<-c(cluster,"Rphenograph")}}
if(length(config$cytofpipe$CLUSTERX)==1){tolower(config$cytofpipe$CLUSTERX);if(config$cytofpipe$CLUSTERX == "yes"){cluster<-c(cluster,"ClusterX")}}
if(length(config$cytofpipe$DENSVM)==1){tolower(config$cytofpipe$DENSVM);if(config$cytofpipe$DENSVM == "yes"){cluster<-c(cluster,"DensVM")}}
if(length(config$cytofpipe$FLOWSOM)==1){tolower(config$cytofpipe$FLOWSOM);if(config$cytofpipe$FLOWSOM == "yes"){cluster<-c(cluster,"FlowSOM");flowsom_num=config$cytofpipe$FLOWSOM_K}}
if(length(cluster) == 0){cluster<-c(cluster,"NULL")}
	
if(length(config$cytofpipe$TSNE)==1){tolower(config$cytofpipe$TSNE);if(config$cytofpipe$TSNE == "yes"){visualization<-c(visualization,"tsne")}}
if(length(config$cytofpipe$PCA)==1){tolower(config$cytofpipe$PCA);if(config$cytofpipe$PCA == "yes"){visualization<-c(visualization,"pca")}}
if(length(visualization) == 0){visualization<-c(visualization,"NULL")}


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
#- Get scaled and norm01 heatmaps for mean and percentage
#------------------------------------------------------------------

## @knitr scaledHeatmaps

exprs <- as.data.frame(analysis_results$expressionData)
clusterData <- analysis_results$clusterRes
ifMultiFCS <- length(unique(sub("_[0-9]*$", "", row.names(exprs)))) > 1

if(!is.null(clusterData) && length(clusterData) > 0){
	for(j in 1:length(clusterData)){
		methodj <- names(clusterData)[j]
		dataj <- clusterData[[j]]
		if(!is.null(dataj)){
                    
			exprs_cluster_sample <- data.frame(exprs, cluster = dataj, check.names = FALSE)
		
			## cluster mean 
			cluster_mean <- cytof_clusterStat(data= exprs_cluster_sample, cluster = "cluster", statMethod = "mean")
			pdf(paste0(outputdir,"/",projectName, "_",methodj, "_cluster_mean_heatmap_scaled.pdf"))
			cytof_heatmap(cluster_mean, scaleMethod="column", paste(projectName, methodj, "\ncluster mean (scaled)", sep = " "))
			dev.off()

			cluster_mean_norm01<-as.data.frame( apply(cluster_mean, 2, range01))
			pdf(paste0(outputdir,"/",projectName, "_",methodj, "_cluster_mean_heatmap_norm01.pdf"))
			cytof_heatmap(cluster_mean_norm01, paste(projectName, methodj, "\ncluster mean (norm01)", sep = " "))
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

