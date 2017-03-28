inputfiles<- as.character(Sys.getenv("INPUT"))
outputdir<- as.character(Sys.getenv("OUTPUT"))
markersFile<- as.character(Sys.getenv("MARKERS"))
configFile<- as.character(Sys.getenv("CONFIG"))

library("cytofkit") 
library(flowCore)
library(ini)
library(hash)


#-----------------
#- Get input data
#-----------------
files <- list.files(inputfiles,pattern='.fcs$', full=TRUE)
config<-read.ini(configFile)
parameters <- as.character(read.table(markersFile, header = FALSE)[,1])

fcs1<-read.FCS(files[1])
markernames<-pData(parameters(fcs1))$name
markerdesc<-pData(parameters(fcs1))$desc
h<-hash()
for(i in 1:length(markerdesc)){
	h[[ strsplit(markerdesc, "[_]")[[i]][2] ]] <- markernames[i]
}
parameters2<-values(h[parameters])

#------------------------------------------------------------------
#- Parse config file
#------------------------------------------------------------------

projectName = "cytofpipe"
transform = "cytofAsinh"
merge ="fixed"
num="10000"
flowsom_num = "15"
cluster<-"phenograph"
visualization<-"tsne"

if(configFile != "${CYTOFPIPE_HOME}/default_config.txt"){
	cluster<-vector()
	visualization<-vector()

	data<-read.ini(configFile)

	transform = data$cytofpipe$TRANSFORM
	num=data$cytofpipe$DOWNSAMPLE
	
	if(length(data$cytofpipe$PHENOGRAPH)==1){tolower(data$cytofpipe$PHENOGRAPH);if(data$cytofpipe$PHENOGRAPH == "yes"){cluster<-c(cluster,"Rphenograph")}}
	if(length(data$cytofpipe$CLUSTERX)==1){tolower(data$cytofpipe$CLUSTERX);if(data$cytofpipe$CLUSTERX == "yes"){cluster<-c(cluster,"ClusterX")}}
	if(length(data$cytofpipe$DENSVM)==1){tolower(data$cytofpipe$DENSVM);if(data$cytofpipe$DENSVM == "yes"){cluster<-c(cluster,"DensVM")}}
	if(length(data$cytofpipe$FLOWSOM)==1){tolower(data$cytofpipe$FLOWSOM);if(data$cytofpipe$FLOWSOM == "yes"){cluster<-c(cluster,"FlowSOM");flowsom_num=data$cytofpipe$FLOWSOM_K}}
	if(length(cluster) == 0){reduction<-c(cluster,"NULL")}
	
	if(length(data$cytofpipe$TSNE)==1){tolower(data$cytofpipe$TSNE);if(data$cytofpipe$TSNE == "yes"){visualization<-c(visualization,"tsne")}}
	if(length(data$cytofpipe$PCA)==1){tolower(data$cytofpipe$PCA);if(data$cytofpipe$PCA == "yes"){visualization<-c(visualization,"pca")}}
	if(length(visualization) == 0){visualization<-c(visualization,"NULL")}

}

#------------------------------------------------------------------
#- Run cytofkit wraper
#------------------------------------------------------------------

analysis_results <- cytofkit(fcsFiles = files, 
                markers = parameters2, 
                projectName = projectName,
                transformMethod = transform, 
                mergeMethod = "fixed",
		fixedNum = as.numeric(num),
                dimReductionMethod = "tsne",
                clusterMethods = cluster,
                visualizationMethods = visualization,
                FlowSOM_k = as.numeric(flowsom_num),
		progressionMethod = "NULL",
                resultDir = outputdir,
                saveResults = TRUE, 
                saveObject = FALSE)


#------------------------------------------------------------------
#- Get scaled and norm01 heatmaps for mean and percentage
#------------------------------------------------------------------

exprs <- as.data.frame(analysis_results$expressionData)
clusterData <- analysis_results$clusterRes
ifMultiFCS <- length(unique(sub("_[0-9]*$", "", row.names(exprs)))) > 1

range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

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

