inputfiles<- as.character(Sys.getenv("INPUT"))
outputfiles<- as.character(Sys.getenv("OUTPUT"))
markersFile<- as.character(Sys.getenv("MARKERS"))

library("cytofkit") 
library(flowCore)


#-----------------
#- Get input data
#-----------------
files <- list.files(inputfiles,pattern='.fcs$', full=TRUE)
parameters <- as.character(read.table("parameters.txt", header = FALSE)[,1])


#------------------------------------------------------------------
#- Extract (and combine) the expression matrix with transformation
#------------------------------------------------------------------
if(markersFile == "NONE"){
	data_transformed <- cytof_exprsMerge(fcsFiles = files, 
					comp=FALSE, 
					transformMethod = "cytofAsinh", 
					mergeMethod ="fixed", 
					fixedNum=10000)
}else{
	data_transformed <- cytof_exprsMerge(fcsFiles = files, 
					comp=FALSE, 
					markers=parameters,
					transformMethod = "cytofAsinh", 
					mergeMethod ="fixed", 
					fixedNum=10000)
}

#-----------
#- Run tSNE
#-----------
data_transformed_tsne <- cytof_dimReduction(data=data_transformed, method = "tsne")

#-----------------
#- Run Phenograph
#-----------------
cluster_PhenoGraph <- cytof_cluster(xdata = data_transformed, method = "Rphenograph")

#------------------------------------------------
#- Combine transformed data, tSNE and Phenograph
#------------------------------------------------
data_all <- cbind(data_transformed, data_transformed_tsne, 
                     PhenoGraph = cluster_PhenoGraph)
data_all <- as.data.frame(data_all)

#-----------------------------------
#- Phenograph clusters plot on tSNE
#-----------------------------------
png(paste0(outputfiles,"/phenograph_tsne1_tsne2.png"))
cytof_clusterPlot(data=data_all, xlab="tsne_1", ylab="tsne_2", 
                  cluster="PhenoGraph", sampleLabel = FALSE)
dev.off()

#-----------------------------
#- PhenoGraph cluster heatmap
#-----------------------------
PhenoGraph_cluster_mean <- aggregate(. ~ PhenoGraph, data = data_all, mean)
write.table(PhenoGraph_cluster_mean, file=paste0(outputfiles,"/PhenoGraph_cluster_mean.txt"), sep="\t",quote=F, row.names=F);

range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
PhenoGraph_cluster_mean_norm01 <- as.data.frame(lapply(PhenoGraph_cluster_mean, range01))
write.table(PhenoGraph_cluster_mean_norm01, file=paste0(outputfiles,"/PhenoGraph_cluster_mean_norm01.txt"), sep="\t",quote=F, row.names=F);

png(paste0(outputfiles,"/heatmap_norm.png"))
cytof_heatmap(PhenoGraph_cluster_mean_norm01, baseName = "PhenoGraph NormMean")
dev.off()


#------------------------------------
#- save analysis results to FCS file
#------------------------------------
cytof_addToFCS(data_all, rawFCSdir=inputfiles, analyzedFCSdir=paste0(outputfiles,"/analysed_FCS"), 
               transformed_cols = c("tsne_1", "tsne_2"), 
               cluster_cols = c("PhenoGraph"))



