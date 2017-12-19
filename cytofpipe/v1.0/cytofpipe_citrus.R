options(stringsAsFactors = F)
rm(list = ls())

jobid <- as.character(Sys.getenv("JOB_ID"))

input <- paste0(jobid, ".txt")
lines <- readLines(input, n = 5)

dataDirectory <- lines[1]
conditionsFile <- lines[2]
outputDirectory <- lines[3]
markersFile <- lines[4]
asinh_cofactor <- lines[5]

library(citrus)
library(hash)



#——————————————————————--------
#- READ FCS AND EXTRACT MARKERS
#——————————————————————--------

files <- list.files(dataDirectory,pattern='.fcs$', full=TRUE)
usermarkers <- as.character(read.table(markersFile, header = FALSE)[,1])
 
fcs1<-read.FCS(files[1])

allMarkerNames<-as.vector(pData(parameters(fcs1))$name)
allMarkerDesc<-as.vector(pData(parameters(fcs1))$desc)
 
UserName2Desc <-hash()
for(i in 1:length(allMarkerDesc)){
	if(is.na(allMarkerDesc[i])){
 		UserName2Desc[[ allMarkerNames[i] ]] <- allMarkerNames[i]
	}else{
		UserName2Desc[[ allMarkerDesc[i] ]] <- allMarkerDesc[i]
	}
}
if (sum(has.key( usermarkers, UserName2Desc )) == 0) {
  	clear(UserName2Desc)
 	clear(Desc2UserName)
  	for(i in 1:length(allMarkerDesc)){
  		if(!is.na(allMarkerDesc[i])){
 			id <- gsub( "^[^_]+_", "", allMarkerDesc[i])
  			UserName2Desc[[ id ]] <- allMarkerDesc[i]
  			Desc2UserName[[ allMarkerNames[i] ]] <- id
		}
  	}
}
Desc2Name <-hash()
for(i in 1:length(allMarkerDesc)){
 	if(is.na(allMarkerDesc[i])){
 		Desc2Name[[ allMarkerNames[i] ]] <- allMarkerNames[i]
	}else{
		Desc2Name[[ allMarkerDesc[i] ]] <- allMarkerNames[i]
	}
}
  
markersDesc <-vector()
for(i in 1:length(usermarkers)){
	markersDesc[i]<-values(UserName2Desc,  keys=usermarkers[i])
}
markersName <-vector()
for(i in 1:length(markersDesc)){
	markersName[i]<-values(Desc2Name, keys=markersDesc[i])
}


#———————————----———
#- PARAMETERS
#——————————————---

clusteringColumns<-markersName
transformColumns<-allMarkerNames[-c(grep("Time|Event|Cell_length|viability", allMarkerNames, ignore.case = TRUE))]
scaleColumns = transformColumns
transformCofactor <- as.numeric(asinh_cofactor) 
family = "classification"

minimumClusterSizePercent = 0.05
modelTypes = c("pamr","glmnet","sam")

fileSampleSize = 10000
nFolds = 1
featureType = c("abundances")

data_conditions<-read.table(conditionsFile, sep="\t",header=F)
fileList = data.frame(defaultCondition=data_conditions[,1])
labels = as.factor(data_conditions[,2])

modelTypes<-vector()
if(length(levels(labels)) > 2){
	modelTypes = c("pamr","sam")
}else{
	modelTypes = c("pamr","glmnet","sam")
}


#———————————----———
#- CITRUS
#——————————————---

# Read Data
citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory,fileList,fileSampleSize,transformColumns,transformCofactor)

# Cluster all the data
citrus.foldClustering = citrus.clusterAndMapFolds(citrus.combinedFCSSet,clusteringColumns,labels,nFolds)

# Make vector of conditions for analysis. If comparing two conditions, should be 
# two elements - first element is baseline condition and second is comparison condition.
conditions = colnames(fileList)[1]

# Build cluster features
citrus.foldFeatureSet = citrus.calculateFoldFeatureSet(citrus.foldClustering,citrus.combinedFCSSet,
                                                         featureType=featureType,
                                                         minimumClusterSizePercent=minimumClusterSizePercent,
                                                         conditions=conditions
                                                         )

# Build regression models for each model type
citrus.regressionResults = mclapply(modelTypes,citrus.endpointRegress,citrus.foldFeatureSet=citrus.foldFeatureSet,labels=labels,family=family)

# Plot Results for each model
lapply(citrus.regressionResults,plot,outputDirectory=outputDirectory,citrus.foldClustering=citrus.foldClustering,citrus.foldFeatureSet=citrus.foldFeatureSet,citrus.combinedFCSSet=citrus.combinedFCSSet,
	theme="white")


sessionInfo()

