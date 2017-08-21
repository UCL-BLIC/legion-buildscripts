options(stringsAsFactors = F)
rm(list = ls())

jobid <- as.character(Sys.getenv("JOB_ID"))

input <- paste0(jobid, ".txt")
lines <- readLines(input, n = 5)

working.dir <- lines[1]
ref.file <- paste0(lines[2], ".clustered.txt")
outputdir <- lines[3]
markersFile <- lines[4]
asinh_cofactor <- lines[5]


library(scaffold)
library(flowCore)
library(tools)
library(igraph)
library(reshape)
library(ggrepel)



#————————————
#- CLUSTERING
#————————————

col.names <- as.character(read.table(markersFile, header = FALSE)[,1])
num.cores <- as.numeric(1)
num_clusters <- as.numeric(200)
num_samples <- as.numeric(50)
asinh.cofactor <- as.numeric(asinh_cofactor) 

scaffold:::cluster_fcs_files_in_dir(working.dir, num.cores, col.names, num_clusters, num_samples, asinh.cofactor)


#——————————————
#- RUN SCAFFOLD
#——————————————


ew_influence<- ceiling(length(col.names) / 3)

files.list <- list.files(path = working.dir, pattern = "*.clustered.txt$")
files.list <- files.list[files.list != ref.file]
print(sprintf("Markers used for SCAFFoLD: %s", paste(col.names, collapse = ", ")))

files.list <- c(ref.file, files.list)
print(paste("Using as reference", files.list[1], sep = " "))

ref.dir <- paste(working.dir, "gated/", sep = "/")
gated_data <- scaffold:::load_attractors_from_gated_data(ref.dir, asinh.cofactor)
tab.attractors <- gated_data$tab.attractors
att.labels <- gated_data$cellType_key$population
G.attractors <- NULL

ret <- list(graphs = list(), clustered.data = list())

for(fi in files.list)
{
    f <- paste0(working.dir,"/",fi)
    print (paste("Processing", f, sep = " "))
    tab <- read.table(f, header = T, sep = "\t", quote = "", check.names = F, comment.char = "", stringsAsFactors = F)
    col.names.inter_cluster <- col.names
    tab <- tab[!apply(tab[, col.names], 1, function(x) {all(x == 0)}),]    
    names(tab) <- gsub("cellType", "groups", names(tab))
    names(tab) <- gsub("^X", "", names(tab))
    print(sprintf("Running with Edge weight: %f", ew_influence))
    res <- scaffold:::process_data(tab, G.attractors, tab.attractors,
         col.names = col.names, att.labels = att.labels, already.clustered = T, ew_influence = ew_influence, 
         col.names.inter_cluster = col.names.inter_cluster, inter.cluster.connections = T, overlap_method = NULL)
    G.complete <- scaffold:::get_highest_scoring_edges(res$G.complete)
    clustered.data <- tab 
    ret$graphs[basename(f)] <- list(G.complete)
    ret$clustered.data[basename(f)] <- list(clustered.data)
        
    G.attractors <- res$G.attractors
}

ret <- c(ret, list(dataset.statistics = scaffold:::get_dataset_statistics(ret)))
ret <- c(list(scaffold.col.names = col.names, landmarks.data = gated_data$downsampled.data), ret)
ref=gsub(pattern = "\\.fcs.*", "", ref.file)
scaffold:::my_save(ret, paste(outputdir, sprintf("%s.scaffold", ref), sep = "/"))



#———————————----———
#- SCAFFOLD PLOTS
#——————————————---

f_name <- paste(outputdir, sprintf("%s.scaffold", ref), sep = "/")

con <- file(f_name, "rb")
data <- unserialize(con)
close(con)

#- Get min and max x and y coordinates
x<-vector()
y<-vector()
for (i in 1:length(data$graphs)) {
	G <- data$graphs[[i]]
	
#	layout<-layout.auto(G)
#	layout<-layout.forceatlas2(G, iterations=1000, plotstep=500)
#	fixed <- rep(FALSE, vcount(G))
#	fixed[1:length(att.labels)] <- TRUE
#	layout<-scaffold:::layout.forceatlas2(G, fixed = fixed)
	
	layout<-cbind(V(G)$x,V(G)$y)  
	x<-c(x, layout[,1])
	y<-c(y, layout[,2])
}
xlim=c(min(x), max(x))
ylim=c(min(y), max(y))

for (i in 1:length(data$graphs)) {
	G <- data$graphs[[i]]
	
	name=names(data$graphs[i])
	name2=gsub(pattern = "\\.fcs.*", "", name)
	
	labels<-c(V(G)$name[1:length(att.labels)], rep(NA, length(V(G)$name)-length(att.labels)))
	colores<-c(rep(rgb(1,0,0,0.9), length(att.labels)), rep(rgb(0,0,1,.1), length(V(G)$name)-length(att.labels)))
	sizes<-c(rep(10, length(att.labels)), rep(5, length(V(G)$name)-length(att.labels)))
#	V(G)$label.cex<-c(rep(1, length(att.labels)), rep(0.8, length(V(G)$name)-length(att.labels)))
    

	#-get the node coordinates
#	plotcord <- data.frame(layout.auto(G))
	plotcord <- data.frame(cbind(V(G)$x,V(G)$y), row.names=V(G)$name)
	colnames(plotcord) = c("X1","X2")
	
	#-get edges, which are pairs of node IDs
	edgelist <- get.edgelist(G)
	
	#-convert to a four column edge data frame with source and destination coordinates
	edges <- data.frame(plotcord[edgelist[,1],], plotcord[edgelist[,2],])
	colnames(edges) <- c("X1","Y1","X2","Y2")
	
	pdf(paste0(outputdir,"/scaffold_map_",name2,".pdf"))
	p<-ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = 0.5, colour="grey") + 
		geom_point(aes(X1, X2), size =sizes, colour=colores, data=plotcord) +
		geom_text_repel(aes(X1, X2),data=plotcord, label = labels) +
		xlim(xlim) +
		ylim(ylim) +
		theme(axis.line=element_blank(),
			axis.text.x=element_blank(),
			axis.text.y=element_blank(),
			axis.ticks=element_blank(),
			axis.title.x=element_blank(),
			axis.title.y=element_blank(),
			legend.position="none",
			panel.background=element_blank(),
			panel.border=element_blank(),
			panel.grid.major=element_blank(),
			panel.grid.minor=element_blank(),
			plot.background=element_blank())

	print(p)
	dev.off()

}


