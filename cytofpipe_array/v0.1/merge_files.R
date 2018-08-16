library(dplyr)
library(stringr)

out <- as.character(Sys.getenv("OUT"))
iter <- as.character(Sys.getenv("ITER"))
clusts <- as.character(Sys.getenv("CLUST"))

cl<-str_split(clusts, ",")

for (i in cl[[1]]){
	for (j in c("median_data", "cell_percentage")){
		filenames <- list.files(paste0(out, "_array_files/",i,"/",j,"/"), pattern="*.csv", full.names=TRUE)
		All <- lapply(filenames,function(i){
			read.csv(i, header=TRUE, stringsAsFactors = FALSE)
		})
		all<-bind_rows(All)
		print(paste0(out, "_array_files/",i,"/",j,"_",iter, "_iters.csv"))
	
		write.table(all, file=paste0(out, "_array_files/",i,"/",j,"_",iter, "_iters.csv"), sep=",", quote=F, row.names=F);
	}
}
