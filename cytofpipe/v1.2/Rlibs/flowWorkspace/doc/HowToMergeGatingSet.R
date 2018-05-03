## ----eval=FALSE----------------------------------------------------------
#  groupByTree(x)
#  checkRedundantNodes(x)
#  dropRedundantNodes(x,toRemove)
#  dropRedundantChannels(gs, ...)

## ----echo=FALSE, message=FALSE, results='hide'---------------------------
library(flowWorkspace)
flowDataPath <- system.file("extdata", package = "flowWorkspaceData")
gs <- load_gs(file.path(flowDataPath,"gs_manual"))
gs1 <- clone(gs)
sampleNames(gs1) <- "1.fcs"

# simply the tree
nodes <- getNodes(gs1)
for(toRm in nodes[grepl("CCR", nodes)])
  Rm(toRm, gs1)

# remove two terminal nodes
gs2 <- clone(gs1)
sampleNames(gs2) <- "2.fcs"
Rm("DPT", gs2)
Rm("DNT", gs2)

# remove singlets gate
gs3 <- clone(gs2)
Rm("singlets", gs3)
add(gs3, getGate(gs2, "CD3+"), parent = "not debris")
for(tsub in c("CD4", "CD8"))
  {
    add(gs3, getGate(gs2, tsub), parent = "CD3+")
    for(toAdd in getChildren(gs2, tsub))
    {
        thisParent <- getParent(gs2[[1]], toAdd,path="auto")
        add(gs3, getGate(gs2, toAdd), parent = thisParent) 
    }
  }
sampleNames(gs3) <- "3.fcs"

# spin the branch to make it isomorphic
gs4 <- clone(gs3)
# rm cd4 branch first
Rm("CD4", gs4)
# add it back
add(gs4, getGate(gs3, "CD4"), parent = "CD3+")
# add all the chilren back
for(toAdd in getChildren(gs3, "CD4"))
{
    thisParent <- getParent(gs3[[1]], toAdd)
    add(gs4, getGate(gs3, toAdd), parent = thisParent)
}
sampleNames(gs4) <- "4.fcs"

gs5 <- clone(gs4)
# add another redundant node
add(gs5, getGate(gs, "CD4/CCR7+ 45RA+")[[1]], parent = "CD4")
add(gs5, getGate(gs, "CD4/CCR7+ 45RA-")[[1]], parent = "CD4")
sampleNames(gs5) <- "5.fcs"

library(knitr)
opts_chunk$set(fig.show = 'hold', fig.width = 4, fig.height = 4, results= 'asis')


## ----echo=FALSE----------------------------------------------------------
plot(gs1)
plot(gs2)

## ----echo=FALSE----------------------------------------------------------
plot(gs2)
plot(gs3)

## ------------------------------------------------------------------------
invisible(setNode(gs2, "singlets", FALSE))
plot(gs2)
plot(gs3)

## ----results='hold'------------------------------------------------------
getNodes(gs2)[5]
getNodes(gs3)[5]

## ----results='hold'------------------------------------------------------
getNodes(gs2, path = "auto")[5]
getNodes(gs3, path = "auto")[5]

## ----echo=FALSE----------------------------------------------------------
#restore gs2
invisible(setNode(gs2, "singlets", TRUE))

## ----echo=FALSE----------------------------------------------------------
plot(gs3)
plot(gs4)

## ------------------------------------------------------------------------
gslist <- list(gs1, gs2, gs3, gs4, gs5)
gs_groups <- groupByTree(gslist)
length(gs_groups)

## ----error=TRUE----------------------------------------------------------
res <- try(checkRedundantNodes(gs_groups), silent = TRUE)
print(res[[1]])

## ------------------------------------------------------------------------
for(gp in gs_groups)
  plot(gp[[1]])

## ------------------------------------------------------------------------
for(i in c(2,4))
  for(gs in gs_groups[[i]])
    invisible(setNode(gs, "singlets", FALSE))

## ------------------------------------------------------------------------
toRm <- checkRedundantNodes(gs_groups)
toRm

## ----results='hide'------------------------------------------------------
dropRedundantNodes(gs_groups, toRm)

## ------------------------------------------------------------------------
GatingSetList(gslist)

## ------------------------------------------------------------------------
dropRedundantChannels(gs1)

