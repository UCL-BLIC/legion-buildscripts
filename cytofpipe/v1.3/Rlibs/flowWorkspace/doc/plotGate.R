## ----eval=FALSE----------------------------------------------------------
#  plotGate(x, y, ...)

## ----echo=FALSE, message=FALSE, results='hide'---------------------------
library(flowWorkspace)
flowDataPath <- system.file("extdata", package = "flowWorkspaceData")
gs <- load_gs(file.path(flowDataPath,"gs_manual"))
# change the CD3 gate as 1d gate
g <- rectangleGate(`<V450-A>` = c(2000, Inf), filterId = "CD3+")
g <- sapply(sampleNames(gs), function(sn)g)
setGate(gs, "CD3+", g)
recompute(gs, "CD3+")

#preset lattice option
old.theme <- flowWorkspace.par.get("theme.novpadding")
new.theme <- lattice:::updateList(old.theme, list(axis.text = list(cex = 1)
                                                  , par.xlab.text = list(cex = 1)
                                                  , par.ylab.text = list(cex = 1)
                                                  )
                                  )

flowWorkspace.par.set(name="theme.novpadding", new.theme)
flowWorkspace.par.set("plotGate", list(arrange = FALSE))
library(knitr)
opts_chunk$set(fig.show = 'hold', fig.width = 3, fig.height = 3, results = 'hide', message = FALSE)


## ----fig.width = 4-------------------------------------------------------
# plot the entire gating tree
plot(gs)
# plot a sub-tree
plot(gs, "CD3+")
plot(gs, "CD4")

## ------------------------------------------------------------------------
# default xbin = 32      
plotGate(gs[[1]], "CD4/38+ DR+", main = "default (xbin = 32)") 
# increase the resolution
plotGate(gs[[1]], "CD4/38+ DR+", xbin = 128, main = "xbin = 128") 
# have the highest resolution by disabling hexbin 
plotGate(gs[[1]], "CD4/38+ DR+", xbin = 0, main = "xbin = 0 (no binning)") 
# speed up plotting by sub-sampling the data points
# most of them it is subsampled data is sufficient to represent the 
# population distribution
plotGate(gs[[1]], "CD4/38+ DR+", sample.ratio = 0.4, main = "sample.ratio = 0.4") 

## ------------------------------------------------------------------------
# hide the strip
plotGate(gs[[1]], "CD4/38+ DR+", strip = FALSE, main = "strip = FALSE")       
# change the gating path length in strip 
plotGate(gs[[1]], "CD4/38+ DR+", path = "full", main = "path = 'full'")
plotGate(gs[[1]], "CD4/38+ DR+", path = 2, main = "path = 2")

## ------------------------------------------------------------------------
# hide the pop stats
plotGate(gs[[1]], "CD4/38+ DR+", stats = FALSE, main = "stats = FALSE")
# adjust the location, size, background of the stats
plotGate(gs[[1]], "CD4/38+ DR+"
                , pos = c(0.6, 0.9)
                , digits = 4
                , par.settings = list(gate.text = list(background = list(fill = "yellow") #stats background
                                                        , cex = 2 #stats font size
                                                        
                                                      )
                                      )
         )

## ------------------------------------------------------------------------
# show the tranformed scale instead of raw scale
plotGate(gs[[1]], "CD4/38+ DR+", raw.scale = FALSE, main = "raw.scale = FALSE")
# use the actual data range instead of instrument range
plotGate(gs[[1]], "CD4/38+ DR+", xlim = "data", ylim = "data", main = "xlim = 'data'")       
# zooming by manually set xlim and ylim
plotGate(gs[[1]], "CD4/38+ DR+", xlim = c(1000, 4000), ylim = c(1000, 4000), main = "xlim = c(1000, 4000)")       
# hide the channel names and only show the marker names in x,y labs
plotGate(gs[[1]], "CD4/38+ DR+", marker.only = TRUE, main = "marker.only = TRUE")      

## ------------------------------------------------------------------------
# change the gate color and width
plotGate(gs[[1]], "CD4/38+ DR+"
                , par.settings = list(gate = list(col = "black"
                                                  , lwd = 3
                                                  , lty = "dotted")
                                      )
         )
# change the panel background and axis label size
plotGate(gs[[1]], "CD4/38+ DR+"
                , par.settings = list(panel.background = list(col = "white")
                                    , par.xlab.text = list(cex = 1.5)
                                    , par.ylab.text = list(cex = 1.5)
                                    , axis.text = list(cex = 1.5)
                                      ) 
         )

## ----fig.width = 5, fig.height= 5----------------------------------------
# display one population as an overlay on top of another population      
plotGate(gs[[1]], "CD4/CCR7- 45RA+", overlay = "CD4/38+ DR+")
#adjust the overlay symbol
plotGate(gs[[1]], "CD4/CCR7- 45RA+"
              , overlay = "CD4/38+ DR+"
              , par.settings = list(overlay.symbol = list(cex = 0.1
                                                          , fill = "black" 
                                                          , bg.alpha = 0.1 #dim the background density
                                                          )
                                    )
         )
#add multiple overlays
plotGate(gs[[1]], "CD4/CCR7- 45RA+"
              , overlay = c("CD4/CCR7+ 45RA+", "CD4/CCR7+ 45RA-")
              , par.settings = list(overlay.symbol = list(cex = 0.05, bg.alpha = 0.1))
         )
# customize the symbol and legend
plotGate(gs[[1]], "CD4/CCR7- 45RA+"
              , overlay = c("CD4/CCR7+ 45RA+", "CD4/CCR7+ 45RA-")
              , par.settings = list(overlay.symbol = list(cex = 0.05, bg.alpha = 0.1))
              , overlay.symbol = list(`CCR7+ 45RA+` = list(fill = "black") #this overwrite the par.settings
                                      , `CD4/CCR7+ 45RA-` = list(fill = "darkgreen")
                                      )
              , key = list(text = list(c("CCR7+ 45RA+", "CCR7+ 45RA-")) # add legend
                          , points = list(col = c("black", "darkgreen"), pch = 19, cex = 0.5)
                          , columns = 2
                           )
         )

## ----echo=FALSE----------------------------------------------------------
flowWorkspace.par.set("plotGate", list(arrange = TRUE))

## ------------------------------------------------------------------------
# plot multiple populations (gates are merged into same panel if they share the same parent and projections)      
plotGate(gs[[1]], c("CD4", "CD8"))      
# flip the projection
plotGate(gs[[1]], c("CD4", "CD8"), projections = list("CD4" = c(x = "CD8", y = "CD4")))      
# dot not merge the gates
plotGate(gs[[1]], c("CD4", "CD8"), merge = FALSE) 
# use gpar to change layout 
# use arrange.main to change the default title
plotGate(gs[[1]], c("CD4", "CD8"), merge = FALSE, gpar = list(nrow = 1), arrange.main = "CD4 vs CD8")      

## ----eval =FALSE---------------------------------------------------------
#  # or return populations as individual trellis objects
#  #and arrange them manually
#  latticeObjs <- plotGate(gs[[1]], c("CD4", "CD8"), merge = FALSE, arrange = FALSE)
#  do.call(grid.arrange, c(latticeObjs, nrow = 1, main = "Tcell subsets"))

## ----echo=FALSE----------------------------------------------------------
# clone 4 samples to prepare lattice plot
gslist <- lapply(1:4, function(i){
              gs_clone <- clone(gs)
              sampleNames(gs_clone) <- paste0(i, ".fcs")
              pData(gs_clone)[["name"]] <- paste0(i, ".fcs")
              gs_clone
              })
gs <- GatingSetList(gslist)
# create study variables
pData(gs)$Stim <- c("Vaccined", "Placebo", "Vaccined", "Placebo")
pData(gs)$PTID <- c("P001", "P001", "P002", "P002")

set.seed(1)
# add noise to CD3 pop of P002
for(i in c(1,2)){
  cd3_ind <- getIndices(gs[[i]], "CD3+")
  nTcell <- length(which(cd3_ind))
  toNoiseInd <- which(cd3_ind)[sort(sample.int(nTcell, 0.8 * nTcell))]
  cd3_sub_sig <- exprs(flowData(gs@data[[i]])[[1]])[toNoiseInd, "<V450-A>"] 
  noise <- cd3_sub_sig - 1000 * pnorm(1:length(toNoiseInd))
  exprs(flowData(gs@data[[i]])[[1]])[toNoiseInd, "<V450-A>"] <- noise
  recompute(gs[[i]], "CD3+")
}

# decrease active Tcells signals from Placebo groups
for(i in c(2,4)){
  act_ind <- getIndices(gs[[i]], "CD4/38+ DR+")
  nActTcell <- length(which(act_ind))
  toNoiseInd <- which(act_ind)[sort(sample.int(nActTcell, 0.5 * nActTcell))]
  sub_sig <- exprs(flowData(gs@data[[i]])[[1]])[toNoiseInd, "<R660-A>"] 
  noise <- sub_sig - 3000 * pnorm(1:length(toNoiseInd))
  exprs(flowData(gs@data[[i]])[[1]])[toNoiseInd, "<R660-A>"] <- noise
  recompute(gs[[i]], "CD4/38+ DR+")
}

## ----fig.width = 5, fig.height= 5----------------------------------------
# condition on `Stim` and `PTID`
plotGate(gs, "CD4/38+ DR+", cond = "Stim+PTID")
# to put the second condition variable to the y axis
latticeExtra::useOuterStrips(plotGate(gs, "CD4/38+ DR+", cond = "Stim+PTID"))

## ----fig.width = 3, fig.height= 3----------------------------------------
plotGate(gs[c(1,3)], "CD3+") 
# change the default y axis
plotGate(gs[c(1,3)], "CD3+", default.y = "<B710-A>") 
# change the plot type
plotGate(gs[c(1,3)], "CD3+", type = "densityplot") 
# fit the 1d gate onto 1d density 
plotGate(gs[c(1,3)], "CD3+", type = "densityplot", fitGate = TRUE) 
# stack them together
plotGate(gs[c(1,3)], "CD3+", type = "densityplot", stack = TRUE) 

## ---- fig.width= 12, fig.height=6----------------------------------------
plotGate(gs[[3]], gpar = list(nrow = 2))

