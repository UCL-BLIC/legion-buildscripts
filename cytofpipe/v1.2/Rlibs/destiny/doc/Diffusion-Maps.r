
library(IRkernel)
library(IRdisplay)
library(repr)
library(base64enc)

suppressPackageStartupMessages({
    library(xlsx)
    library(destiny)
    library(Biobase)
})

options(device = function(...) png('/dev/null', 7, 6, 'in', res = 120))
options(repr.plot.width = 7, repr.plot.height = 6)
options(jupyter.plot_mimetypes = c('application/pdf', 'image/png'))

setHook('on.rgl.close', function(...) {
    name <- tempfile()
    par3d(windowRect = c(0, 0, 1200, 1200))
    Sys.sleep(1)
    
    rgl.snapshot(  filename = paste0(name, '.png'))
   #rgl.postscript(filename = paste0(name, '.pdf'), fmt='pdf')  # doesn’t work with spheres
    
    res <- getOption('repr.plot.res')
    
    publish_mimebundle(list(
        'image/png'       = base64encode(paste0(name, '.png'))
     #, 'application/pdf' = base64encode(paste0(name, '.pdf'))
    ), list(
        width  = res * getOption('repr.plot.width'),
        height = res * getOption('repr.plot.height')
    ))
}, 'replace')

library(xlsx)
raw_ct <- read.xlsx('mmc4.xls', sheetName = 'Sheet1')

raw_ct[1:9, 1:9]  #preview of a few rows and columns

library(destiny)
library(Biobase)

ct <- as.ExpressionSet(raw_ct)
ct

num_cells <- gsub('^(\\d+)C.*$', '\\1', ct$Cell)
ct$num_cells <- as.integer(num_cells)

# cells from 2+ cell embryos
have_duplications <- ct$num_cells > 1
# cells with values ≤ 28
normal_vals <- apply(exprs(ct), 2, function(smp) all(smp <= 28))

cleaned_ct <- ct[, have_duplications & normal_vals]

housekeepers <- c('Actb', 'Gapdh')  # houskeeper gene names

normalizations <- colMeans(exprs(cleaned_ct)[housekeepers, ])

guo_norm <- cleaned_ct
exprs(guo_norm) <- exprs(guo_norm) - normalizations

library(destiny)
#data(guo_norm)
dm <- DiffusionMap(guo_norm)

plot(dm)

palette(cube_helix(6)) # configure color palette

plot(dm, pch = 20,         # pch for prettier points
     col_by = 'num_cells', # or “col” with a vector or one color
     legend_main = 'Cell stage')

plot(dm, 1:2, pch = 20, col_by = 'num_cells',
     legend_main = 'Cell stage')

library(rgl)
plot3d(eigenvectors(dm)[, 1:3],
       col = log2(guo_norm$num_cells),
       type = 's', radius = .01)
view3d(theta = 10, phi = 30, zoom = .8)
# now use your mouse to rotate the plot in the window
rgl.close()

library(ggplot2)
qplot(DC1, DC2, data = dm, colour = factor(num_cells)) +
    scale_color_cube_helix()
# or alternatively:
#ggplot(dif, aes(DC1, DC2, colour = factor(num.cells))) + ...

plot(eigenvalues(dm), ylim = 0:1, pch = 20,
     xlab = 'Diffusion component (DC)', ylab = 'Eigenvalue')

oh <- options('repr.plot.height')
options(repr.plot.height = 3)

par(mfrow = c(1, 2), mar = c(2,2,2,2), pch = 20)

plot(dm, 3:4,   col_by = 'num_cells', draw_legend = FALSE)
plot(dm, 19:20, col_by = 'num_cells', draw_legend = FALSE)

options(oh)

hist(exprs(cleaned_ct)['Aqp3', ], breaks = 20,
     xlab = 'Ct of Aqp3', main = 'Histogram of Aqp3 Ct',
     col = palette()[[4]], border = 'white')

dilutions <- read.xlsx('mmc6.xls', 1L)
dilutions$Cell <- NULL  # remove annotation column

get_lod <- function(gene) gene[[max(which(gene != 28))]]

lods <- apply(dilutions, 2, get_lod)
lod <- ceiling(median(lods))
lod

lod_norm <- ceiling(median(lods) - mean(normalizations))
max_cycles_norm <- ceiling(40 - mean(normalizations))

list(lod_norm = lod_norm, max_cycles_norm = max_cycles_norm)

guo <- guo_norm
exprs(guo)[exprs(cleaned_ct) >= 28] <- lod_norm

thresh_dm <- DiffusionMap(guo,
                          censor_val = lod_norm,
                          censor_range = c(lod_norm,
                                           max_cycles_norm),
                          verbose = FALSE)

plot(thresh_dm, 1:2, col_by = 'num_cells', pch = 20,
     legend_main = 'Cell stage')

# remove rows with divisionless cells
ct_w_missing <- ct[, ct$num_cells > 1L]
# and replace values larger than the baseline
exprs(ct_w_missing)[exprs(ct_w_missing) > 28] <- NA

housekeep <- colMeans(exprs(ct_w_missing)[housekeepers, ],
                      na.rm = TRUE)

w_missing <- ct_w_missing
exprs(w_missing) <- exprs(w_missing) - housekeep

exprs(w_missing)[is.na(exprs(ct_w_missing))] <- lod_norm

dif_w_missing <- DiffusionMap(w_missing,
                              censor_val = lod_norm,
                              censor_range = c(lod_norm,
                                               max_cycles_norm),
                              missing_range = c(1, 40),
                              verbose = FALSE)

plot(dif_w_missing, 1:2, col_by = 'num_cells', pch = 20,
     legend_main = 'Cell stage')

ct64 <- guo[, guo$num_cells == 64]

dm64 <- DiffusionMap(ct64)

ct32 <- guo[, guo$num_cells == 32]
pred32 <- dm_predict(dm64, ct32)

par(mar = c(2,2,1,5), pch = 20)
plot(dm64,    1:2,     col     = palette()[[6]],
     new_dcs = pred32, col_new = palette()[[4]])
colorlegend(c(32L, 64L), palette()[c(4,6)])
