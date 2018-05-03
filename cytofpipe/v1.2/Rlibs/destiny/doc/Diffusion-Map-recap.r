
options(repr.plot.width = 7, repr.plot.height = 6)
options(jupyter.plot_mimetypes = c('application/pdf', 'image/png'))

library(destiny)
data(guo)

dm_guo <- DiffusionMap(guo, verbose = FALSE,
                       censor.val = 15, censor.range = c(15, 40))
dm_guo

plot(dm_guo, pch = 20)

palette(cube_helix(6))
plot(dm_guo, col_by = 'num_cells', pch = 20,
     legend_main = 'Cell stage')
