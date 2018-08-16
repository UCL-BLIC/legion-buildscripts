
options(repr.plot.width = 7, repr.plot.height = 6)
options(jupyter.plot_mimetypes = c('application/pdf', 'image/png'))

library(destiny)
data(guo_norm)

sigmas <- find_sigmas(guo_norm, verbose = FALSE)
optimal_sigma(sigmas)

par(pch = 20, mfrow = c(2, 2), mar = c(3,2,2,2))
palette(cube_helix(6))

for (sigma in list('local', 5, round(optimal_sigma(sigmas), 2), 100))
    plot(DiffusionMap(guo_norm, sigma), 1:2,
         main = substitute(sigma == s, list(s = sigma)),
         col_by = 'num_cells', draw_legend = FALSE)
