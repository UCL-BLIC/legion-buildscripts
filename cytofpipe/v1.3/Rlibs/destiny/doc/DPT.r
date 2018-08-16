
options(repr.plot.width = 7, repr.plot.height = 6)
options(jupyter.plot_mimetypes = c('application/pdf', 'image/png'))
set.seed(1)

library(destiny)  # load destiny…
data(guo)         # …and sample data

old <- options(repr.plot.height = 2)
par(mar = rep(0, 4))
graph <- igraph::graph_from_literal(
    data -+ 'transition probabilities' -+ DiffusionMap,
    'transition probabilities' -+ DPT)
plot(
    graph, layout = igraph::layout_as_tree,
    vertex.size = 50,
    vertex.color = 'transparent',
    vertex.frame.color = 'transparent',
    vertex.label.color = 'black')
options(old)

dm <- DiffusionMap(guo)
dpt <- DPT(dm)

plot(dpt, pch = 20)  # “pch” for prettier points

old <- options(repr.plot.height = 3)

par(mfrow = c(1,2), pch = 20, mar = c(2,2,0,1))
plot(dpt, col_by = 'DPT3')
plot(dpt, col_by = 'Gata4', pal = viridis::magma)

options(old)

plot(dpt, root = 2, paths_to = c(1,3), col_by = 'branch', pch = 20)

plot(dpt, col_by = 'branch', divide = 3, dcs = c(-1,3,-2), pch = 20)
