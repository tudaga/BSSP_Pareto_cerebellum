make_3d_pc_plot_ident = function(data_sce, arc, ident, line_size, data_size) {
  p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                   which_dimensions = 1:3, line_size = line_size,
                   data_lab = data_sce@colData$ident==ident,
                   text_size = 60, data_size = data_size) 
  return(plotly::layout(p_pca, title = paste0("Colored by ",
                                              ident)))
}

make_3d_pc_plot_gene = function(data_sce, arc, gene, line_size, data_size) {
  p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                   which_dimensions = 1:3, line_size = line_size,
                   data_lab = as.numeric(data_sce@assays@data$logcounts[gene,]),
                   text_size = 60, data_size = data_size) 
  return(plotly::layout(p_pca, title = paste0("Colored by logcounts of ",gene)))
}

make_3d_pc_plot_region = function(data_sce, arc, region, line_size, data_size) {
  p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                   which_dimensions = 1:3, line_size = line_size,
                   data_lab = as.factor(data_sce@colData$region==region),
                   text_size = 60, data_size = data_size) 
  return(plotly::layout(p_pca, title = paste0("Colored by region ",
                                              region)))
}

make_3d_pc_plot_sex = function(data_sce, arc, sex, line_size, data_size) {
  p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                   which_dimensions = 1:3, line_size = line_size,
                   data_lab = data_sce@colData$orig.ident==sex,
                   text_size = 60, data_size = data_size) 
  return(plotly::layout(p_pca, title = paste0("Colored by sex ",
                                              sex)))
}

make_3d_pc_plot_gene_arc = function(data_sce, arc, gene, 
                                    line_size, data_size,
                                    archetype) {
  p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                   which_dimensions = 1:3, line_size = line_size,
                   data_lab = as.numeric(data_sce@assays@data$logcounts[gene, ]),
                   text_size = 60, data_size = data_size) 
  return(plotly::layout(p_pca, title = paste0(archetype,
                                              ", logcounts of ", gene)))
}


make_3d_pc_plot_dist_cl = function(data_sce, arc, cl, color_by, line_size, data_size) {
  p_pca = plot_arc(arc_data = arc, data = PCs4arch[,data_sce@colData$ident==cl], 
                   which_dimensions = 1:3, line_size = line_size,
                   data_lab = log(color_by[data_sce@colData$ident==cl]),
                   text_size = 60, data_size = data_size) 
  return(plotly::layout(p_pca, title = paste0("Log dist to point, cl ", cl)))
}

