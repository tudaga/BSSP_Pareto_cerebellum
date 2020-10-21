# setting up packages and conda environment
reticulate::use_condaenv("reticulate_PCHA", conda = "auto", required = TRUE) 
library(ggplot2)
library(scater)
library(scran)
library(SingleCellExperiment)
library(parallel)
library(Seurat)
library(ParetoTI)

# PLIs Purkinje Layer Interneurons (inhibitory)
# PLI1 Globular
# PLI2 Candellabrum
# PLI3 Lugaro

# paths that need to be adapted to local paths
utils_path = "/Users/tudaga/Documents/github/BSSP_Pareto_cerebellum/ParetoTI_utils.R"
rds_path = '/Users/tudaga/Documents/github/BSSP_Pareto_cerebellum/data/plis.RDS'
save_path = '/Users/tudaga/Documents/github/BSSP_Pareto_cerebellum/Lug_Glob_Cand_analysis'

data_s = readRDS(rds_path)
data_s = UpdateSeuratObject(data_s)
hvgs = VariableFeatures(data_s)
data_sce = as.SingleCellExperiment(data_s)
rm(data_s)
data_sce

# look at nuclear-encoded MT genes (find those genes using GO annotations)
# taxonomy id is still 10090 (mouse)
go_annot = map_go_annot(taxonomy_id = 10090, keys = rownames(data_sce),
                        columns = c("GOALL"), keytype = "ALIAS",
                        ontology_type = c("CC"))

## Filtering

# only keep hvgs
length(hvgs)
hvgs_bool = rownames(data_sce) %in% hvgs
sum(hvgs_bool)
data_sce = data_sce[hvgs_bool, ]
data_sce


# Normalize gene expression by cell sum factors and log-transform
# scaling factor to normalize the counts
# estimate of the total number of UMIs per cell
data_sce = scran::computeSumFactors(data_sce) 
# plus pseudo count of 1, divide by the normalizing factor, and take log
data_sce = scater::logNormCounts(data_sce) 

# Find principal components
data_sce = scater::runPCA(data_sce, ncomponents = 7, scale = T, 
                          exprs_values = "logcounts")

# Plot PCA colored by identity (cell subtype)
scater::plotReducedDim(data_sce, ncomponents = 3, dimred = "PCA",
                       colour_by = "ident")

# Plot PCA colored by sex
scater::plotReducedDim(data_sce, ncomponents = 3, dimred = "PCA",
                       colour_by = "orig.ident")

# Plot PCA colored by spatial region
scater::plotReducedDim(data_sce, ncomponents = 3, dimred = "PCA",
                       colour_by = "region")


PCs4arch = t(reducedDim(data_sce, "PCA"))
# Fit k=2:8 polytopes to cells to find which k best describes the data ####
# find archetypes
arc_ks = k_fit_pch(PCs4arch, ks = 2:5, check_installed = T,
                   bootstrap = T, bootstrap_N = 200, maxiter = 1000,
                   bootstrap_type = "m",
                   seed = 2543,
                   volume_ratio = "t_ratio", # set to "none" if too slow
                   delta=0, conv_crit = 1e-04, order_type = "align",
                   sample_prop = 0.75, 
                   clust_options = c(cluster_type = "FORK"))

# Show variance explained by a polytope with each k (cumulative)
plot_arc_var(arc_ks, type = "varexpl", point_size = 2, line_size = 1.5) + theme_bw()

# Show variance explained by k-vertex model on top of k-1 model (each k separately)
plot_arc_var(arc_ks, type = "res_varexpl", point_size = 2, line_size = 1.5) + theme_bw()

# Show variance in position of vertices obtained using bootstraping 
# - use this to find largest k that has low variance
plot_arc_var(arc_ks, type = "total_var", point_size = 2, line_size = 1.5) +
  theme_bw() +
  ylab("Mean variance in position of vertices")

# Show t-ratio
plot_arc_var(arc_ks, type = "t_ratio", point_size = 2, line_size = 1.5) + theme_bw()

# Examine the polytope with best k & look at known markers of subpopulations ####
# fit a polytope with bootstrapping of cells to see stability of positions
noc = 4 # number of vertices has to be chosen for each analysis
arc4 = fit_pch_bootstrap(PCs4arch, n = 200, sample_prop = 0.75, seed = 235,
                         noc = noc, delta = 0, conv_crit = 1e-04)

p_pca = plot_arc(arc_data = arc4, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 2.5,
                 data_lab = data_sce@colData$ident,
                 text_size = 60, data_size = 3) 
plotly::layout(p_pca, title = paste0("Colored by subtype ident"))
#for (ident in unique(data_sce@colData$ident)) {
#  print(make_3d_pc_plot_ident(data_sce, arc4, ident, 2.5, 2))
#}

source(utils_path)
gene_vec = c("Aldh1a3", "Slc6a5", "Htr2a", "Nxph1", "Cdh22") # "Edil",
print(gene_vec)
for (gene in gene_vec) {
  print(make_3d_pc_plot_gene(data_sce, arc4, gene, 2.5, 3))
}

for (region in unique(data_sce@colData$region)) {
  print(make_3d_pc_plot_region(data_sce, arc4, region, 2.5, 2))
}

p_pca = plot_arc(arc_data = arc4, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 2.5,
                 data_lab = data_sce@colData$orig.ident,
                 text_size = 60, data_size = 2) 
plotly::layout(p_pca, title = paste0("Colored by sex"))
#for (sex in unique(data_sce@colData$orig.ident)) {
#  print(make_3d_pc_plot_sex(data_sce, arc4, sex, 2.5, 2))
#}


# find archetypes on all data (allows using archetype weights to describe cells)
arc_all = fit_pch(PCs4arch, volume_ratio = "t_ratio", maxiter = 500,
                  noc = noc, delta = 0,
                  conv_crit = 1e-04)
# check that positions are similar to bootstrapping average from above
p_pca = plot_arc(arc_data = arc_all, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 2.5,
                 data_lab = data_sce@colData$ident,
                 text_size = 60, data_size = 2) 
plotly::layout(p_pca, title = paste0("Colored by subtype ident"))

print(gene_vec)
for (gene in gene_vec) {
  print(make_3d_pc_plot_gene(data_sce, arc_all, gene, 2.5, 2))
}

for (region in unique(data_sce@colData$region)) {
  print(make_3d_pc_plot_region(data_sce, arc_all, region, 2.5, 2))
}

p_pca = plot_arc(arc_data = arc_all, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 2.5,
                 data_lab = data_sce@colData$orig.ident,
                 text_size = 60, data_size = 2) 
plotly::layout(p_pca, title = paste0("Colored by sex"))


# Find genes and gene sets enriched near vertices ####
# Map GO annotations and measure activities
#library(doRNG)
#library(doParallel)
logcounts(data_sce) = as.matrix(logcounts(data_sce))
counts(data_sce) = as.matrix(counts(data_sce))

activ = measure_activity(data_sce, # row names are assumed to be gene identifiers
                         which = "BP", return_as_matrix = F,
                         taxonomy_id = 10090, keytype = "ALIAS",
                         lower = 20, upper = 1000,
                         aucell_options = list(aucMaxRank = nrow(data_sce) * 0.1,
                                               binary = F, nCores = 1,
                                               plotStats = FALSE))

# Merge distances, gene expression and gene set activity into one matrix
data_attr = merge_arch_dist(arc_data = arc_all, data = PCs4arch, 
                            feature_data = as.matrix(logcounts(data_sce)),
                            colData = activ,
                            dist_metric = c("euclidean", "arch_weights")[1],
                            colData_id = "cells", rank = F) 


# Use Wilcox test to find genes maximally expressed in 10% closest to each vertex
enriched_genes = find_decreasing_wilcox(data_attr$data, data_attr$arc_col,
                                        features = data_attr$features_col,
                                        bin_prop = 0.1, method = "BioQC")
enriched_genes

enriched_sets = find_decreasing_wilcox(data_attr$data, data_attr$arc_col,
                                       features = data_attr$colData_col,
                                       bin_prop = 0.1, method = "BioQC")

enriched_sets

# Take a look at top genes and functions for each archetype
labs = get_top_decreasing(summary_genes = enriched_genes, summary_sets = enriched_sets,
                          cutoff_genes = 0.01, cutoff_sets = 0.05, 
                          cutoff_metric = "wilcoxon_p_val", 
                          p.adjust.method = "fdr",
                          order_by = "mean_diff", order_decreasing = T,
                          min_max_diff_cutoff_g = 0.4, min_max_diff_cutoff_f = 0.03)

print(labs[["enriched_sets"]][, 1:2])

enriched_genes_arc = list()
i=1
for (archetype in unique(labs[["enriched_genes"]]$arch_name)) {
  select_rows = labs[["enriched_genes"]]$arch_name==archetype
  enriched_genes_arc[[i]] = labs[["enriched_genes"]]$genes[select_rows]
  gene_vec = head(enriched_genes_arc[[i]], 10)
  print(gene_vec)
  for (gene in gene_vec) {
    print(make_3d_pc_plot_gene_arc(data_sce, arc_all, gene, 2.5, 2, i))
  }
  i = i + 1
}


#p_pca = plot_arc(arc_data = arc, data = PCs4arch,
#                 which_dimensions = 1:3, line_size = 2.5,
#                 data_lab = activ$ribosomal_large_subunit_biogenesis,
#                 text_size = 60, data_size = 2)
#plotly::layout(p_pca, title = "ribosomal_large_subunit_biogenesis activity")

# rename the save path below depending on the dataset being analyzed
save(list = ls(all.names = TRUE),
     file=paste0(save_path, '/hvgs_lug_glob_cand',
                 noc,'.Rdata'), envir = .GlobalEnv)
