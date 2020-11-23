# setting up packages and conda environment
reticulate::use_condaenv("reticulate_PCHA", conda = "auto",
                         required = TRUE) 

library(ggplot2)
library(scater)
library(scran)
library(SingleCellExperiment)
library(parallel)
library(Seurat)
library(ParetoTI)
library(BiocManager)
library(remotes)
library(doParallel)
library(doRNG)
library(reticulate)



#########################
#### MLI Cells  ####
#########################


utils_path = "/Users/graceluettgen/Desktop/Github/BSSP_Pareto_cerebellum/ParetoTI_utils.R"
rds_path = '/Users/graceluettgen/Desktop/Github/BSSP_Pareto_cerebellum/ints_seurat.RDS'
save_path = '/Users/graceluettgen/Desktop/Github/BSSP_Pareto_cerebellum/gracesMLI1_analysis'


mli_s = readRDS(rds_path)
mli_s = UpdateSeuratObject(mli_s)
allNeurons = as.SingleCellExperiment(mli_s)
allMLIs = allNeurons[,allNeurons$ident %in% c("Grm8_MLI1", "Npas3_MLI1", "Nxph1_MLI2")]
mli1s = allMLIs[,allMLIs$ident %in% c("Grm8_MLI1", "Npas3_MLI1")]
mli2 = allMLIs[,allMLIs$ident %in% c("Nxph1_MLI2")]


  # look at mitochondrial-encoded MT genes for mli1s
  mito.genes_mli1s = grep(pattern = "^mt-",
                      x = rownames(mli1s), 
                      value = TRUE)
  
  # fraction of genes that are nuclear-encoded MT genes for each cell
  mli1s$perc.mito = colSums(counts(mli1s[mito.genes_mli1s, ])) / colSums(counts(mli1s))
  qplot(mli1s$perc.mito, geom = "histogram")
  hist(mli1s$perc.mito)
  
  # look at nuclear-encoded MT genes (find those genes using GO annotations)
  # taxonomy id is still 10090 (mouse)
  go_annot_mli1s = map_go_annot(taxonomy_id = 10090, keys = rownames(mli1s),
                          columns = c("GOALL"), keytype = "ALIAS",
                          ontology_type = c("CC"))
  
  # search for those genes annotated with GO:0005739 ("mitochondrion")
  mitochondria_located_genes_mli1s = unique(go_annot_mli1s$annot_dt[GOALL == "GO:0005739", ALIAS])
  # fraction of genes that are mitochondrial genes for each cell
  mli1s$all_mito_genes = colSums(counts(mli1s[mitochondria_located_genes_mli1s, ])) / colSums(counts(mli1s))
 
   ## Filtering
  # remove cells with more less than 1000 or more than 30000 UMI
  mli1s = mli1s[, colSums(counts(mli1s)) > 1000 &
                      colSums(counts(mli1s)) < 30000]
  
  # remove genes with too many zeros (> 95% cells) 
  mli1s = mli1s[rowMeans(counts(mli1s) > 0) > 0.05,]
  # remove cells with too many zeros (> 85%)
  mli1s = mli1s[,colMeans(counts(mli1s) > 0) > 0.15]
  
  # Normalise gene expression by cell sum factors and log-transform
  mli1s = scran::computeSumFactors(mli1s) 
  # scaling factor to normalize the counts
  # estimate of the total number of UMIs per cell
  mli1s = scater::logNormCounts(mli1s) 
  
  
  # Find principal components
  mli1s = scater::runPCA(mli1s, ncomponents = 7,
                           scale = T, exprs_values = "logcounts")
  
  # just checking--should ncomponents be 3 instead of 2? we are coloring by ident...
  # but also technically we are treating mli1s as 3 components
  # Plot PCA colored by identity 
  scater::plotReducedDim(mli1s, ncomponents = 2, dimred = "PCA",
                         colour_by = "ident")
  

  PCs4arch_mli1s = t(reducedDim(mli1s, "PCA"))
  mli1s
  # Fit k=2:8 polytopes to MLI1s to find which k best describes the data ####
  # find archetypes
  ?k_fit_pch
  library(parallel)
  arc_ks_mli1s = k_fit_pch(PCs4arch_mli1s, ks = 2:8, check_installed = T,
                     bootstrap = T, bootstrap_N = 200, maxiter = 1000,
                     #bootstrap_type = "m", 
                     seed = 2543, 
                     volume_ratio = "t_ratio", # set to "none" if too slow
                     delta=0, conv_crit = 1e-04, order_type = "align",
                     sample_prop = 0.75)


# Examine the polytope with best k & look at known markers of subpopulations ####
# fit a polytope with bootstrapping of cells to see stability of positions
arc_mli1s = fit_pch_bootstrap(PCs4arch_mli1s, n = 200, sample_prop = 0.75, seed = 235,
                        noc = 3, delta = 0, conv_crit = 1e-04)


p_pca = plot_arc(arc_data = arc_mli1s, data = PCs4arch_mli1s, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(mli1s["Grm8",])),
                 text_size = 60, data_size = 2) 
plotly::layout(p_pca, title = "MLIs colored by Grm8")


# find archetypes on all data (allows using archetype weights to describe cells)
arc_1_mli1s = fit_pch(PCs4arch_mli1s, volume_ratio = "t_ratio",  maxiter = 500,
                noc = 3, delta = 0,
                conv_crit = 1e-04)

# check that positions are similar to bootstrapping average from above
p_pca = plot_arc(arc_data = arc_mli1s, data = PCs4arch_mli1s, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(mli1s["Grm8",])),
                 text_size = 60, data_size = 2) 
plotly::layout(p_pca, title = "MLIs colored by Grm8")

# Find genes and gene sets enriched near vertices ####
# Map GO annotations and measure activities

logcounts(mli1s) = as.matrix(logcounts(mli1s))
counts(mli1s) = as.matrix(counts(mli1s))

activ_mli1s = measure_activity(mli1s, # row names are assumed to be gene identifiers
                         which = "BP", return_as_matrix = F,
                         taxonomy_id = 10090, keytype = "ALIAS",
                         lower = 20, upper = 1000,
                         aucell_options = list(aucMaxRank = nrow(mli1s) * 0.1,
                                               binary = F, nCores = 1,
                                               plotStats = FALSE))


# Merge distances, gene expression and gene set activity into one matrix
data_attr_mli1s = merge_arch_dist(arc_data = arc_1_mli1s, data = PCs4arch_mli1s, 
                            feature_data = as.matrix(logcounts(mli1s)),
                            colData = activ_mli1s,
                            dist_metric = c("euclidean", "arch_weights")[1],
                            colData_id = "cells", rank = F) 


# Use Wilcox test to find genes maximally expressed in 10% closest to each vertex
enriched_genes_mli1s = find_decreasing_wilcox(data_attr_mli1s$data, data_attr_mli1s$arc_col,
                                        features = data_attr_mli1s$features_col,
                                        bin_prop = 0.1, method = "BioQC")

enriched_sets_mli1s = find_decreasing_wilcox(data_attr_mli1s$data, data_attr_mli1s$arc_col,
                                       features = data_attr_mli1s$colData_col,
                                       bin_prop = 0.1, method = "BioQC")

# Take a look at top genes and functions for each archetype
labs_mli1s = get_top_decreasing(summary_genes = enriched_genes_mli1s, 
                                summary_sets = enriched_sets_mli1s,
                          cutoff_genes = 0.01, cutoff_sets = 0.05, 
                          cutoff_metric = "wilcoxon_p_val", 
                          p.adjust.method = "fdr",
                          order_by = "mean_diff", order_decreasing = T,
                          min_max_diff_cutoff_g = 0.4, min_max_diff_cutoff_f = 0.03)


# plotting by ribosomal_large_subunit_biogenesis
p_pca = plot_arc(arc_data = arc_mli1s, data = PCs4arch_mli1s,
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = activ_mli1s$ribosomal_large_subunit_biogenesis,
                 text_size = 60, data_size = 6)
plotly::layout(p_pca, title = "ribosomal_large_subunit_biogenesis activity")

noc = 3
save(list = ls(all.names = TRUE),
     file=paste0(save_path, '/gracesMLI1_pareto',
                 noc,'.Rdata'), envir = .GlobalEnv)

