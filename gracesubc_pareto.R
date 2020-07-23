# setting up packages and conda environment
reticulate::use_condaenv("reticulate_PCHA", conda = "auto",
                         required = TRUE) 
library(Seurat)
library(ggplot2)
library(SingleCellExperiment)
library(scater)
library(scran)
library(BiocManager)
BiocManager::install("vitkl/ParetoTI", dependencies = c("Depends", "Imports", "LinkingTo"))
library(ParetoTI)
library(remotes)
library(doParallel)
library(doRNG)
library(reticulate)




#########################
#### UBC Cells  ####
#########################

ubc_s = readRDS('/Users/graceluettgen/Downloads/BSSP_Pareto_cerebellum-master/ubc_seurat2_3.RDS')
ubc_s = UpdateSeuratObject(ubc_s)
ubc_sce = as.SingleCellExperiment(ubc_s)



# look at mitochondrial-encoded MT genes
mito.genes = grep(pattern = "^mt-",
                  x = rownames(ubc_sce), 
                  value = TRUE)

# fraction of genes that are nuclear-encoded MT genes for each cell (what is the
# difference between this and mitochondria_located_genes?)
ubc_sce$perc.mito = colSums(counts(ubc_sce[mito.genes, ])) / colSums(counts(ubc_sce))
qplot(ubc_sce$perc.mito, geom = "histogram")
hist(ubc_sce$perc.mito)

# look at nuclear-encoded MT genes (find those genes using GO annotations)
# taxonomy id is still 10090 (mouse)
go_annot = map_go_annot(taxonomy_id = 10090, keys = rownames(ubc_sce),
                        columns = c("GOALL"), keytype = "ALIAS",
                        ontology_type = c("CC"))

# search for those genes annotated with GO:0005739 ("mitochondrion")
mitochondria_located_genes = unique(go_annot$annot_dt[GOALL == "GO:0005739", ALIAS])
# fraction of genes that are mitochondrial genes for each cell
ubc_sce$all_mito_genes = colSums(counts(ubc_sce[mitochondria_located_genes, ])) / colSums(counts(ubc_sce))
## Filtering

# remove cells with more less than 1000 or more than 30000 UMI
ubc_sce = ubc_sce[, colSums(counts(ubc_sce)) > 1000 &
                            colSums(counts(ubc_sce)) < 30000]

# remove genes with too many zeros (> 95% cells) (removes Albumine here)
ubc_sce = ubc_sce[rowMeans(counts(ubc_sce) > 0) > 0.05,]
# remove cells with too many zeros (> 85%)
ubc_sce = ubc_sce[,colMeans(counts(ubc_sce) > 0) > 0.15]

# Normalise gene expression by cell sum factors and log-transform
ubc_sce = scran::computeSumFactors(ubc_sce) 
# scaling factor to normalize the counts
# estimate of the total number of UMIs per cell
ubc_sce = scater::logNormCounts(ubc_sce) 

# plus pseudo count of 1, divide by the normalizing factor, and take log
ubc_sce = scater::logNormCounts(ubc_sce, log = FALSE) # just normalize


# Find principal components
ubc_sce = scater::runPCA(ubc_sce, ncomponents = 7,
                             scale = T, exprs_values = "logcounts")


# Plot PCA colored by identity (originally "batch", which does not exist)
scater::plotReducedDim(ubc_sce, ncomponents = 3, dimred = "PCA",
                       colour_by = "ident")

PCs4arch = t(reducedDim(ubc_sce, "PCA"))


# Fit k=2:8 polytopes to UBCs to find which k best describes the data ####
# find archetypes
?k_fit_pch
library(parallel)
arc_ks = k_fit_pch(PCs4arch, ks = 2:8, check_installed = T,
                   bootstrap = T, bootstrap_N = 200, maxiter = 1000,
                   #bootstrap_type = "m", 
                   seed = 2543, 
                   volume_ratio = "t_ratio", # set to "none" if too slow
                   delta=0, conv_crit = 1e-04, order_type = "align",
                   sample_prop = 0.75)

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
plot_arc_var(arc_1, type = "t_ratio", point_size = 2, line_size = 1.5) + theme_bw()

# Examine the polytope with best k & look at known markers of subpopulations ####
# fit a polytope with bootstrapping of cells to see stability of positions
arc = fit_pch_bootstrap(PCs4arch, n = 200, sample_prop = 0.75, seed = 235,
                        noc = 4, delta = 0, conv_crit = 1e-04)

#, type = "m")


# I was able to eliminate the errors that kept popping up by adding this line of
# code--I think that somewhere in the lines above PCs4arch was indirectly modified,
# so I'm just re-making the object here
PCs4arch = t(reducedDim(ubc_sce, "PCA"))

# vertex 1
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Ckb",])),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "UBCs colored by Ckb")

# vertex 2
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Ubb",])),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "UBCs colored by Ubb")


# vertex 3
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Epha6",])),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "UBCs colored by Epha6")

# vertex 4
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Plcb4",])),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "UBCs colored by Plcb4")

#archetype 1
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce[labs$lab[1,2],])),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "UBCs colored by archetype 1")


#Identity plots
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(ubc_sce@colData$ident == "Brinp2_On_UBCs"),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "UBCs colored by on UBCs")
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(ubc_sce@colData$ident == "Calb2_Off_UBCs"),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "UBCs colored by off UBCs")
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(ubc_sce@colData$ident == "Intermediate_UBCs"),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "UBCs colored by intermediate UBCs")


# You can also check which cells have high entropy of logistic regression 
# predictions when classifying all cells in a tissue into cell types. 
# These could have been misclassified by the method and wrongly assigned 
# to UBCs, or these could be doublets.

# find archetypes on all data (allows using archetype weights to describe cells)
arc_1 = fit_pch(PCs4arch, volume_ratio = "t_ratio", maxiter = 500,
                noc = 4, delta = 0,
                conv_crit = 1e-04)
arc_1
# check that positions are similar to bootstrapping average from above
p_pca = plot_arc(arc_data = arc_1, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5, 
                 data_lab = as.numeric(logcounts(ubc_sce["Plcb1",])),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "UBCs colored by Plcb1")

# Find genes and gene sets enriched near vertices ####
# Map GO annotations and measure activities

logcounts(ubc_sce) = as.matrix(logcounts(ubc_sce))
counts(ubc_sce) = as.matrix(counts(ubc_sce))

activ = measure_activity(ubc_sce, # row names are assumed to be gene identifiers
                         which = "BP", return_as_matrix = F,
                         taxonomy_id = 10090, keytype = "ALIAS",
                         lower = 20, upper = 1000,
                         aucell_options = list(aucMaxRank = nrow(ubc_sce) * 0.1,
                                               binary = F, nCores = 1,
                                               plotStats = FALSE))

# Merge distances, gene expression and gene set activity into one matrix
data_attr = merge_arch_dist(arc_data = arc_1, data = PCs4arch, 
                            feature_data = as.matrix(logcounts(ubc_sce)),
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

labs$lab
labs
p_pca = plot_arc(arc_data = arc, data = PCs4arch,
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = activ$ribosomal_large_subunit_biogenesis,
                 text_size = 60, data_size = 6)
plotly::layout(p_pca, title = "ribosomal_large_subunit_biogenesis activity")

# Randomise variables to measure goodness of observed fit ####
# To measure goodness of observed fit I compare observed tetrahedron shape to shape 
# of data with no relationships between variables. This is done by comparing the 
# ratio of tertahedron volume to volume of convex hull, a complex shape that 
# contains all of the data. Empirical p-value is fraction of random t-ratios 
# that are at least as high as the observed t-ratio.

# use permutations within each dimension - this is only possible for less than 8 
# vertices because computing convex hull gets exponentially slower with more dimensions
start = Sys.time()
pch_rand = randomise_fit_pch(PCs4arch, arc_data = arc_1,
                             n_rand = 1000,
                             replace = FALSE, bootstrap_N = NA,
                             volume_ratio = "t_ratio",
                             maxiter = 500, delta = 0, conv_crit = 1e-4,
                             #type = "m", 
                             clust_options = list(cores = 3))
# use type m to run on a single machine or cloud
# type = "m", clust_options = list(cores = 3))
# use clustermq (type cmq) to run as jobs on a computing cluster (higher parallelisation)
# type = "cmq", clust_options = list(njobs = 10)) 

# This analysis took:
Sys.time() - start

# # plot background distribution of t-ratio and show p-value
#plot(pch_rand, type = c("t_ratio"), nudge_y = 5)

pch_rand 
pch_rand$


