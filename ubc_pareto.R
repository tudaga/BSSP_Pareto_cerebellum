# Install ParetoTI package, this should also install reticulate package, if not - install manually.
install.packages("BiocManager") # for installing BioConductor dependencies
install.packages("remotes")
BiocManager::install("vitkl/ParetoTI", dependencies = c("Depends", "Imports", "LinkingTo"))
BiocManager::install(c("scran", "scater"))
install.packages("doParallel")
install.packages("doRNG")

# ParetoTI package was supposed to load here
library(ParetoTI)
# If package does not load because "py_pcha is not found" make sure you do not have
# a python environment already loaded in this R session (e.g. restart R and try loading again).

# Install python dependencies (like py_pcha) into python conda environment,
# and (optionally) install *extra_packages*.
ParetoTI::install_py_pcha(method = "conda", 
                          extra_packages = c("tensorflow", "tensorflow-probability",
                                             "pandas", "keras", "h5py",
                                             "geosketch", "pydot", "scikit-learn==0.20",
                                             "umap-learn"))
reticulate::py_discover_config("py_pcha")

# Run separately if you want to change the conda param
# To make sure R uses the correct conda enviroment you can run this when you start R:
reticulate::use_condaenv("reticulate_PCHA", conda = "auto",
                         required = TRUE) # set TRUE to force R to use reticulate_PCHA

# Removed tutorial piece from here

#####################################################################
#### NOT Liver Tutorial, but instead UBC data in liver tutorial  ####
#####################################################################

# No need for GEOquery since we're loading in our own data as a Seurat object
R.version



### Load more (essential) libraries prior to usage

library(Seurat)
library(SingleCellExperiment)
library(scran)
library(scater)
library(ggplot2)
library(ggpubr)
library(data.table)
library(parallel)
library(doRNG)
library(doParallel)

### The real start of this program

###### <- This start indicates an error somewhere near those lines. There is 
# none here though

# Read in UBC data (stored one directory level above)
ubc_s <- readRDS('../ubc_seurat2_3.RDS')

ubc_s <- UpdateSeuratObject(ubc_s)

# Print out meta data
head(ubc_s@meta.data)

# Convert to single cell experiment
ubc_sce = as.SingleCellExperiment(ubc_s)

### Replaced hepatocytes with ubc_sce

# look at mitochondrial-encoded MT genes
mito.genes = grep(pattern = "^mt-",
                  x = rownames(ubc_sce), 
                  value = TRUE)
ubc_sce$perc.mito = colSums(counts(ubc_sce[mito.genes, ])) / colSums(counts(ubc_sce))
qplot(ubc_sce$perc.mito, geom = "histogram")
hist(ubc_sce$perc.mito)

# look at nuclear-encoded MT genes (find those genes using GO annotations)
go_annot = map_go_annot(taxonomy_id = 10090, keys = rownames(ubc_sce),
                        columns = c("GOALL"), keytype = "ALIAS",
                        ontology_type = c("CC"))

mitochondria_located_genes = unique(go_annot$annot_dt[GOALL == "GO:0005739", ALIAS])
ubc_sce$all_mito_genes = colSums(counts(ubc_sce[mitochondria_located_genes, ])) / colSums(counts(ubc_sce))


##########################
### Filtering the data ###
##########################

# subset to the celltype of interest
# remove batches of different cells (probably non-hepatocytes)
ubc_sce = ubc_sce[, !ubc_sce$batch %in% c("AB630", "AB631")]
table(ubc_sce$batch)

# remove cells with more less than 1000 or more than 30000 UMI
ubc_sce = ubc_sce[, colSums(counts(ubc_sce)) > 1000 &
                            colSums(counts(ubc_sce)) < 30000]

# remove cells that express less than 1% of albumine
alb_perc = counts(ubc_sce)["Alb",] / colSums(counts(ubc_sce))

ubc_sce = ubc_sce[, alb_perc > 0.01]

# remove genes with too many zeros (> 95% cells)
ubc_sce = ubc_sce[rowMeans(counts(ubc_sce) > 0) > 0.05,]

# remove cells with too many zeros (> 85%)
ubc_sce = ubc_sce[,colMeans(counts(ubc_sce) > 0) > 0.15]

# Normalise gene expression by cell sum factors and log-transform
ubc_sce = scran::computeSumFactors(ubc_sce) 

# scaling factor to normalize the counts
# estimate of the total number of UMIs per cell
ubc_sce = scater::logNormCounts(ubc_sce) 

?scater::logNormCounts

# plus pseudo count of 1, divide by the normalizing factor, and take log
ubc_sce = scater::logNormCounts(ubc_sce, log = FALSE) # just normalize

##### Error here -- mostly because of sparse matrix

# Find principal components
ubc_sce = scater::runPCA(ubc_sce,
                             scale = T, exprs_values = "logcounts")
# Plot PCA colored by batch
scater::plotReducedDim(ubc_sce, ncomponents = 3, dimred = "PCA",
                       colour_by = "batch")

PCs4arch = t(reducedDim(ubc_sce, "PCA"))

# Fit k=2:8 polytopes to Hepatocytes to find which k best describes the data ####
# find archetypes
?k_fit_pch

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

# Examine the polytope with best k & look at known markers of subpopulations ####
# fit a polytope with bootstrapping of cells to see stability of positions
arc = fit_pch_bootstrap(PCs4arch, n = 200, sample_prop = 0.75, seed = 235,
                        noc = 4, delta = 0, conv_crit = 1e-04)
#, type = "m")

p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Alb",])),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "UBCs colored by Alb (Albumine)")

p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Cyp2e1",])),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "UBCs colored by Cyp2e1")

p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Gpx1",])),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "UBCs colored by Gpx1")

p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Apoa2",])),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "UBCs colored by Apoa2")

# You can also check which cells have high entropy of logistic regression 
# predictions when classifying all cells in a tissue into cell types. 
# These could have been misclassified by the method and wrongly assigned 
# to Hepatocytes, or these could be doublets.

# find archetypes on all data (allows using archetype weights to describe cells)
arc_1 = fit_pch(PCs4arch, volume_ratio = "t_ratio", maxiter = 500,
                noc = 4, delta = 0,
                conv_crit = 1e-04)
# check that positions are similar to bootstrapping average from above
p_pca = plot_arc(arc_data = arc_1, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5, 
                 data_lab = as.numeric(logcounts(ubc_sce["Alb",])),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "UBCs colored by Alb")

# Find genes and gene sets enriched near vertices ####
# Map GO annotations and measure activities
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

enriched_sets = find_decreasing_wilcox(data_attr$data, data_attr$arc_col,
                                       features = data_attr$colData_col,
                                       bin_prop = 0.1, method = "BioQC")

# Take a look at top genes and functions for each archetype
labs = get_top_decreasing(summary_genes = enriched_genes, summary_sets = enriched_sets,
                          cutoff_genes = 0.01, cutoff_sets = 0.05, 
                          cutoff_metric = "wilcoxon_p_val", 
                          p.adjust.method = "fdr",
                          order_by = "mean_diff", order_decreasing = T,
                          min_max_diff_cutoff_g = 0.4, min_max_diff_cutoff_f = 0.03)

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
plot(pch_rand, type = c("t_ratio"), nudge_y = 5)

pch_rand


