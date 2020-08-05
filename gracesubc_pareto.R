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

ubc_s = readRDS('/Users/graceluettgen/Desktop/Github/BSSP_Pareto_cerebellum/ubc_seurat2.RDS')
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
## start
ubc_sce = scater::runPCA(ubc_sce, ncomponents = 7,
                             scale = T, exprs_values = "logcounts")


# Plot PCA colored by identity 
scater::plotReducedDim(ubc_sce, ncomponents = 3, dimred = "PCA",
                       colour_by = "ident")


PCs4arch = t(reducedDim(ubc_sce, "PCA"))
ubc_sce

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
##end

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
                        noc = 2, delta = 0, conv_crit = 1e-04)


p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Grm1",])),
                 text_size = 60, data_size = 4) 
plotly::layout(p_pca, title = "UBCs colored by Grm1")
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Plcb4",])),
                 text_size = 60, data_size = 4) 
plotly::layout(p_pca, title = "UBCs colored by Plcb4")

p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Calb2",])),
                 text_size = 60, data_size = 4) 
plotly::layout(p_pca, title = "UBCs colored by Calb2")
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Plcb1",])),
                 text_size = 60, data_size = 4) 
plotly::layout(p_pca, title = "UBCs colored by Plcb1")

p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Epha6",])),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "UBCs colored by Epha6")


# arch 4--top 2 enriched genes (according to labs) and Ephb1
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Kcnip4",])),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "UBCs colored by Kcnip4")
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Dgkb",])),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "UBCs colored by Dgkb")
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Cacna1a",])),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "UBCs colored by Cacna1a")

#Identity plots
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(ubc_sce@colData$ident == "Brinp2_On_UBCs"),
                 text_size = 60, data_size = 2) 
plotly::layout(p_pca, title = "UBCs colored by on UBCs")

p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(ubc_sce@colData$ident == "Calb2_Off_UBCs"),
                 text_size = 60, data_size = 2) 
plotly::layout(p_pca, title = "UBCs colored by off UBCs")

p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(ubc_sce@colData$ident == "Intermediate_UBCs"),
                 text_size = 60, data_size = 4) 
plotly::layout(p_pca, title = "UBCs colored by intermediate UBCs")

# by region

p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(ubc_sce@colData$region == "X"),
                 text_size = 60, data_size = 4) 
plotly::layout(p_pca, title = "UBCs colored by X")

p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(ubc_sce@colData$region == "IX"),
                 text_size = 60, data_size = 4) 
plotly::layout(p_pca, title = "UBCs colored by IX")

p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(ubc_sce@colData$region == "PRM"),
                 text_size = 60, data_size = 4) 
plotly::layout(p_pca, title = "UBCs colored by PRM")

p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Ntng1",])),
                 text_size = 60, data_size = 4) 
plotly::layout(p_pca, title = "UBCs colored by Ntng1")


# You can also check which cells have high entropy of logistic regression 
# predictions when classifying all cells in a tissue into cell types. 
# These could have been misclassified by the method and wrongly assigned 
# to UBCs, or these could be doublets.

# find archetypes on all data (allows using archetype weights to describe cells)

arc_1 = fit_pch(PCs4arch, volume_ratio = "t_ratio",  maxiter = 500,
                noc = 2, delta = 0,
                conv_crit = 1e-04)

# check that positions are similar to bootstrapping average from above
p_pca = plot_arc(arc_data = arc_1, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5, 
                 data_lab = as.numeric(ubc_sce@colData$ident == "Intermediate_UBCs"),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "UBCs colored by intermediate UBCs (Checking above prediction)")

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
labs$enriched_genes[arch_name == "archetype_1",]

p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Kcnip4",])),
                 text_size = 60, data_size = 2) 
plotly::layout(p_pca, title = "UBCs colored by Kcnip4")
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Ntng1",])),
                 text_size = 60, data_size = 2) 
plotly::layout(p_pca, title = "UBCs colored by Ntng1")
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Sgcd",])),
                 text_size = 60, data_size = 2) 
plotly::layout(p_pca, title = "UBCs colored by Sgcd")
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Auts2",])),
                 text_size = 60, data_size = 2) 
plotly::layout(p_pca, title = "UBCs colored by Auts2")




p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Kcnip4",])),
                 text_size = 60, data_size = 2) 
plotly::layout(p_pca, title = "UBCs colored by Kcnip4")
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Dgkb",])),
                 text_size = 60, data_size = 2) 
plotly::layout(p_pca, title = "UBCs colored by Dgkb")
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Ntn1",])),
                 text_size = 60, data_size = 5) 
plotly::layout(p_pca, title = "UBCs colored by Unc5c")
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Dgkg",])),
                 text_size = 60, data_size = 4) 
plotly::layout(p_pca, title = "UBCs colored by Dgkg")
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Arpp21",])),
                 text_size = 60, data_size = 2) 
plotly::layout(p_pca, title = "UBCs colored by Arpp21")





### Glutamate receptors (metabotropic)
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Grm1",])),
                 text_size = 60, data_size = 2) 
plotly::layout(p_pca, title = "UBCs colored by Grm1")
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Grm3",])),
                 text_size = 60, data_size = 2) 
plotly::layout(p_pca, title = "UBCs colored by Grm3")
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Grm8",])),
                 text_size = 60, data_size = 2) 
plotly::layout(p_pca, title = "UBCs colored by Grm8")

### Glutamate receptors (ionotropic)
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Grid1",])),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "UBCs colored by Grid1")
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Grid2",])),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "UBCs colored by Grid2")
p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(ubc_sce["Gria2",])),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "UBCs colored by Gria2")






# plotting by ribosomal_large_subunit_biogenesis
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





###### MODIFIED CEREBELLUM SPATIAL PLOTTING ######

# used for spatial plotting
arc_by_cells = bin_cells_by_arch(data_attr$data, data_attr$arc_col, bin_prop = 0.25, return_names = T)
data_attr$data$Ntng1
data_attr$data$sample_id
arc_by_cells$archetype_1
# the .csv files are located in the google bucket:
# gs://cerebellum_data/vkozarev/coords/
#reticulate::use_condaenv("reticulate_PCHA", conda = "auto",
#    required = TRUE) 
install.packages("sp")
install.packages("viridisLite")
install.packages("RColorBrewer")
# set up plotting coordinates and Polygons
library(sp)
library(viridisLite)
library(Seurat)
library(RColorBrewer)

coordsI<-read.csv("/Users/graceluettgen/Downloads/BSSP_Pareto_cerebellum-master/coords/coordsI_half.csv")
coordsII<-read.csv("/Users/graceluettgen/Downloads/BSSP_Pareto_cerebellum-master/coords/coordsII_half.csv")
coordsIII<-read.csv("/Users/graceluettgen/Downloads/BSSP_Pareto_cerebellum-master/coords/coordsIII_half.csv")
coordsIX<-read.csv("/Users/graceluettgen/Downloads/BSSP_Pareto_cerebellum-master/coords/coordsIX_half.csv")
coordsVI<-read.csv("/Users/graceluettgen/Downloads/BSSP_Pareto_cerebellum-master/coords/coordsVI_half.csv")
coordsVII<-read.csv("/Users/graceluettgen/Downloads/BSSP_Pareto_cerebellum-master/coords/coordsVII_half.csv")
coordsVIII<-read.csv("/Users/graceluettgen/Downloads/BSSP_Pareto_cerebellum-master/coords/coordsVIII_half.csv")
coordsX<-read.csv("/Users/graceluettgen/Downloads/BSSP_Pareto_cerebellum-master/coords/coordsX_half.csv")
coordsCul<-read.csv("/Users/graceluettgen/Downloads/BSSP_Pareto_cerebellum-master/coords/coordscul_half.csv")

coordsan12<-read.csv("/Users/graceluettgen/Downloads/BSSP_Pareto_cerebellum-master/coords/coordsan12.csv")
coordsan22<-read.csv("/Users/graceluettgen/Downloads/BSSP_Pareto_cerebellum-master/coords/coordsan22.csv")
coordsprm2<-read.csv("/Users/graceluettgen/Downloads/BSSP_Pareto_cerebellum-master/coords/coordsprm2.csv")
coordssim2<-read.csv("/Users/graceluettgen/Downloads/BSSP_Pareto_cerebellum-master/coords/coordssim2.csv")
coordscop2<-read.csv("/Users/graceluettgen/Downloads/BSSP_Pareto_cerebellum-master/coords/coordscop2.csv")
coordsf2<-read.csv("/Users/graceluettgen/Downloads/BSSP_Pareto_cerebellum-master/coords/coordsf2.csv")
coordspfl2<-read.csv("/Users/graceluettgen/Downloads/BSSP_Pareto_cerebellum-master/coords/coordspf2.csv")

IP<-Polygon(coordsI,hole=F)
IIP<-Polygon(coordsII,hole=F)
IIIP<-Polygon(coordsIII,hole=F)
CulP<-Polygon(coordsCul,hole=F)
VIP<-Polygon(coordsVI,hole=F)
VIIP<-Polygon(coordsVII,hole=F)
VIIIP<-Polygon(coordsVIII,hole=F)
IXP<-Polygon(coordsIX,hole=F)
XP<-Polygon(coordsX,hole=F)

an12<-Polygon(coordsan12,hole=F)
an22<-Polygon(coordsan22,hole=F)
prm2<-Polygon(coordsprm2,hole=F)
sim2<-Polygon(coordssim2,hole=F)
f2<-Polygon(coordsf2,hole=F)
pfl2<-Polygon(coordspfl2,hole=F)

# created object cop2 (not in orig. file)
cop2 <- Polygon(coordscop2,hole=F)

IP2<-Polygons(list(IP), "I")
IIP2<-Polygons(list(IIP), "II")
IIIP2<-Polygons(list(IIIP), "III")
CulP2<-Polygons(list(CulP), "CUL")
VIP2<-Polygons(list(VIP), "VI")
VIIP2<-Polygons(list(VIIP), "VII")
VIIIP2<-Polygons(list(VIIIP), "VIII")
IXP2<-Polygons(list(IXP), "IX")
XP2<-Polygons(list(XP), "X")
# do only one half of the whole thing
AN1P<-Polygons(list(an12), "AN1")
AN2P<-Polygons(list(an22), "AN2")
PRMP<-Polygons(list(prm2), "PRM")
SIMP<-Polygons(list(sim2), "SIM")
COPP<-Polygons(list(cop2), "COP")
FLOCP<-Polygons(list(f2), "F")
PARP<-Polygons(list(pfl2), "PF")

SPs = SpatialPolygons(list(IP2,
                           IIP2,
                           IIIP2,
                           CulP2,
                           VIP2,
                           VIIP2,
                           VIIIP2,
                           IXP2,
                           XP2,
                           AN1P,
                           AN2P,
                           PRMP,
                           SIMP,
                           COPP,
                           FLOCP,
                           PARP))

# modified code--to generate the plots comparing the archetypes, set cluster
# to one of the four inputs (3 are commented out)

arc_by_cells$archetype_1 = rownames(ubc_s@meta.data)[ubc_s@meta.data$named_liger_ident == "Brinp2_On_UBCs"]
arc_by_cells$archetype_2 = rownames(ubc_s@meta.data)[ubc_s@meta.data$named_liger_ident == "Calb2_Off_UBCs"]
arc_by_cells$archetype_3 = rownames(ubc_s@meta.data)[ubc_s@meta.data$named_liger_ident == "Intermediate_UBCs"]

length(arc_by_cells$archetype_1)
length(arc_by_cells$archetype_2)
length(arc_by_cells$archetype_3)
ubc_s@meta.data$named_liger_ident 
rownames(ubc_s@meta.data)[ubc_s@meta.data$named_liger_ident == "Brinp2_On_UBCs"]
seurat =  readRDS('/Users/graceluettgen/Desktop/Github/BSSP_Pareto_cerebellum/ubc_seurat2.RDS')
cluster = "Ratio of 3 to all" # other archetypes: cluster = "Spatial distribution of archetype 2" 
            # cluster = "Ratio of archetypes 2 : 1" cluster = "Ratio of archetypes 1 : 2"
quantile.p = 0.5
use.pos.expr = T
use.raw = T
order_regions = NULL
do.print = T
color_divergent = T
return_df = F


  # right now this assumes all regions are represented in all cell types -- fix later
  #if (cluster == "Ratio of 2 : 1 divided by the # of cells in the region"|cluster == "Ratio of 1 : 2 divided by the # of cells in the region")
  {
    region_prop = 1
  } #else
  #{
   # region_prop = table(seurat@meta.data$region) / ncol(seurat@raw.data)
  #}

  # watch out for partial matching!
  levels_include = levels(factor(seurat@meta.data$region))
  if (is.null(order_regions)) {
    order_regions = c("I","II", "III", "CUL", "VI", "VII", "VIII",  "IX","X","AN1", "AN2",
                      "PRM", "SIM", "COP", "F",  "PF")
  }
  if (!is.null(cluster)) {
    if (cluster == "Spatial distribution of archetype 1")
      high_express = rownames(seurat@meta.data)[rownames(seurat@meta.data) %in% arc_by_cells$archetype_1]
    else if (cluster == "Spatial distribution of archetype 2")
      high_express = rownames(seurat@meta.data)[rownames(seurat@meta.data) %in% arc_by_cells$archetype_2]
    else if (cluster == "Ratio of 1 to all"|cluster == "Ratio of 2 to all"| cluster == "Ratio of 3 to all")
    {
      high_express1 = rownames(seurat@meta.data)[rownames(seurat@meta.data) %in% arc_by_cells$archetype_1]
      high_express2 = rownames(seurat@meta.data)[rownames(seurat@meta.data) %in% arc_by_cells$archetype_2]
      high_express3 = rownames(seurat@meta.data)[rownames(seurat@meta.data) %in% arc_by_cells$archetype_3]
    }
    else
      high_express = names(seurat@ident)[which(seurat@ident == cluster)]
    
    title = paste0(cluster)
  } else {
    data.use = seurat@raw.data[gene, ]
    if (!use.raw) {
      data.use = seurat@scale.data[gene, ]
      use.scaled = T
    } else {
      use.scaled = F
    }
    if (use.pos.expr) {
      data.use = data.use[data.use > 0]
    }
    cutoff = quantile(data.use, probs = c(quantile.p))
    if (do.print) { print(cutoff) }
    high_express = WhichCells(seurat, subset.name = gene, accept.low = cutoff, 
                              use.raw = use.raw, use.scaled = use.scaled)
    title = paste0(gene, ", quantile = ", quantile.p)
  }
  table(factor(seurat@meta.data[high_express2,]$region, levels = levels_include))
  table(factor(seurat@meta.data[high_express1,]$region, levels = levels_include))
  if (do.print) { print(table(seurat@meta.data[high_express,]$region)) }
  if(cluster == "Ratio of 2 to all")
  {
    #gene_prop = (table(factor(seurat@meta.data[high_express2,]$region, levels = levels_include)) / 
     #table(factor(seurat@meta.data[high_express1,]$region, levels = levels_include))) 
    ratio1 = table(factor(seurat@meta.data[high_express1,]$region, levels = levels_include))/length(high_express1)
    ratio2 = table(factor(seurat@meta.data[high_express2,]$region, levels = levels_include))/length(high_express2)
    ratio3 = table(factor(seurat@meta.data[high_express3,]$region, levels = levels_include))/length(high_express3)
    
    gene_prop = 2*ratio2/(ratio1+ratio2+ratio3)
  } else if (cluster == "Ratio of 1 to all")
  {
    #gene_prop = (table(factor(seurat@meta.data[high_express1,]$region, levels = levels_include)) / 
                        #table(factor(seurat@meta.data[high_express2,]$region, levels = levels_include))) 
    ratio1 = table(factor(seurat@meta.data[high_express1,]$region, levels = levels_include))/length(high_express1)
    ratio2 = table(factor(seurat@meta.data[high_express2,]$region, levels = levels_include))/length(high_express2)
    ratio3 = table(factor(seurat@meta.data[high_express3,]$region, levels = levels_include))/length(high_express3)
    
    gene_prop = 2*ratio1/(ratio1+ratio2+ratio3)
    
  } else if (cluster == "Ratio of 3 to all")
  {
    #gene_prop = (table(factor(seurat@meta.data[high_express1,]$region, levels = levels_include)) / 
    #table(factor(seurat@meta.data[high_express2,]$region, levels = levels_include))) 
    ratio1 = table(factor(seurat@meta.data[high_express1,]$region, levels = levels_include))/length(high_express1)
    ratio2 = table(factor(seurat@meta.data[high_express2,]$region, levels = levels_include))/length(high_express2)
    ratio3 = table(factor(seurat@meta.data[high_express3,]$region, levels = levels_include))/length(high_express3)
    
    gene_prop = 2*ratio3/(ratio1+ratio2+ratio3)
    
  } else if (cluster == "Spatial distribution of archetype 1" | cluster == "Spatial distribution of archetype 2")
  {
    
  } else if (cluster == "Prop2")
  {
    gene_prop = 2.5*(table(factor(seurat@meta.data[high_express2,]$region, levels = levels_include)) / 
                   (table(factor(seurat@meta.data[high_express2,]$region, levels = levels_include))+
                      table(factor(seurat@meta.data[high_express1,]$region, levels = levels_include))))
    
  } else
  {
    gene_prop = table(factor(seurat@meta.data[high_express,]$region, levels = levels_include))/ 
     sum(table(factor(seurat@meta.data[high_express,]$region, levels = levels_include)))
  }
  
  if (do.print) {
    barplot((gene_prop / region_prop)[order_regions], las = 2, main = title, 
            ylab = "Relative Proportion of high expressing cells")
    abline(h = 1.0)
  }
  
  values = data.frame((gene_prop / region_prop)[order_regions])
  row.names(values) = order_regions
  colnames(values) = c('region', 'avg')
  values[['log_avg']] = log(values$avg, base = 2)
  
  
  
  if (return_df) {
    
    return(values)
  } else {
    # look in cerebellum analysis functions to see how these are made 
    Sps.df<-SpatialPolygonsDataFrame(SPs, values, match.ID = TRUE)
    
    # nn<-data.frame(Sps.df@data)
    # nn<-round(nn, 7)
    
    if (color_divergent) {
      
      col_names = RColorBrewer:::brewer.pal(11,"RdBu")
      # cols = colorRampPalette(colors = rev(col_names), space="Lab")(50)
      # palpos = gplots::colorpanel(sum(values$log_avg>0),high = 'Red')
      # palneg = gplots::colorpanel(sum(values$log_avg<0), 'White', 'Blue')
      
      #palpos = colorRampPalette(c('white', 'red'), space = 'Lab')(floor((max(values$avg)-1) * 50))
      #palneg = colorRampPalette(c('white', 'blue'), space = 'Lab')(floor((1-min(values$avg)) * 50))
      
      palneg = colorRampPalette(c('white', 'blue'), space = 'Lab')(1)
      palpos = colorRampPalette(c('white', 'red'), space = 'Lab')(1.5)
      
      #floor((max(values$avg)-1) * 50)
      #max(values$avg)-1
      palette <- c(rev(palneg),palpos)
      
      # final attempt
      cols <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(50))
      max_abs <- max(abs(values$avg - 1))
 
      
      # add an extra 0.1 plus because the top color can get lost
      brk <- lattice::do.breaks(c(0.99-max_abs, 1.01 + max_abs), 50)
      # print(brk)
      first_true <- which.max(brk > min(values$avg))
      last_true <- which.max(brk > max(values$avg))
      # print(first_true)
      # print(last_true)
      # print(brk > max(values$avg))
      # print(length(brk))
      # print(length(cols))
      brk <- brk[(first_true -1):min(last_true, length(brk))]
      cols <- cols[(first_true -1):min(last_true, length(cols))]
      
      spplot(Sps.df, zcol="avg",
             main = title, col.regions= cols,
             at = brk)
   
      # colorkey = list(col = cols, 
      #                 at = brk))
    } else {
      spplot(Sps.df, zcol="avg",
             main = title, col.regions= viridis(50, option = "C"))
    }
  }

  # return ((gene_prop / region_prop)[order_regions])

#}

values
  # If you want to print the color map, you have to take away the method header
  # and set the parameters manually for each case.
  # ubc_s =  readRDS('/Users/graceluettgen/Downloads/BSSP_Pareto_cerebellum-master/ubc_seurat2_3.RDS')
  #spatial_enrichment_map(ubc_s, cluster = 'Archetype 1', color_divergent = T)
  #spatial_enrichment_map(ubc_s, cluster = 'Archetype 2', color_divergent = T)
  #spatial_enrichment_map(ubc_s, cluster = 'Archetype 3', color_divergent = T)
  #spatial_enrichment_map(ubc_s, cluster = 'Archetype 4', color_divergent = T)
  

  # Spatial prevalence of archetypes

