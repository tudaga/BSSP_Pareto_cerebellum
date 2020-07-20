# Install ParetoTI package, this should also install reticulate package, if not - install manually.
install.packages("BiocManager") # for installing BioConductor dependencies
install.packages("remotes")
BiocManager::install("vitkl/ParetoTI", dependencies = c("Depends", "Imports", "LinkingTo"))

# Load package
library(ParetoTI)
# If package does not load because "py_pcha is not found" make sure you do not have
# a python environment already loaded in this R session (e.g. restart R and try loading again).

# Install python dependencies (like py_pcha) into python conda environment,
# and (optionally) install *extra_packages*.
ParetoTI::install_py_pcha(method = "conda", conda = "/Users/tudaga/anaconda3/bin/conda",
                          extra_packages = c("tensorflow", "tensorflow-probability",
                                             "pandas", "keras", "h5py",
                                             "geosketch", "pydot", "scikit-learn==0.20",
                                             "umap-learn"))
reticulate::py_discover_config("py_pcha")

reticulate::use_condaenv("reticulate_PCHA", conda = "/Users/tudaga/anaconda3/bin/conda",
                         required = TRUE) # set TRUE to force R to use reticulate_PCHA

library(ParetoTI)
library(ggplot2)

# Generate random data that fits into the triangle (3D)
set.seed(4355)
archetypes = generate_arc(arc_coord = list(c(5, 0, 4), c(-10, 15, 0), c(-30, -20, -5)),
                          mean = 0, sd = 1)
data = generate_data(archetypes$XC, N_examples = 1e4, jiiter = 0.04, size = 0.9)

install_py_pcha()

# Fit polytope to those data
arc_data = fit_pch(data, noc = as.integer(3), delta = 0)

# Show results as interactive 3D scatterplot using plotly
plot_arc(arc_data = arc_data, data = data,
         which_dimensions = 1:3)
# Plot static 2D scatterplot using ggplot2
plot_arc(arc_data = arc_data, data = data,
         which_dimensions = 1:2) +
  theme_bw()

install.packages('umap')
# Plot data as 2D density rather than scatterplot
plot_arc(arc_data = arc_data, data = data,
         which_dimensions = 1:2, geom = ggplot2::geom_bin2d) +
  theme_bw()

# Project to UMAP coordinates (3D -> 2D)
umap = reticulate::import('umap')
arc_umap = arch_to_umap(arc_data, data, which_dimensions = 1:2,
                        method ="umap") # umap is the only option

plot_arc(arc_data = arc_umap$arc_data, data = arc_umap$data,
         which_dimensions = 1:2) + theme_bw()

# Project to tSNE coordinates (3D -> 2D, requires Rtsne package)
arc_tsne = arch_to_tsne(arc_data, data, which_dimensions = 1:2)
plot_arc(arc_data = arc_tsne$arc_data, data = arc_tsne$data,
         which_dimensions = 1:2) + theme_bw()

#########################
#### Liver Tutorial  ####
#########################
install.packages("GEOquery")
#Warning in install.packages 
#  package GEOquery is not available (for R version 4.0.2)
R.version

# library(GEOquery)


# uncomment to load data -------------------------------------------------------
# gse = GEOquery::getGEO("GSE84498", GSEMatrix = TRUE)
# filePaths = GEOquery::getGEOSuppFiles("GSE84498", fetch_files = T, baseDir = "./processed_data/")

# filePaths = c("./processed_data/GSE84498/GSE84498_experimental_design.txt.gz",
#              "./processed_data/GSE84498/GSE84498_umitab.txt.gz")
# design = fread(filePaths[1], stringsAsFactors = F)
# data = fread(filePaths[2], stringsAsFactors = F, header = T)

# data = as.matrix(data, rownames = "gene")

# this is just my working directory -- you can skip this line
ubc_s = readRDS('../ubc_seurat2_3.RDS')
# these are some of the important metadata used
# including sex, cerebellar region (lobule), and named cluster
head(ubc_s@meta.data)

# convert to single cell experiment
hepatocytes = as.SingleCellExperiment(ubc_s)

# look at mitochondrial-encoded MT genes
mito.genes = grep(pattern = "^mt-",
                  x = rownames(data), 
                  value = TRUE)
hepatocytes$perc.mito = colSums(counts(hepatocytes[mito.genes, ])) / colSums(counts(hepatocytes))
#qplot(hepatocytes$perc.mito, geom = "histogram")

# look at nuclear-encoded MT genes (find those genes using GO annotations)
go_annot = map_go_annot(taxonomy_id = 10090, keys = rownames(hepatocytes),
                        columns = c("GOALL"), keytype = "ALIAS",
                        ontology_type = c("CC"))










