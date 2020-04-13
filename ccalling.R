library(Biobase)
library(cluster) # here pam?
library(limma) # avearrays?
library(stats) # stays::cor
library(pheatmap)
library(rUtils) # normalize.rows
library(igraph)
library(GAM.networks)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(gatom) # load_all("/home/octopus/R-studio/repo/gatom/")
library(dplyr)
library(igraph)
rUtils::messagef()
# or
sprintf()
makeENetworkWithScore() # здесь используется reflink и куча из гатома и GAMа
solver.sb()
getCenter()
processModule() #  # install @neato, take `s` outside?
solver.gmwcs()
solver.gmwcs(net1) # it's also in /home/octopus/R-studio/GAM_files/sgmwcs/sgmwcs-solver.jar
sgmwcs.batchSolver() # path.expand()
getHeatmaps() # reflink!!!
getGeneTables()
GAM::get.edge.attributes()
data.table()
annotateModules() # there we download
gseaReactome # this source function is used here
# make shure you have installed these libs (simpl list file):
Packages <- c("BiocManager", "BioNet", "Biobase", "devtools", "parallel", "logging",
              "cluster", "data.table", "plyr", "pryr", "dplyr", "ggplot2", "cowplot",
              "pheatmap", "igraph", "RColorBrewer", "BiocParallel", "DESeq2", "limma",
              "org.Mm.eg.db", "scales", "gatom", "rUtils", "GAM", "GAM.db", "GAM.networks")
# devtools::install_github("ctlab/gatom")
# devtools::install_github("ctlab/GAM")
# devtools::install_github("assaron/rUtils")
