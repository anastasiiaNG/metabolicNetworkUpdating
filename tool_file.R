# load data obtained like in ...R (Formal requirments to data structure and properties: ...)
# load("~/Documents/immGen/scell/tms/iniSeuratDataFromCHPC/myelPartAfterFirstFiltOfini250Myel/es_files/es4_woALL_62.rda")
load("~/Documents/immGen/final_objects/243_es.top12k.Rda")
Biobase::pData(es.top12k)$sample <- colnames(es.top12k)

# prepearing working environment:
# work.dir <- "~/Documents/immGen/scell/tms/iniSeuratDataFromCHPC/myelPartAfterFirstFiltOfini250Myel/gam/gam_clustering_es4_91_32_woALL/"
source("csourse.R")
work.dir <- "~/Desktop/"
dir.create(work.dir, showWarnings = F)

### If you don't have entrez IDs in your annotation
entrez <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db,
                                column = "ENTREZID",
                                keytype = "SYMBOL",
                                keys=as.character(fData(es.top12k)$symbol))


### Initial clustering:

#' Initial clustering. Sets the number of initital clusters.
#'
#' @param es value: ExpressionSet
#' @param organism value: choose either `mouse`, either `human`
#' @param repeats e.g. you may use gsub() ti collapce replicas or use any custon annotation
#' @param initialNumber
#' @param showInitialClustering show heatmap of initial clusters
#' @return The table sum of \code{x} and \code{y}.
#'
#' # Currently availble for Mus musculus and GAM method only
preCluster <- preClustering(es = es.top12k,
                            repeats = Biobase::pData(es.top12k)$sample,
                            entrezIDs = Biobase::fData(es.top12k)$entrez,
                            initialNumber = 3,
                            showInitialClustering = T)
str(preCluster)
str(preCluster$network)
str(preCluster$gene.exprs_orig)
str(preCluster$gene.exprs)
str(preCluster$curCenters)

#' GAM-clustering per se
#'
#' @param gene.exprs value: ExpressionSet
#' @param base reflects the degree of retention of edges in the active module
#' @return The table sum of \code{x} and \code{y}.
gamCluster <- gamClustering(network = preCluster$network,
                            gene.exprs = preCluster$gene.exprs,
                            curCenters = preCluster$curCenters,
                            repeats = Biobase::pData(es.top12k)$sample,
                            base = 0.4,
                            work.dir = work.dir,
                            batch.script = "/home/octopus/R-studio/nclust/sgmwcs-slurm-batch-0.9.5",
                            showIntermediateClustering = TRUE,
                            saveSession = TRUE)
gamCluster$k
gamCluster$revs
gamCluster$curRev
gamCluster$curRev$modules
session::restore.session(file=paste0(work.dir, "/session.RDa"))


### Visualizing gam-clustering results:

#' getGraphs
#'
#' @param curRev value:
#' @return results of these functions calling can be seen in work.dir
getGraphs(curRev = gamCluster$curRev,
          work.dir = work.dir)

#' getGraphs
#'
#' @param curRev value:
#' @return results of these functions calling can be seen in work.dir
getHeatmaps(curRev = gamCluster$curRev,
            work.dir = work.dir)

#' getGraphs
#'
#' @param curRev value:
#' @return results of these functions calling can be seen in work.dir
getGeneTables(curRev = gamCluster$curRev,
              work.dir = work.dir,
              organism = "mouse")

### Annotation of gam-clusters:
annotateModules(modules)



