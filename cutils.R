preClustering <- function(es = es.top12k,
                          organism = "mouse", # human
                          method = "GAM", # gatom
                          repeats = Biobase::pData(es.top12k)$sample,
                          entrezIDs = Biobase::fData(es.top12k)$entrez,
                          initialNumber = 32,
                          showInitialClustering = TRUE){

  if(method == "GAM"){

    if(organism == "mouse") {
      network <- get(data(kegg.mouse.network, package = "GAM.networks"))
      network$graph.raw <- network$graph.raw[-which(network$graph.raw$met.x == "C00288"), ] # HCO3-
      network$graph.raw <- network$graph.raw[-which(network$graph.raw$met.y == "C00288"), ] # kegg.db$mets2mask?

      } else if (organism == "human"){
        network <- get(data(kegg.human.network, package = "GAM.networks"))
        network$graph.raw <- network$graph.raw[-which(network$graph.raw$met.x == "C00288"), ] # HCO3-
        network$graph.raw <- network$graph.raw[-which(network$graph.raw$met.y == "C00288"), ]

        } else {
          stop("Invalid `organism` value: choose either `mouse` or `human`")}

    } else if(method == "gatom"){
      load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.rda"))
      met.to.filter <- data.table::fread(system.file("mets2mask.lst", package="gatom"))$ID
      # met.to.filter <- data.table::fread("/home/octopus/R-studio/repo/gatom/inst/mets2mask.lst")$ID
      # load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/met.kegg.db.rda"))

      if(organism == "mouse") {
        load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rda"))
        # org.Mm.eg.gatom.anno <- gatom::makeOrgGatomAnnotation(org.Mm.eg.db)

        `%>%` <- dplyr::`%>%`
        gatom.net <- "atom"
        if(gatom.net == "atom"){
          globalEdgeTable_pre <- dplyr::tibble(
            rpair = network$reaction2rpair$rpair,
            reaction = network$reaction2rpair$reaction) %>%
            dplyr::inner_join(., network$rpair2align,
                              by = "rpair") %>%
            dplyr::inner_join(., network$enzyme2reaction,
                              by = "reaction") %>%
            dplyr::inner_join(.,  org.Mm.eg.gatom.anno$gene2enzyme,
                              by = "enzyme") %>%
            # rpair | reaction | atom.x | atom.y | enzyme |  gene
            dplyr::rename(from = atom.x, to = atom.y) %>%
            dplyr::select(c(from, to, gene)) %>%
            dplyr::filter(!duplicated(.)) %>%
            dplyr::filter(!gsub("_.*", "", to) %in% met.to.filter) %>%
            dplyr::filter(!gsub("_.*", "", from) %in% met.to.filter) %>%
            dplyr::filter(.$gene %in% entrezIDs)
        } else {
          globalEdgeTable_pre <- dplyr::tibble(
            rpair = network$reaction2rpair$rpair,
            reaction = network$reaction2rpair$reaction) %>%
            dplyr::inner_join(., network$rpairs,
                              by = "rpair") %>%
            dplyr::inner_join(., network$enzyme2reaction,
                              by = "reaction") %>%
            dplyr::inner_join(.,  org.Mm.eg.gatom.anno$gene2enzyme,
                              by = "enzyme") %>%
            # rpair | reaction | compound.x | compound.y | enzyme |  gene
            dplyr::rename(from = compound.x, to = compound.y) %>%
            dplyr::select(c(from, to, gene)) %>%
            dplyr::filter(!duplicated(.)) %>%
            dplyr::filter(!to %in% met.to.filter) %>%
            dplyr::filter(!from %in% met.to.filter) %>%
            dplyr::filter(.$gene %in% entrezIDs)
        }   # length(unique(globalEdgeTable_pre$gene))

        globalGraph <- BioNet::largestComp(
          igraph::graph_from_edgelist(as.matrix(globalEdgeTable_pre[, -3]), directed=FALSE))

        x.1p <- do.call("paste", globalEdgeTable_pre[, -3])
        x.2p <- do.call("paste", igraph::as_data_frame(globalGraph))
        globalEdgeTable <- globalEdgeTable_pre[x.1p %in% x.2p, ]
        # length(unique(globalEdgeTable$gene))

        } else if (organism == "human"){
          load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Hs.eg.gatom.anno.rda"))

          } else {
            stop("Invalid `organism` value: choose either `mouse` or `human`")}
    } else{
        stop("Invalid `method` value: choose either `GAM` or `gatom`")}



  rownames(es) <- entrezIDs
  t <- Biobase::exprs(es)

  gene.exprs <- t

  # GAM
  gene.exprs <- gene.exprs[rownames(gene.exprs) %in% network$rxn2gene$gene, ]
  # gatom
  # gene.exprs <- gene.exprs[rownames(gene.exprs) %in% org.Mm.eg.gatom.anno$gene2enzyme$gene, ]
  # gene.exprs <- gene.exprs[rownames(gene.exprs) %in% globalEdgeTable$gene, ]

  sprintf("There are %.f metabolic genes in the current dataset", dim(gene.exprs)[1])

  zeroSDgenes <- (apply(gene.exprs, 1, sd, na.rm=T) == 0)
  gene.exprs[zeroSDgenes,] <- gene.exprs[zeroSDgenes,] + rnorm(sum(zeroSDgenes) * ncol(gene.exprs), sd = 0.1)

  gene.cor <- cor(t(gene.exprs), use="pairwise.complete.obs")
  gene.cor.dist <- as.dist(1 - gene.cor)
  # gene.cor.dist.m <- as.matrix(gene.cor.dist)

  set.seed(42)
  gene.pam <- cluster::pam(gene.cor.dist, k=initialNumber)

  if(showInitialClustering){
    pheatmap::pheatmap(
      rUtils::normalize.rows(gene.exprs[gene.pam$medoids, ]),
      cluster_rows=F, cluster_cols=F,
      show_rownames=F, show_colnames=T)
  }

  curCenters <- gene.exprs[gene.pam$medoids,]
  curCenters <- limma::avearrays(curCenters, ID=repeats)
  curCenters <- curCenters[, repeats]

  list(
    network = network,
    gene.exprs_orig = t,
    gene.exprs = gene.exprs,
    curCenters = curCenters
  )
}



gamClustering <- function(network = preCluster$network,
                          gene.exprs = preCluster$gene.exprs,
                          curCenters = preCluster$curCenters,
                          repeats = Biobase::pData(es.top12k)$sample,
                          base = 0.4,
                          work.dir = work.dir,
                          method = "GAM", # gatom
                          batch.script = "/home/octopus/R-studio/nclust/sgmwcs-slurm-batch", # synonimus-file-producing
                          # batch.script = "/home/octopus/R-studio/nclust/sgmwcs-slurm-batch-0.9.5", # signal-file-producing # use sgmwcs.batchSolver2 instead of sgmwcs.batchSolver
                          showIntermediateClustering = TRUE,
                          saveSession = TRUE){

  # solver.sb <- sgmwcs.batchSolver2(path.expand(batch.script),
  solver.sb <- sgmwcs.batchSolver(path.expand(batch.script),
                                  nthreads=4,
                                  edges.group.by = "origin", # in GAM, # "gene" in gatom
                                  nodes.group.by = NULL,
                                  group.only.positive = T,
                                  c.size = 50,
                                  timeLimit = 300) # previously used 600

  k <- 1
  revs <- list()

  while (T) {
    while (T) {

      rev <- new.env()
      rev$modules <- list()
      gK1 <- nrow(curCenters)

      rev$centers.pos <- matrix(nrow=gK1,
                                ncol=ncol(gene.exprs),
                                dimnames = list(
                                  paste0("c.pos", seq_len(gK1)),
                                  colnames(gene.exprs)))
      rev$centers.all <- matrix(nrow=gK1,
                                ncol=ncol(gene.exprs),
                                dimnames = list(
                                  paste0("c.all", seq_len(gK1)),
                                  colnames(gene.exprs)))

      dist.to.centers <- 1-cor(t(curCenters), y=t(gene.exprs))
      dist.to.centers[dist.to.centers < 1e-10] <- 0
      idxs <- seq_len(gK1)

      nets <- lapply(idxs, function(i) {

        rUtils::messagef("Processing cluster center %s", i)
        minOther <- pmin(apply(dist.to.centers[-i, ], 2, min), base)
        score <- log2(minOther) - log2(dist.to.centers[i, ])
        rUtils::messagef("Number of genes scored: %s", length(score))

        score[score == Inf] <- 0
        score <- pmax(score, -1000)
        # View(sort(score)) # 84092

        # ----------------------------------------------------------------------
        # if(method == "GAM"){
          es.re.scored <- makeENetworkWithScore(score=score, base=0, network=network)
          es.re.scored$subnet.scored
          # View(igraph::as_data_frame( es.re.scored$subnet.scored))
          # }
          # str(as_data_frame(es.re.scored$subnet.scored))

        # if(method == "gatom"){
        #   EdgeTable <- globalEdgeTable %>%
        #     dplyr::mutate(score = score[gene]) %>%
        #     dplyr::arrange(desc(score)) %>%
        #     dplyr::filter(!duplicated(paste0(from, "__", to)))
        #     # left gene called itself (NB: gene is signal)
        #   igraph::graph_from_data_frame(EdgeTable, directed = F)}
        # -----------------------------------------------------------------------
      })

      cat("Calling solver.sb(nets)\n")
      ms <- solver.sb(nets)
      cat("Done: solver.sb(nets)\n")
      rev$modules <- ms

      for (i in idxs) {

        module <- ms[[i]]
        # for debugging:
        # capture.output(processModule(module, work.dir, sprintf("%s.gene.%s.refined.b%s", "m", i, base)))

        center.pos <- if (ulength(igraph::E(module)[score > 0]$origin) >= 3) { # ORIGIN -> GENE (gatom)
          getCenter(gene.exprs, unique(igraph::E(module)[score > 0]$origin)) # ORIGIN -> GENE (gatom)
        } else {
          curCenters[i, ]
        }
        center.all <- if (ulength(igraph::E(module)$origin) >= 3) { # ORIGIN -> GENE (gatom)
          getCenter(gene.exprs, unique(igraph::E(module)$origin)) # ORIGIN -> GENE (gatom)
        } else {
          curCenters[i, ]
        }
        rev$centers.pos[i, ] <- center.pos
        rev$centers.all[i, ] <- center.all
      }

      rev$centers.pos <- limma::avearrays(rev$centers.pos, ID=repeats)[, repeats]
      rev$centers.all <- limma::avearrays(rev$centers.all, ID=repeats)[, repeats]

      if (showIntermediateClustering) {
        heatmapTable <- rbind(curCenters, rev$centers.pos)[rbind(
          seq_len(gK1),
          seq_len(gK1) + gK1), ]
        pheatmap::pheatmap(
          normalize.rows(heatmapTable),
          cluster_rows=F, cluster_cols=F,
          show_rownames=T, show_colnames=T)
      }

      revs[[k]] <- rev
      revsToCheck <- revs[
        sapply(revs[seq_len(k-1)], function(rev) nrow(rev$centers.pos)) == nrow(rev$centers.pos)]

      diff <- max(abs(rev$centers.pos - curCenters))
      if (length(revsToCheck) > 0) {
        diff <- min(sapply(revsToCheck, function(prevRev) max(abs(rev$centers.pos - prevRev$centers.pos))))
      }

      curCenters <- rev$centers.pos ### AS offered in Aug 2018: TRY rev$centers.all (1)
      rUtils::messagef("Max diff: %s", diff)

      if (diff < 0.01) {
        break
      }
      k <- k + 1
    }

    m.sizes <- sapply(revs[[k]]$modules, function(m) ulength(igraph::E(m)$origin)) # ORIGIN -> GENE (gatom)

    if (all(m.sizes >= 5)) {
      m.sizes <- sapply(revs[[k]]$modules,
                        purrr::compose(diameter, purrr::partial(dualGraph, what="gene")))
      if (all(m.sizes >= 4)) {
        break
      }
    }

    bad <- which(m.sizes == min(m.sizes))
    centersCor <- cor(t(curCenters))
    diag(centersCor) <- NA
    toRemove <- bad[which.max(apply(centersCor, 1, max, na.rm=T)[bad])]
    curCenters <- curCenters[-toRemove, ]
    m.sizes <- m.sizes[-toRemove]
    curCenters <- curCenters[order(m.sizes, decreasing = T),]
  }

  if(saveSession){
    session::save.session(file=paste0(work.dir, "/session.RDa"))
  }

  list(k = k,
       revs = revs,
       curRev = revs[[k]])
}




getGraphs <- function(curRev = gamCluster$curRev,
                      work.dir = work.dir,
                      method = "GAM",
                      solver.paths = "/home/octopus/R-studio/GAM_files/sgmwcs/sgmwcs-0.9.5"){

  solver.gmwcs <- gatom::sgmwcs.solver(solver.paths,
                                       nthreads = detectCores(),
                                       timeLimit=600, # 60 in gatom
                                       # minimize.size = T, # if gatom
                                       group.only.positive = T,
                                       edges.group.by = "origin", # ORIGIN -> GENE (gatom)
                                       # edges.group.by = "gene", # ORIGIN -> GENE (gatom)
                                       nodes.group.by = NULL)

  for (i in seq_along(curRev$modules)) {
    message(i)
    net1 <- curRev$modules[[i]]
    igraph::V(net1)$score <- -1e-2

    m1 <- solver.gmwcs(net1)

    if(method == "gatom"){ # at the moment it has from | to | gene | score only
      igraph::E(m1)$pval <- 0.01
      igraph::E(m1)$label <- org.Mm.eg.gatom.anno$genes$symbol[match(igraph::E(m1)$gene, org.Mm.eg.gatom.anno$genes$gene)]
      igraph::E(m1)$label <- org.Mm.eg.gatom.anno$genes[igraph::E(m1)$gene]$symbol
      igraph::E(m1)$log2FC <- igraph::E(m1)$score
      igraph::V(m1)$label <- met.kegg.db$metabolites$metabolite_name[match(igraph::V(m1)$name, met.kegg.db$metabolites$metabolite)]
      # V(m1)$label <- met.kegg.db$metabolites$metabolite[V(m1)$name]$metabolite_name
      # met.kegg.db$metabolites$metabolite)]
    }
    processModule(module = m1, dir = work.dir, s=1, name = sprintf("%s.%s", "m", i))
    # check "\" in makePdf() here
    # nt -> et (см слак)
  }
}


getHeatmaps <- function(curRev = gamCluster$curRev,
                        work.dir = work.dir){

  # annotation_colors <- list()
  # annotation_colors$geneVar <- setNames(colors, levels(annotation$anno))

  out <- pheatmap::pheatmap(
    normalize.rows(curRev$centers.pos), #[, ord],
    cluster_rows=F, cluster_cols=F,
    file=sprintf("%s/%s.centers.pdf", work.dir, "m"),
    width=25, height=10,
    show_rownames=T, show_colnames=T)
    # annotation = annotation,
    # annotation_colors = annotation_colors)

  # colnames(curRev$centers.pos[, out$tree_col[["order"]]])
  # col_in_order <- colnames(curRev$centers.pos[, out$tree_col[["order"]]])

  for (i in seq_along(curRev$modules)) {
    gs <- unique(igraph::E(curRev$modules[[i]])[order(score)]$origin) # ORIGIN -> GENE (gatom)
    heatmap <- gene.exprs[gs, , drop=F]
    # GAM
    rownames(heatmap) <- reflink[match(rownames(heatmap), Entrez), symbol]
    # gatom
    rownames(heatmap) <- org.Mm.eg.gatom.anno$genes$symbol[
      match(rownames(heatmap), org.Mm.eg.gatom.anno$genes$gene)]

    df <- rUtils::normalize.rows(heatmap)
    # ord <- order(apply(df, 2, mean)) # medium

    pheatmap::pheatmap(
      df,
      cluster_rows=F, cluster_cols=F,
      file=sprintf("%s/%s.%s.genes.png", work.dir, "m", i),
      width=20, height=10,
      show_rownames=T, show_colnames=T)
  }
}


getGeneTables <- function(curRev = gamCluster$curRev,
                          work.dir = work.dir,
                          organism = "mouse"){

  for (i in seq_along(curRev$modules)) {

    t <- GAM::get.edge.attributes(curRev$modules[[i]])[, c("origin", "symbol", "score")] # ORIGIN -> GENE (gatom)
    t <- t[!duplicated(t$origin), ] # ORIGIN -> GENE (gatom)
    colnames(t)[1] <- "Entrez"
    t$cor <- cor(curRev$centers.pos[i, ], t(gene.exprs[t$Entrez, ]))[1,]
    t <- t[order(t$cor, decreasing = T),]

    rUtils::write.tsv(t, file=sprintf("%s/%s.%s.genes.tsv", work.dir, "m", i))

    m1 <- curRev$modules[[i]]
    net <- nets[[i]]

    notInModule <- data.table(Entrez=setdiff(igraph::E(net)$origin, igraph::E(m1)$origin)) # ORIGIN -> GENE (gatom)
    notInModule[, score := igraph::E(net)[match(Entrez, origin)]$score] # ORIGIN -> GENE (gatom)
    notInModule[, symbol := igraph::E(net)[match(Entrez, origin)]$symbol] # ORIGIN -> GENE (gatom), SYMBOL -> GENE (gatom)
    notInModule[, cor := cor(curRev$centers.pos[i, ],
                             t(gene.exprs[Entrez, ]))[1,]]
    notInModule <- notInModule[order(cor, decreasing=T), ]
    rUtils::write.tsv(notInModule[score > 0], file=sprintf("%s/%s.%s.notInModule.genes.tsv", work.dir, "m", i))



    orgdb <- ifelse(organism == "mouse", org.Mm.eg.db::org.Mm.eg.db, org.Hs.eg.db::org.Hs.eg.db)

    notInModule <- data.table(Entrez=rownames(gene.exprs2))
    symbols <- AnnotationDbi::mapIds(orgdb,
                                     keys=notInModule$Entrez,
                                     column="SYMBOL",
                                     keytype="ENTREZID")
    notInModule[, symbol := unname(symbols)]
    notInModule[, cor := cor(curRev$centers.pos[i, ],
                             t(gene.exprs2[Entrez, ]))[1,]]
    notInModule <- notInModule[order(cor, decreasing=T), ]
    rUtils::write.tsv(notInModule[1:300],
                      file=sprintf("%s/%s.%s.complete.genes.tsv", work.dir, "m", i))

    }
}


getGraphs(curRev = gamCluster$curRev,
          work.dir = work.dir)
getHeatmaps(curRev = gamCluster$curRev,
            work.dir = work.dir)
getGeneTables(curRev = gamCluster$curRev,
              work.dir = work.dir,
              organism = "mouse")



annotateModules <-function(universe = Biobase::fData(es.top12k)$entrez,
                           work.dir = work.dir,
                           organism = "mouse"){

  m_files <- list.files(work.dir, "m\\.[0-9]+\\.genes\\.tsv", full.names = T)
  m_filenames <- list.files(work.dir, "m\\.[0-9]+\\.genes\\.tsv")
  if(length(m_files) == 0){
    print("Use `annotateModules()` after `getGeneTables()` only")}

  # str(reactome.db::reactome.db)
  reactomepath <- na.omit(AnnotationDbi::select(reactome.db::reactome.db, universe, "PATHID", "ENTREZID"))
  reactomepath <- split(reactomepath$ENTREZID, reactomepath$PATHID)
  # head(reactomepath)

  keggmodule <- KEGGREST::keggLink("mmu", "module")
  keggmodule <- gsub("mmu:", "", keggmodule)
  names(keggmodule) <- gsub("md:", "", names(keggmodule))
  keggmodule <- split(keggmodule, names(keggmodule))
  # head(keggmodule)
  # keggmodule <- lapply(listoflists, unname)
  # View(keggmodule)

  # keggmdnames <- KEGGREST::keggList("module", "mmu") # 404 after September, 2019
  keggmdnames <- KEGGREST::keggList("module")
  keggmd2name <- data.table::as.data.table(keggmdnames, keep.rownames=T)
  keggmd2name$rn <- gsub("md:", "", keggmd2name$rn)
  data.table::setnames(keggmd2name, c("rn","keggmdnames"), c("PATHID","PATHNAME"))
  keggmd2name$PATHID <- paste0("mmu_", keggmd2name$PATHID)
  # head(keggmd2name)

  keggpathway <- KEGGREST::keggLink("mmu", "pathway")
  keggpathway <- gsub("mmu:", "", keggpathway)
  names(keggpathway) <- gsub("path:", "", names(keggpathway))
  keggpathway <- split(keggpathway, names(keggpathway))
  keggpathway <- lapply(keggpathway, unname)
  # head(keggpathway)
  keggpathnames <- KEGGREST::keggList("pathway", "mmu")
  # head(keggpathnames)

  keggpath2name <- data.table::as.data.table(keggpathnames, keep.rownames=T)
  keggpath2name$rn <- gsub("path:", "", keggpath2name$rn)
  keggpath2name$keggpathnames <- gsub(" - Mus musculus \\(mouse\\)", "",
                                      keggpath2name$keggpathnames)
  data.table::setnames(keggpath2name, c("rn","keggpathnames"), c("PATHID","PATHNAME"))
  # head(keggpath2name)

  reactomepathway2name <- data.table::as.data.table(na.omit(
    AnnotationDbi::select(reactome.db::reactome.db,
                          names(reactomepath),
                          c("PATHNAME"), 'PATHID')))
  reactomepathway2name[, PATHNAME := sub("^[^:]*: ", "", PATHNAME)] # Remove organism prefix
  # head(reactomepathway2name)

  # combine kegg modules and pathways with reactome data
  pathways <- c(reactomepath, keggmodule, keggpathway)
  pathways <- pathways[sapply(pathways, length) >= 10]

  pathways$`5991024` <- NULL
  pathways$`R-MMU-1430728` <- NULL # Metabolism
  pathways$`mmu01100` <- NULL # Metabolic pathways
  pathways$`mmu01200` <- NULL # Carbon metabolism
  pathways$`mmu01230` <- NULL # Biosynthesis of amino acids
  # pathways$`R-MMU-71387` <- NULL # Metabolism of carbohydrates

  pathway2name <- do.call("rbind", list(reactomepathway2name,
                                        keggmd2name,
                                        keggpath2name))

  # pathway2name <- do.call("rbind", list(reactomepathway2name, keggpath2name))

  # fgseaResMain <- fgseaRes[match(mainPathways, pathway)]
  # fgseaResMain[, leadingEdge := lapply(leadingEdge, mapIds, x=org.Mm.eg.db,
  #                                      keytype="ENTREZID", column="SYMBOL")]
  # fwrite(fgseaResMain, file="fgseaResMain.txt", sep="\t", sep2=c("", " ", ""))

  gseaReactome <- function(genes) {
    overlaps <- data.frame(
      q=sapply(sapply(pathways, intersect, genes), length),
      m=sapply(pathways, length),
      n=length(universe)-sapply(pathways, length),
      k=length(genes))

    igenes <- sapply(pathways, intersect, genes)

    for(i in seq_along(igenes)){
      overlaps[i, 5] <- paste(igenes[[i]], collapse=" ")}

    # q-1 because we want probability of having >=q white balls
    pathways.pvals <- with(overlaps,
                           mapply(phyper, q-1, m, n, k, lower.tail = FALSE))
    res <- data.table::data.table(PATHID=names(pathways),
                                  pval=pathways.pvals,
                                  k=overlaps$q,
                                  K=overlaps$m,
                                  genes = overlaps$V5)

    res[, padj := p.adjust(pval, method="BH")]

    res <- merge(res, pathway2name, by="PATHID")
    res <- res[order(pval),]
    res <- res[, c(1, 2, 3, 4, 6, 7, 5)]
  }

  load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rda"))
  # dplyr::glimpse(org.Mm.eg.gatom.anno)
  # it's simple annotation: $genes, $baseId, $gene2enzyme, $mapFrom

  for (i in 1:length(m_files)){
    file <- data.table::fread(m_files[i])
    name <- basename(m_files[i])
    name <- sapply(strsplit(name, ".", fixed = T),"[[", 2)
    out <- gseaReactome(file$Entrez)

    out <- out[padj < 0.05,]

    pz <- sapply(out$genes, function(x) strsplit(x, " "))
    new_pz <- vector("list", length = nrow(out))

    for(e in seq_along(pz) ){
      for(j in seq_along(pz[[e]]) ){
        new_pz[[e]] <-
          append(new_pz[[e]],
          org.Mm.eg.gatom.anno$genes$symbol[which(org.Mm.eg.gatom.anno$genes$gene
                                                  == pz[[e]] [[j]] )] ) } }

    for(a in seq_along(new_pz)){
      out[a, 7] <- paste(new_pz[[a]], collapse=" ")}

    write.tsv(out, file=sprintf("%s/m.%s.pathways_mod.tsv", work.dir, name))
    print(basename(m_files[i]))
  }
}
