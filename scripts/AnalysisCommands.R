# Gene Expression Distribution By cellType
densityPlots(yeast$meta, color_by = 'cellType', dtype = 'sensitivity', main = 'Gene Detection Distribution', save='sensitivityFilteredDensity')
#PCA/tSNE

makePCA_UMAP <- function(cds, cellTypes, counts = F, output_dir='./', color_by = 'cellType')
{
  umap_config= list(n_neighbors = 15, min_dist = 0.5, metric='pearson')
  samples <- row.names(cds$meta[cds$meta$cellType %in% cellTypes,])
  meta <- cds$meta[samples,]
  mat <- cds$tpm$bio
  #mat <- mat[rowSums(mat > 0) > ncol(mat)*0.1,]
  #mat <- mat[rowMeans(mat) >= quantile(rowMeans(mat)[rowMeans(mat) > 0], 0.025) & rowMeans(mat) <= quantile(rowMeans(mat)[rowMeans(mat) > 0], 0.975), ]
  mat <- log(mat+1)
  if(counts){
    ct <- cds$ct$bio
    #ct <- cds$ct$bio[setdiff(row.names(cds$ct$bio), c('Gm26917', 'Lars2')),]
    mat <- DESeq2::varianceStabilizingTransformation(as.matrix(ct))
    #mat <- log(t(t(mat[,samples])/edgeR::calcNormFactors(mat[,samples]), )+1)
    #mat <- log(edgeR::cpm(mat)+1)
    #mat <- mat[rowSums(mat > 0) > ncol(mat)*0.1,]
  }
  plts <- list()
  print(dim(mat))
  #mat <- mat[rowMeans(mat) >= quantile(rowMeans(mat)[rowMeans(mat) > 0], 0.025), ]
  sample_PCA(mat, meta, umap =  T, dimension=2 ,main = "UMAP of Log Expression", labeling = T, point_size = 2, color_by = color_by, umap.config = umap_config)
  ggsave(paste(output_dir, "_uMAP_label.tiff", sep=''), units="in", width=5, height=4, dpi=300, compression = 'lzw')
  plts[['umap']] <- sample_PCA(mat, meta, umap =  T, dimension=2 ,main = "UMAP of Log Expression", labeling = T, point_size = 4, color_by = color_by, umap.config = umap_config)
  ggsave(paste(output_dir, "_uMAP.tiff", sep=''), units="in", width=5, height=4, dpi=300, compression = 'lzw')
  sample_PCA(mat, meta, dimension=2,main = "PCA of Log Expression", labeling = T, point_size = 2, color_by = color_by)
  ggsave(paste(output_dir,"_PCA_label.tiff", sep=''), units="in", width=5, height=4, dpi=300, compression = 'lzw')
  plts[['pca']] <- sample_PCA(mat, meta, dimension=2,main = "PCA of Log Expression", point_size = 4, color_by = color_by)
  ggsave(paste(output_dir,"_PCA.tiff", sep=''), units="in", width=5, height=4, dpi=300, compression = 'lzw')
  return(plts)
}





#sample heatmaps
makeHmap <- function(cds, cellTypes = NULL, output_dir='./', legend = T)
{
  conds <- unique(cds$meta$cellType)
  if(!is.null(cellTypes)){conds <- cellTypes}
  samples <- row.names(cds$meta[cds$meta$cellType %in% conds,])
  samples <- sample(samples)
  meta <- cds$meta[samples,]
  mat <- cds$tpm$bio[,samples]
  #mat <- t(t(mat)/estimateSizeFactorsForMatrix(mat))
  mat <- mat[rowMeans(mat) >= quantile(rowMeans(mat)[rowMeans(mat) > 0], 0.025) & rowMeans(mat) <= quantile(rowMeans(mat)[rowMeans(mat) > 0], 0.975), ]
  
  anno <- data.frame(cellType = as.character(CONDITIONS_NAMES[as.character(meta$cellType)]))
  row.names(anno) <- colnames(mat)
  png(paste(output_dir,'heatmap.png'), units="in", width=5.5, height=5, res=300)
  plts <- sample_heatmap(log(mat+1), distmethod = 'pearson', main='', annotations = anno, fontsize = 12, legend = legend)
  dev.off()
  return(plts)
}

makeDEGmaps <- function(cds, cellTypes = 'all', res, main = '', output_dir='./', individual = T)
{
  if(cellTypes != 'all'){
    samples <- do.call("c", lapply(cellTypes, FUN = function(x) {row.names(cds$meta[cds$meta$cellType %in% x,])}))
  }else{
    samples <- row.names(cds$meta[order(cds$meta$cellType),])
  }
  meta <- cds$meta[samples,]
  mat_cpm <- cpm(cds$ct$bio[,samples])
  mat_ct <- t(t(cds$ct$bio[,samples])/estimateSizeFactorsForMatrix(cds$ct$bio[,samples]))
  mat_tpm <- cds$tpm$bio[,samples]
  mat_tmm <- t(t(mat_ct)/edgeR::calcNormFactors(mat_ct))
  mat_wilcox <-  t(t(mat_tpm)/edgeR::calcNormFactors(mat_tpm))
  maps <- list()
  for(i in names(res$gene_names)){
    if(i %in% c('monocle', 'mast', 'bimod', 'roc', 't', 'LR')){
      mat <- mat_tpm
    }else if(i %in% c('scde')){
      mat <- mat_ct
    }else if(i %in% c('bpsc', 'limma', 'scdd', 'basics')){
      mat <- mat_cpm
    }else if(i == 'wilcox'){
      mat <- mat_wilcox
    }else{
      mat <- mat_tmm
    }
    #mat <- mat_ct
    genes <- res$gene_names[[i]]
    if(length( genes) == 0 ){
      next
    }
    anno <- data.frame(cellType = as.character(CONDITIONS_NAMES[as.character(meta$cellType)]))
    if(length(intersect(genes, row.names(mat))) == 1){
      if(individual){
        fonts = 12
        png(paste(output_dir,i,'.png', sep = '_'), units="in", width=2, height=4, res=300)
      }else{
        fonts = 22
      }
      maps[[i]] <- geneHeatmap(log2(mat[intersect(genes, row.names(mat)),, drop=F]+1), distmethod = 'pearson', clusterGenes = F, annotations = anno, fontsize = fonts, main = i)$gtable
      if(individual){
      dev.off()
      }
      next
    }
    mat <- mat[intersect(genes, row.names(mat)),,drop=F]
    mat <- mat[rowSums(mat) > 0,]
    if(nrow(mat) == 0 ){next}
    if(i != 'esr'){
      mat <- mat[order(res$tpmLog2FC[row.names(mat)]),]
    }
  row.names(anno) <- colnames(mat)
  if(individual){
    fonts = 12
    png(paste(output_dir, i,'png', sep = '.'), units="in", width=1.5, height=4, res=300)
  }else{
    fonts = 22
  }
  maps[[i]] <- geneHeatmap(log2(mat+1), distmethod = 'pearson', clusterGenes = F,annotations = anno, fontsize = fonts, main = i )$gtable
  if(individual){
    dev.off()
  }
  }
  if(!individual){
  plot_g <- cowplot::plot_grid(plotlist = maps, ncol=ceiling(sqrt(length(maps))), nrow = ceiling(sqrt(length(maps))), labels=NULL, label_size = 22)
  title <- cowplot::ggdraw() + 
    cowplot::draw_label(
      main,
      fontface = 'bold', size = 44, 
      x = 0,
      hjust = 0
    ) 
  
  cowplot::plot_grid(
    title, plot_g,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  )+theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )
  ggsave(paste(output_dir, '.png', sep = ''), width = ceiling(sqrt(length(maps)))*8, height = ceiling(sqrt(length(maps)))*9, dpi = 300)
  }
}



#diff_test scde
performPairwiseSCDE <- function(cds, cellTypes='all', batch=NULL, result_name='SCDE', write=FALSE, knn=F, output_dir='./')
{
  out <- tryCatch(
  {
    #PairwiseSCDE <- list()
    #comparisons <- combn(cellTypes, m = 2, simplify = FALSE)
    #for (c in comparisons){
     # scde_results <- scde_diff(cds$ct$bio, cds$meta, compare = c)
      #PairwiseSCDE[[paste(c[1],'_v_',c[2], sep = "")]] <- scde_results
    #}
    knn.model = NULL
    scde.model = NULL
    if(!is.null(cds[[result_name]])){
      if(!is.null(cds[[result_name]][['scde.error.model']])) {scde.model = cds[[result_name]][['scde.error.model']]}
      if(!is.null(cds[[result_name]][['knn.error.model']])) {knn.model = cds[[result_name]][['knn.error.model']]}
    }
    if(knn){result_name <- paste(result_name,'_knn', sep='')}
    cds[[result_name]] <- scde_diff(cds$ct$bio, cds$meta, knn = knn, write_txt=F, compare = cellTypes, knn.model = knn.model, scde.model = scde.model)
    if(write)
    {
      if(knn){knn <- 'knn_'}else{knn <- ''}
      for(c in names(cds[[result_name]])){
      scde_result=cds[[result_name]][[c]]
      write.table(scde_result[scde_result$padj < 0.05, ], file = paste(output_dir,"scde_", knn, c ,".tsv", sep=""), quote = FALSE, sep = "\t")
      }
    }
  return(cds)
  }, error=function(cond){
    message(cond)
    return(cds)
  })
  return(out)
}

yeast <- performPairwiseSCDE(yeast, write = TRUE)
yeast_bc <- performPairwiseSCDE(yeast_bc, write = TRUE)


#figure 3.6
makePlots3.6 <- function(cds, output_file = "./chen_3.6_single_new.pdf", tpm=T, method = 'glm')
{
  
  somePDFPath = output_file
  pdf(file=somePDFPath)
  par(mfrow=c(2,2), cex.lab=1.5, font.main = 2,font.lab = 2,  cex.main = 2, cex.axis=1.5)
  for ( cond in unique(cds$meta$cellType)){
    cells <- row.names(subset(cds$meta, cellType == cond))
  if(tpm){
    failW = scde.failure.probability(models = cds$SCDE$scde.error.model[cells,], counts = cds$ct$bio[,cells])
    ct = cds$tpm$bio[,row.names(subset(cds$meta, cellType == cond))]
    if(method == 'glm'){
      chen_fig_36_1(ct, main=cellTypeS_NAMES[cond],dot_col=DOT_COLOR[cond], median.norm = F, quantile.norm = F)
    }else{
    chen_fig_36_1(ct, dot_col=DOT_COLOR[cond], line_col="#FF000080", main=cellTypeS_NAMES[cond], quantile.norm = F, failW = failW, method = 'scde')
    }
  }else{
    ct = cds$ct$bio[,row.names(subset(cds$meta, cellType == cond))]
    failW = scde.failure.probability(models = cds$SCDE$scde.error.model[cells,], counts = cds$ct$bio[,cells])
    if(method == 'glm'){
      chen_fig_36_1(ct, main=cellTypeS_NAMES[cond],dot_col=DOT_COLOR[cond], median.norm = T, quantile.norm = F)
    }else{
    chen_fig_36_1(ct, dot_col=DOT_COLOR[cond], line_col="#FF000080", main=cellTypeS_NAMES[cond], quantile.norm = T, failW = failW, method = method)
    }
   }
  
  }
  dev.off(); 
}

env2dataframe <- function(env){
  x <- c()
  for (p in ls(env)){

  }
}
#figure3.7
makePlots3.7 <- function(cds, cellTypes,resid, output_dir='./')
{
  env2dataframe <- function(env){
    
  }
  x <- chen_fig_37(ct = cds$tpm$bio, meta = cds$meta, 
              gene2pathway = cds$biocyc_gene2path, 
              cellTypes = cellTypes,
              residual = resid, dir = output_dir)
  dev.off(); 
}

#figure3.8
makePlots3.8 <- function(cds)
{
  chen_fig_38(cds$tpm$bio, cds$meta, cds$PagodaGE$pagoda_go, go=TRUE)
  chen_fig_38(cds$tpm$bio, cds$meta, cds$PagodaGE$pagoda_kegg, kegg=TRUE)
  chen_fig_38(cds$tpm$bio, cds$meta, cds$PagodaGE$pagoda_biocyc, biocyc=TRUE)
  dev.off(); 
}


#Enrichment with pagoda
performPagodaGE <- function(cds, cellTypes='all', batch=NULL, write=FALSE, result_name='PagodaGE',write_only=FALSE, output_dir='./', model=NULL, method = 1)
{
  if(write_only == FALSE)
  {
    PagodaGE <- list()
    PagodaGE[["pagoda_biocyc"]] <- pagoda(ct = cds$ct$bio, meta = cds$meta, env = cds$biocyc.env,  err.model = model, compare = cellTypes, method = method)
    PagodaGE[["pagoda_go"]] <- pagoda(ct = cds$ct$bio, meta = cds$meta, env = cds$go.env,  err.model = model, compare = cellTypes, method = method)
    PagodaGE[["pagoda_kegg"]] <- pagoda(ct = cds$ct$bio, meta = cds$meta, env = cds$kegg.env,  err.model = model, compare = cellTypes, method = method)
    all_env <- new.env()
    #keggenv <- new.env()
    #for(n in names(cds$kegg.env)){
      #keggenv[[n]] <- cds$kegg.env[[n]][[2]]
      #all_env[[paste(n, cds$kegg.env[[n]][[1]])]] <- cds$kegg.env[[n]][[2]]
    #}
    #pagoda_kegg_res <- pagoda_enrichment_result(cds$ct$bio, cds$meta, cellTypes,PagodaGE[["pagoda_kegg"]], keggenv)
    #pagoda_kegg_res <- pagoda(cds$ct$bio,cds$meta, cellTypes,PagodaGE[["pagoda_kegg"]], cds$kegg.env)
    #pagoda_biocyc_res <- pagoda_enrichment_result(cds$ct$bio,cds$meta, cellTypes,PagodaGE[["pagoda_biocyc"]], cds$biocyc.env)
    #agoda_go_res <- pagoda_enrichment_result(cds$ct$bio, cds$meta, cellTypes,PagodaGE[["pagoda_go"]], cds$go.env)
    for (e in c('go.env','biocyc.env', 'kegg.env')){
      for (n in names(cds[[e]])){
      all_env[[n]] <- cds[[e]][[n]]
      }
    }
    PagodaGE[['pagoda_all']] <- pagoda(ct = cds$ct$bio, meta = cds$meta, env = all_env,  err.model = model, compare = cellTypes, method = method)
    cds[[result_name]] <- PagodaGE
  }
  else{
    write=TRUE
  }
  if(write == TRUE)
  {
    write_out_enrichment(cds[[result_name]]$pagoda_go, file = paste(output_dir,"pagoda_go_GE", sep = ''))
    write_out_enrichment(cds[[result_name]]$pagoda_kegg, file = paste(output_dir,"pagoda_kegg_GE", sep = ''))
    write_out_enrichment(cds[[result_name]]$pagoda_biocyc, file = paste(output_dir,"pagoda_biocyc_GE", sep = ''))
    write_out_enrichment(cds[[result_name]]$pagoda_all, file = paste(output_dir,"pagoda_all_GE", sep = ''))
    if(write_only==TRUE)
    {return(cds)}
  }
  return(cds)
}

yeast <- performPagodaGE(yeast, write=T, output_dir = './method_1_', model = yeast_bc$SCDE$knn.error.model, result_name = 'PagodaGE_1', method = 1)
yeast <- performPagodaGE(yeast, write=T, output_dir = './method_2_', model = yeast$SCDE$knn.error.model, result_name = 'PagodaGE_2',method = 2)
yeast <- performPagodaGE(yeast, write=T, output_dir = './method_3_', model = yeast$SCDE$knn.error.model, result_name = 'PagodaGE_3',method = 3)

yeast_ribo <- performPagodaGE(yeast_ribo, write= T, model = yeast_ribo$SCDE$knn.error.model, method = 3, output_dir = './method_3_', result_name = 'PagodaGE_3')
yeast_ribo <- performPagodaGE(yeast_ribo, write= T, model = yeast_ribo$SCDE$knn.error.model, method = 2, output_dir = './method_2_', result_name = 'PagodaGE_2')

#enrichment with goseq/DAVID
performGoseqGE <- function(cds, write=FALSE, DE="scde", output_dir='./')
{
  if (DE == "scde")
  {
    DEG <- "SCDE"
    if(!is.null(cds[["GoseqGE"]]))
    {
      cds[["GoseqGE"]][["scde"]] <- list()
    }
    else
    {
      cds[["GoseqGE"]] <- list(scde=list())
    }
  }
  else if( DE == "monocle")
  {
    DEG <- "monocleDE"
    if(!is.null(cds[["GoseqGE"]]))
    {
      cds[["GoseqGE"]][["monocle"]] <- list()
    }
    else
    {
      cds[["GoseqGE"]] <- list(monocle=list())
    }
  }
  else{
    message("this is not correct")
    exit()
  }
  cds[["GoseqGE"]][[DE]][["goseq_go"]] <- list()
  cds[["GoseqGE"]][[DE]][["goseq_kegg"]] <- list()
  cds[["GoseqGE"]][[DE]][["goseq_biocyc"]] <- list()
  cds[["GoseqGE"]][[DE]][["goseq_all"]] <- list()
  env2df <- function(env)
  {
    do.call(rbind, lapply(names(env), function(x) cbind(env[[x]]$genes, rep(x, length(env[[x]]$genes)))))
  }
  keggcats <- as.data.frame(env2df(cds$kegg.env))
  gocats <- as.data.frame(env2df(cds$go.env))
  biocyccats <- as.data.frame(env2df(yeast$biocyc.env))
  allcats <- rbind(keggcats, gocats, biocyccats)
  for (item in ls(cds[[DEG]]))
  {
    enrichedGO <- gores(cds[[DEG]][[item]], rowMeans(cds$len$bio), gocats)
    enrichedKEGG <- gores(cds[[DEG]][[item]], rowMeans(cds$len$bio), keggcats)
    enrichedBiocyc <- gores(cds[[DEG]][[item]], rowMeans(cds$len$bio), biocyccats)
    enrichedAll <- gores(cds[[DEG]][[item]], rowMeans(cds$len$bio), allcats)
    cds[["GoseqGE"]][[DE]][["goseq_go"]][[item]] <- enrichedGO
    cds[["GoseqGE"]][[DE]][["goseq_kegg"]][[item]] <- enrichedKEGG
    cds[["GoseqGE"]][[DE]][["goseq_biocyc"]][[item]] <- enrichedBiocyc
    cds[["GoseqGE"]][[DE]][["goseq_all"]][[item]] <- enrichedAll
  }
  if(write)
  {
    write_out_enrichment(cds$GoseqGE[[DE]]$goseq_go, file = paste(output_dir,DE,"_goseq_go_GE.txt", sep = ""))
    write_out_enrichment(cds$GoseqGE[[DE]]$goseq_kegg, kegg = cds$kegg.env, file = paste(output_dir,DE,"_goseq_kegg_GE.txt", sep = ""))
    write_out_enrichment(cds$GoseqGE[[DE]]$goseq_biocyc, file = paste(output_dir,DE,"_goseq_biocyc_GE.txt", sep = ""))
    write_out_enrichment(cds$GoseqGE[[DE]]$goseq_all, file = paste(output_dir,DE,"_goseq_all_GE.txt", sep = ""))
  }
  return(cds)
}


#Network Analysis of Gene Covariances
performWGCNA <- function(cds, cellType, covariance=FALSE, enrichment=FALSE)
{
  matr <- yeast$tpm$sep.tpm[[cellType]]
  powers <- c(c(1:10), seq(from = 12, to=50, by=1))
  test=cov(t(yeast$tpm$sep.tpm$-glucose[rowMeans(yeast$tpm$sep.tpm$-glucose) > 0.1,]))
  sft=pickSoftThreshold(test, dataIsExpr = TRUE, powerVector = powers, corFnc = cor, corOptions = list(use = 'p'), networkType = "signed")
  softPower=30
  adj= adjacency(test,type = "signed", power = softPower);
  TOM=TOMsimilarityFromExpr(test,networkType = "signed", TOMType = "signed", power = softPower);
  colnames(TOM) =rownames(TOM) =row.names(test)
  dissTOM=1-TOM
  geneTree = flashClust(as.dist(dissTOM),method="average");
  geneTree$height <- round(geneTree$height, 6)
  plot(geneTree, xlab="", sub="",cex=0.1);
  minModuleSize = 20;
  dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
  table(dynamicMods)
}

#Covariance_Correlation
performCovCor <- function(cds)
{
  testList <- list()
  testList <- noiseHeatmaps(testList, cds, "All_samples",size = 1000)
  testList <- noiseHeatmaps(testList, cds, "isotonic",size = 1000)
  testList <- noiseHeatmaps(testList, cds, "-glucose",size = 1000)
  testList <- noiseHeatmaps(testList, cds, "hypotonic",size = 1000)
  testList <- noiseHeatmaps(testList, cds, "-aa",size = 1000)
  cds$gene_relations <- testList
}


#Gene Noise Enrichment Fisher Test
performGeneNREnrichment <- function(cds, cellTypes=c('-glucose','-aa','hypotonic','isotonic'), result_name='NREnrichement', write_result=FALSE, output_dir='./', diff=T)
{
  NRenrichment <- list()
  keggenv <- new.env()
  all_env <- new.env()
  for(n in names(cds$kegg.env))
  {
    all_env[[n]] <- cds$kegg.env[[n]]
  }
  for (n in names(cds$go.env)){
    all_env[[n]] <- cds$go.env[[n]]
  }
  for (n in names(cds$biocyc.env)){
    all_env[[n]] <- cds$biocyc.env[[n]]
  }
  NRenrichment[['biocycNR']] <- GeneNoiseEnrichment(cds$tpm$bio, cds$meta,cellTypes, cds$biocyc.env, differential = diff)
  NRenrichment[['goNR']] <- GeneNoiseEnrichment(cds$tpm$bio, cds$meta,cellTypes, cds$go.env, differential = diff)
  NRenrichment[['keggNR']] <- GeneNoiseEnrichment(cds$tpm$bio, cds$meta,cellTypes, cds$kegg.env, differential = diff)
  NRenrichment[['allNR']] <- GeneNoiseEnrichment(cds$tpm$bio, cds$meta,cellTypes, all_env, differential = diff)
  cds[[result_name]] <- NRenrichment
  if(write_result){
      for (i in ls(cds[[result_name]])){
        for (item in ls(cds[[result_name]][[i]]))
        {
          df <- cds[[result_name]][[i]][[item]]
          if(nrow(df) > 0){
          if(diff){
            df <- df[order(df$Wilcox.SR.adj),]
            df$Wilcox.SR <- formatC(df$Wilcox.SR, format = "e", digits = 2)
            df$Wilcox.RS <- formatC(df$Wilcox.RS, format = "e", digits = 2)
            df$Wilcox.SR.adj <- formatC(df$Wilcox.SR.adj, format = "e", digits = 2)
            df$Wilcox.RS.adj <- formatC(df$Wilcox.RS.adj, format = "e", digits = 2)
          }else{
            df <- df[order(rowMins(as.matrix(df[,c('negPadj', 'posPadj')]))),]
            df$negPval <- formatC(df$negPval, format = "e", digits = 2)
            df$posPval <- formatC(df$posPval, format = "e", digits = 2)
            df$negPadj <- formatC(df$negPadj, format = "e", digits = 2)
            df$posPadj <- formatC(df$posPadj, format = "e", digits = 2)
          }
          }
          write.table(df, file = paste(output_dir, i,'_',item,'_', result_name,'.tsv', sep=''),  quote = FALSE, sep = "\t",  row.names = T)
        }
      }
  }
  return(cds)
}



CompareDEMethods <- function(SCDE, MONOCLE, MAST, pthreshold = 0.05){
  for (c in c('isotonic_v_-glucose', 'isotonic_v_-aa', '-glucose_v_-aa', 'isotonic_v_hypotonic')){
    print(c)
    scde_res <- row.names(subset(SCDE[[c]], padj <= pthreshold ))
    monocle_res <- row.names(subset(MONOCLE[[c]], qval <= pthreshold ))
    mast_res <- subset(MAST[[c]][['DEfull']], fdr <= pthreshold )$primerid
    print(paste('Number of significant genes using', 'SCDE:', length(scde_res)))
    print(paste('Number of significant genes using', 'Monocle:', length(monocle_res)))
    print(paste('Number of significant genes using', 'MAST:', length(mast_res)))
    print(paste('Shared Genes between', 'Monocle and SCDE:', length(intersect(scde_res, monocle_res))))
    print(paste('Shared Genes between', 'Monocle and MAST:', length(intersect(mast_res, monocle_res))))
    print(paste('Shared Genes between', 'MAST and SCDE:', length(intersect(mast_res, scde_res))))
    print(paste('Shared Genes between', 'Monocle and SCDE and MAST:', length(intersect(intersect(scde_res, monocle_res), mast_res))))
    cat('\n')
  }
  
}











