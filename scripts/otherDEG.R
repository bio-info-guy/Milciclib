
SeuratDEG <- function(obj, DEmethod = 'wilcox', normalization = 'LogNormalize', reference = 'isotonic', min.pct = 0.5, min.cells = 5, min.features = 1000, latent.vars = 'sensitivity'){
  mat <- obj$tpm$bio
  meta <- obj$meta[colnames(mat),]
  meta$sensitivity <- colSums(mat[,row.names(meta)] != 0)
  meta <- meta[colnames(mat),]
  seurat_obj <- CreateSeuratObject(counts = mat, project = "sce", min.cells = min.cells, min.features = min.features, meta.data = meta[colnames(mat),])
  seurat_obj@assays[['raw_ct']] <- CreateAssayObject(obj$ct$bio, min.cells = min.cells, min.features = min.features)
  seurat_obj@assays$raw_ct@key <- "raw_ct_"
  Idents(seurat_obj) <- meta[names(Idents(seurat_obj)),]$cellType
  if(normalization == 'sctransform'){
    seurat_obj <- SCTransform(
      seurat_obj, seed.use = NULL)
  }else{
    if(DEmethod == 'wilcox'){
        seurat_obj@assays[['RNA']]@data <- t(t(seurat_obj@assays[['RNA']]@data)/edgeR::calcNormFactors(seurat_obj@assays[['RNA']]@data))
        #seurat_obj <- NormalizeData(seurat_obj, normalization.method = normalization, scale.factor = 1e6, assay = 'RNA', verbose = F)
    }else{
        seurat_obj <- NormalizeData(seurat_obj, normalization.method = normalization, scale.factor = 1e6, assay = 'RNA', verbose = F)
    }
  }
  assay2use = ifelse(DEmethod == 'DESeq2', 'raw_ct', 'RNA')
  if(DEmethod == 'DESeq2'){
     names_underscore <-row.names(obj$ct$bio) 
    }else{
      names_underscore <-row.names(mat)
  }
  
  names(names_underscore)  <- str_replace(names_underscore, pattern = '_', replacement = '-')
  test <- FindMarkers(seurat_obj, 
                      assay = assay2use , 
                      logfc.threshold = 0, 
                      min.pct = min.pct , 
                      ident.1 = setdiff(unique(Idents(seurat_obj)), c(reference))[1], 
                      ident.2 = reference,  
                      test.use = DEmethod, 
                      slot = 'data', 
                      latent.vars = latent.vars, 
                      verbose = F)

  #print(names_underscore[row.names(test)])
  row.names(test) <- names_underscore[row.names(test)]
  if(DEmethod == 'roc'){
    test$Log2FC <- test$avg_log2FC
    return(subset(test, power >= 0))
  }else if(DEmethod == 'wilcox'){
    mat <- seurat_obj@assays[['RNA']]@data
    MeanExp <- tapply(colnames(mat), as.character(meta[colnames(mat),]$cellType), function(x){rowMeans(mat[,x])})
    test$Log2FC <- log2((MeanExp[[setdiff(unique(Idents(seurat_obj)), c(reference))[1]]]+1e-16)/(MeanExp[[reference]]+1e-16))[row.names(test)]
    return(subset(test, p_val_adj <= 1))
  }else{
    test$p_val_adj <- p.adjust(test$p_val, 'BH')
    test$Log2FC <- test$avg_log2FC
    return(subset(test, p_val_adj <= 1))
  }
}

deg_basics <- function(obj, control='isotonic', min.pct = 0.5){
  target <- setdiff(unique(obj$meta$cellType), c(control))
  ctrl_mat <- obj$ct$bio[,row.names(subset(obj$meta, cellType == control))]
  ctrl_mat0 <- ctrl_mat[rowSums(ctrl_mat != 0) > ncol(ctrl_mat)*min.pct,]
  trgt_mat <- obj$ct$bio[,row.names(subset(obj$meta, cellType == target))]
  trgt_mat0 <- trgt_mat[rowSums(trgt_mat != 0) > ncol(trgt_mat)*min.pct,]
  genes <- union(row.names(ctrl_mat0), row.names(trgt_mat0))
  ctrl_mat <- ctrl_mat[genes,]
  trgt_mat <- trgt_mat[genes,]
  sce_control <- SingleCellExperiment(assays = list(counts = ctrl_mat), colData = data.frame(BatchInfo = sample(c(rep(1, floor(ncol(ctrl_mat)/2)), rep(2, ncol(ctrl_mat)-floor(ncol(ctrl_mat)/2))), replace = F)))
  sce_target <- SingleCellExperiment(assays = list(counts = trgt_mat), colData = data.frame(BatchInfo = sample(c(rep(1, floor(ncol(trgt_mat)/2)), rep(2, ncol(trgt_mat)-floor(ncol(trgt_mat)/2))), replace = F)))
  #sce_control <- SingleCellExperiment(assays = list(counts = ctrl_mat), colData = data.frame(BatchInfo = colSums(ctrl_mat > 0)))
  #sce_target <- SingleCellExperiment(assays = list(counts = trgt_mat), colData = data.frame(BatchInfo = colSums(trgt_mat > 0)))
  
  ctrl_chain <- BASiCS_MCMC(Data = sce_control, N = 20000, Thin = 20, Burn = 10000, WithSpikes = F, 
                       PrintProgress = FALSE, Regression = TRUE)
  trgt_chain <- BASiCS_MCMC(Data = sce_target, N = 20000, Thin = 20, Burn = 10000, WithSpikes = F,
                       PrintProgress = FALSE, Regression = TRUE)
  Test <- BASiCS_TestDE(Chain1 = trgt_chain, Chain2 = ctrl_chain,
                        GroupLabel1 = target, GroupLabel2 = control,
                        EpsilonM = log2(1.25), EpsilonD = log2(1.25),
                        EFDR_M = 0.05, EFDR_D = 0.05,
                        Offset = TRUE, PlotOffset = F, Plot = F)

  res <- Test@Results$Mean@Table
  res$diff <- !Test@Results$Mean@Table$ResultDiffMean %in% c('NoDiff','ExcludedLowESS')
  res <- data.frame(res, row.names = 1)
  res$Log2FC <- res$MeanLog2FC
  return(res)
}

deg_basics <- function(obj, control='isotonic', min.pct = 0.5, spike = T){
  target <- setdiff(unique(obj$meta$cellType), c(control))
  ctrl_mat <- obj$ct$bio[,row.names(subset(obj$meta, cellType == control))]
  ctrl_mat0 <- ctrl_mat[rowSums(ctrl_mat != 0) > ncol(ctrl_mat)*min.pct,]
  trgt_mat <- obj$ct$bio[,row.names(subset(obj$meta, cellType == target))]
  trgt_mat0 <- trgt_mat[rowSums(trgt_mat != 0) > ncol(trgt_mat)*min.pct,]
  genes <- union(row.names(ctrl_mat0), row.names(trgt_mat0))
  ctrl_mat <- ctrl_mat[genes,]
  trgt_mat <- trgt_mat[genes,]
  if(is.null(spike)){
    sce_control <- SingleCellExperiment(assays = list(counts = ctrl_mat), colData = data.frame(BatchInfo = sample(c(rep(1, floor(ncol(ctrl_mat)/2)), rep(2, ncol(ctrl_mat)-floor(ncol(ctrl_mat)/2))), replace = F)))
    sce_target <- SingleCellExperiment(assays = list(counts = trgt_mat), colData = data.frame(BatchInfo = sample(c(rep(1, floor(ncol(trgt_mat)/2)), rep(2, ncol(trgt_mat)-floor(ncol(trgt_mat)/2))), replace = F)))
    withSpike = F
  }
  else{
    spike <- obj$ct$spike
    spike <- spike[rowSums(spike > 0) > 10,]
    ctrl_cts <- as.matrix(rbind(ctrl_mat, spike[,colnames(ctrl_mat)]))
    trgt_cts <- as.matrix(rbind(trgt_mat, spike[,colnames(trgt_mat)]))
    Tech = c(rep(F, nrow(ctrl_mat)), rep(T, nrow(spike)))
    SpikeInfo = read.csv('./dataset/ERCC_conc.tsv', header = T, sep = '\t')
    colnames(SpikeInfo) <- c('SpikeID', 'SpikeInput')
    SpikeInfo$SpikeInput <- round(SpikeInfo$SpikeInput*5)
    SpikeInfo <- SpikeInfo[match(row.names(spike), SpikeInfo[,'SpikeID']),]
    
    sce_control <- newBASiCS_Data(ctrl_cts, Tech, SpikeInfo)
    sce_target <- newBASiCS_Data(trgt_cts, Tech, SpikeInfo)
    withSpike = T
  }
  #sce_control <- SingleCellExperiment(assays = list(counts = ctrl_mat), colData = data.frame(BatchInfo = colSums(ctrl_mat > 0)))
  #sce_target <- SingleCellExperiment(assays = list(counts = trgt_mat), colData = data.frame(BatchInfo = colSums(trgt_mat > 0)))
  
  ctrl_chain <- BASiCS_MCMC(Data = sce_control, N = 20000, Thin = 20, Burn = 10000, WithSpikes = withSpike, 
                            PrintProgress = FALSE, Regression = TRUE)
  trgt_chain <- BASiCS_MCMC(Data = sce_target, N = 20000, Thin = 20, Burn = 10000, WithSpikes = withSpike,
                            PrintProgress = FALSE, Regression = TRUE)
  Test <- BASiCS_TestDE(Chain1 = trgt_chain, Chain2 = ctrl_chain,
                        GroupLabel1 = target, GroupLabel2 = control,
                        EpsilonM = log2(1.5), EpsilonD = log2(1.5),
                        EFDR_M = 0.05, EFDR_D = 0.05,
                        Offset = TRUE, PlotOffset = F, Plot = T, MinESS = 10)
  
  res <- Test@Results$Mean@Table
  res$diff <- !Test@Results$Mean@Table$ResultDiffMean %in% c('NoDiff','ExcludedLowESS')
  res <- data.frame(res, row.names = 1)
  res$Log2FC <- res$MeanLog2FC
  return(res)
}



deg_scdd <- function(obj, control='isotonic', permutations = 0, min.pct = 0.5){
  target <- setdiff(unique(obj$meta$cellType), c(control))
  ctrl_mat <- obj$ct$bio[,row.names(subset(obj$meta, cellType == control))]
  ctrl_mat <- ctrl_mat[rowSums(ctrl_mat != 0) > ncol(ctrl_mat)*min.pct,]
  trgt_mat <- obj$ct$bio[,row.names(subset(obj$meta, cellType == target))]
  trgt_mat <- trgt_mat[rowSums(trgt_mat != 0) > ncol(trgt_mat)*min.pct,]
  genes <- union(row.names(ctrl_mat), row.names(trgt_mat))
  mat <-  obj$ct$bio[genes,]
  
  sce <- SingleCellExperiment(assays=list(counts=mat),
                              colData=data.frame(condition = obj$meta[colnames(mat),]$cellType))
  sce <- preprocess(sce, zero.thresh = min.pct, median_norm = T)
  prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
  
  if(permutations > 0){
    res <- list(); i <- 1; j <- 1
    pparam <- BiocParallel::MulticoreParam(workers = 7, progressbar = TRUE)
    while(i <= nrow(sce)){
      print(paste('running', i, 'to', i+999, 'genes', sep = ' '))
      scDatExSim <- scDD(sce[i:min((i+999), nrow(sce)),], prior_param=prior_param, testZeroes=T, adjust.perms = T, permutations = permutations, param=pparam)
      gc()
      res[[j]] <- scDD::results(scDatExSim)
      i <- i+1000; j <- j+1
    }
    final <- do.call(rbind, res)
  }else{
    final <- scDD::results(scDD(sce, prior_param=prior_param, testZeroes=T,  permutations = 0, param = BiocParallel::MulticoreParam(workers=1)))
    gc()
  }
  MeanExp <- tapply(colnames(mat), as.character(obj$meta[colnames(mat),]$cellType), function(x){rowMeans(mat[,x])})
  final$Log2FC <- log2((MeanExp[[target]]+1e-16)/(MeanExp[[control]]+1e-16))[row.names(final)]
  final$diff <- final$combined.pvalue.adj < 0.05 & final$DDcategory != 'NS'
  return(final)
}

deg_bpsc <- function(obj, control='isotonic' , min.pct = 0.5){
  mat <- edgeR::cpm(obj$ct$bio)
  target <- setdiff(unique(obj$meta$cellType), c(control))
  ctrl_mat <- mat[,row.names(subset(obj$meta, cellType == control))]
  ctrl_mat <- ctrl_mat[rowSums(ctrl_mat != 0) > ncol(ctrl_mat)*min.pct,]
  trgt_mat <- mat[,row.names(subset(obj$meta, cellType == target))]
  trgt_mat <- trgt_mat[rowSums(trgt_mat != 0) > ncol(trgt_mat)*min.pct,]
  genes <- union(row.names(ctrl_mat), row.names(trgt_mat))
  mat <- mat[genes,]
  meta <- obj$meta[colnames(mat),]
  meta$sensitivity <- colSums(mat[,row.names(meta)] != 0)
  meta <- meta[colnames(mat),]
  group <- as.character(meta$cellType)
  controlIds <- which(group == control)
  design <- model.matrix(~group+meta$sensitivity)
  coef <- 2
  res=BPglm(data=mat, controlIds=controlIds, design=design, coef=coef, estIntPar=F, useParallel=TRUE) 
  results <- cbind(res$PVAL, res$TVAL, p.adjust(res$PVAL, method = 'BH'))
  row.names(results) <- names(res$PVAL)
  colnames(results) <- c('pvalue', 'tvalue', 'p_adj_val')
  res <- as.data.frame(results)
  res$diff <- res$p_adj_val < 0.05
  MeanExp <- tapply(colnames(mat), as.character(meta[colnames(mat),]$cellType), function(x){rowMeans(mat[,x])})
  res$Log2FC <- log2((MeanExp[[target]]+1e-16)/(MeanExp[[control]]+1e-16))[row.names(res)]
  return(res)
}

deg_limma <- function(obj, control='isotonic', trend= F, min.pct = 0.5,min_cell_grp=2, min_cell = 5){
  ct <- obj$ct$bio
  meta <- obj$meta
  ct <- ct[rowSums(ct > 0) > ncol(ct)*0,]
  
  target <- setdiff(unique(meta$cellType), c(control))
  ctrl_mat <- ct[,row.names(subset(meta, cellType == control))]
  ctrl_mat <- ctrl_mat[rowSums(ctrl_mat > 0) > max(ncol(ctrl_mat)*min.pct, min_cell_grp),]
  trgt_mat <- ct[,row.names(subset(meta, cellType == target))]
  trgt_mat <- trgt_mat[rowSums(trgt_mat > 0) > max(ncol(trgt_mat)*min.pct, min_cell_grp),]
  genes <- union(row.names(ctrl_mat), row.names(trgt_mat))
  
  
  genes <- intersect(row.names(ct)[which(rowSums(ct > 0) > min_cell)], genes)
  filt_genes <- setdiff(row.names(ct), genes)
  
  mat <- ct[genes,]
  meta <- obj$meta[colnames(mat),]
  meta$sensitivity <- colSums(mat[,row.names(meta)] != 0)
  meta <- meta[colnames(mat),]
  group <- as.character(meta$cellType)
  group <- relevel(as.factor(group), control)
  dge <- DGEList(counts=mat)
  design <- model.matrix(~group+meta$sensitivity)
  if(trend){
    logCPM <- cpm(dge, log=TRUE, prior.count=3)
    fit <- lmFit(logCPM, design)
    fit <- eBayes(fit, trend=TRUE)
  }else{
    dge <- calcNormFactors(dge)
    v <- voom(dge, design, plot=TRUE)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
  }
  res <- topTable(fit, coef=2, number=30000, p.value = 1)
  res$diff <- res$adj.P.Val < 0.05
  res$Log2FC <- res$logFC
  return(res)
}

deg_edgeR <- function(obj, control='isotonic', min.pct = 0.5, min_cell_grp = 2, min_cell =5){
  ct <- obj$ct$bio
  meta <- obj$meta
  ct <- ct[rowSums(ct > 0) > ncol(ct)*0,]
  
  target <- setdiff(unique(meta$cellType), c(control))
  ctrl_mat <- ct[,row.names(subset(meta, cellType == control))]
  ctrl_mat <- ctrl_mat[rowSums(ctrl_mat > 0) > max(ncol(ctrl_mat)*min.pct, min_cell_grp),]
  trgt_mat <- ct[,row.names(subset(meta, cellType == target))]
  trgt_mat <- trgt_mat[rowSums(trgt_mat > 0) > max(ncol(trgt_mat)*min.pct, min_cell_grp),]
  genes <- union(row.names(ctrl_mat), row.names(trgt_mat))
  
  
  genes <- intersect(row.names(ct)[which(rowSums(ct > 0) > min_cell)], genes)
  filt_genes <- setdiff(row.names(ct), genes)
  
  mat <- ct[genes,]
  meta <- obj$meta[colnames(mat),]
  meta$sensitivity <- colSums(mat[,row.names(meta)] != 0)
  meta <- meta[colnames(mat),]
  group <- as.character(meta$cellType)
  group <- relevel(as.factor(group), control)
  dge <- DGEList(counts=mat, group = group)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~group+meta$sensitivity)
  dge <- estimateDisp(dge,design)
  fit <- glmQLFit(dge,design)
  qlf <- glmQLFTest(fit,coef=2) 
  res <- topTags(qlf, n = 30000, p.value = 1)$table
  res$diff <- res$FDR < 0.05
  res$Log2FC <- res$logFC
  return(res)
}

mast_diff <- function(obj = NULL, plot=F, ct = NULL, meta = NULL, FCThresh = log2(1.25), normFactor = NULL, control = 'mouseEgg', tpm=T, 
                      freq = 0.5, max_thres = 3, bin_by = 'median', nbins = 20, min_per_bin = 30, 
                      correct_wild_coef = T, corr_det = T, min_cell_grp = 3, min_cell = 6, include_filt_as_NA = F){
  #ct <- if(tpm){obj$tpm$bio}else{t(t(obj$tpm$bio)/edgeR::calcNormFactors(obj$tpm$bio))}
  res = list()
  if(!is.null(obj)){
    ct <- if(tpm){obj$tpm$bio}else{
      obj$ct$bio
    }
    meta <- obj$meta[colnames(ct),]
  }
  if(!tpm){
    if(!is.null(normFactor)){
      ct <- t(t(ct)/normFactor)
    }else{
      ct <- edgeR::cpm(ct)
    }
    
  }
  ct <- ct[rowSums(ct > 0) > ncol(ct)*0,]
  
  target <- setdiff(unique(meta$cellType), c(control))
  ctrl_mat <- ct[,row.names(subset(meta, cellType == control))]
  ctrl_mat <- ctrl_mat[rowSums(ctrl_mat > 0) > max(ncol(ctrl_mat)*freq, min_cell_grp),]
  trgt_mat <- ct[,row.names(subset(meta, cellType == target))]
  trgt_mat <- trgt_mat[rowSums(trgt_mat > 0) > max(ncol(trgt_mat)*freq, min_cell_grp),]
  genes <- union(row.names(ctrl_mat), row.names(trgt_mat))
  
  
  genes <- intersect(row.names(ct)[which(rowSums(ct > 0) > min_cell)], genes)
  filt_genes <- setdiff(row.names(ct), genes)
  ct <- ct[genes,]
  print(length(filt_genes))

  print(dim(ct))
  gene_f<- data.frame(row.names = row.names(ct), features = row.names(ct))
  freq_expressed <- freq
  FCTHRESHOLD <- FCThresh
  sca0 <- FromMatrix(as.matrix(log2(ct+1)), meta, gene_f)
  if(nbins > 0){
    thres <- thresholdSCRNACountMatrix(assay(sca0)[rowMedians(assay(sca0)) < max_thres,], conditions = meta$cellType, nbins = nbins, min_per_bin = min_per_bin, bin_by = bin_by)
    thres_ct <- assay(sca0)
    thres_ct[rowMedians(assay(sca0)) < max_thres,] <- thres$counts_threshold
    #thres_ct[thres$original_data > 2] <- thres$original_data[thres$original_data > 2]
    assays(sca0) <- list(thresh=thres_ct, tpm=assay(sca0))
    if(plot){
      par(mfrow=c(ceiling(sqrt(nbins)), ceiling(sqrt(nbins))))
      plot(thres)
    }
  }else{
    assays(sca0) <- list(thresh=assay(sca0), tpm=assay(sca0))
  }
  
  #assays(sca0) <- list(thresh=assay(sca0), tpm=assay(sca0))
  #expressed_genes <- freq(sca0) > freq_expressed
  #sca0 <- sca0[expressed_genes,]
  #print(dim(sca0))
  #res['sca'] <- sca0
  comparisons = list(comp = c(control, target))
  for (cnd in comparisons){
    sca <- sca0[,colData(sca0)$cellType %in% cnd]
    cond<-factor(colData(sca)$cellType)
    if("mouseEgg" %in% cnd){
      cond<-relevel(cond,"mouseEgg")
      other <- paste('cellType',cnd[which(cnd != "mouseEgg")], sep = '')
    }else{
      cond<-relevel(cond,"ratEgg")
      if(cnd[2] == 'ratEgg'){cnd <- cnd[c(2,1)]}
      other <- paste('cellType', cnd[which(cnd != "ratEgg")], sep = '')
    }
    colData(sca)$cellType<-cond
    if(corr_det){
      zlmCond <- zlm(~cellType + sensitivity, sca, useContinuousBayes=TRUE)
    }else{
      zlmCond <- zlm(~cellType, sca, useContinuousBayes=TRUE)
    }
    summaryCond <- summary(zlmCond, doLRT = other)
    summaryDt <- summaryCond$datatable
    #return(summaryDt)
    fcHurdle <- merge(summaryDt[contrast==other & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast==other & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
    #print(dim(fcHurdle))
    fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    fcHurdleSig <- merge(fcHurdle, data.table::as.data.table(mcols(sca)), by='primerid')
    data.table::setorder(fcHurdleSig, fdr)
    row.names(fcHurdleSig) <- fcHurdleSig$primerid
    fcHurdleSig <- data.frame(fcHurdleSig, row.names = 1)
    if(correct_wild_coef){
      
      con <- subset(data.frame(summaryDt), component == 'C' & contrast == summaryDt$contrast[1])
      con <- data.frame(con, row.names = con$primerid)
      fc <- subset(data.frame(summaryDt), component == 'logFC' & contrast == summaryDt$contrast[1])
      fc <- data.frame(fc, row.names = fc$primerid)
      sub_fc <- data.frame(row.names = fc$primerid, primerid = fc$primerid, coef_diff = abs(con$coef - fc$coef), z = abs(con$z) > abs(fc$z))
      
      congenes <- row.names(subset(sub_fc, coef_diff > 0.01 & z))
      fcHurdleSig[, 'orig_logfc'] <- fc[row.names( fcHurdleSig), 'coef']
      fcHurdleSig[, 'con_logfc'] <- con[row.names( fcHurdleSig), 'coef']
      fcHurdleSig[congenes,c('coef', 'ci.hi', 'ci.lo')] <- con[congenes, c('coef', 'ci.hi', 'ci.lo')]
    }
    
    fcHurdleSig$Log2FC <- fcHurdleSig$coef
    fcHurdle <- data.frame(fcHurdle, row.names = 1)
    fcHurdle$Log2FC <- fcHurdle$coef
    if(include_filt_as_NA){
      filt_gene_res <- matrix(nrow = length(filt_genes), ncol = ncol(fcHurdleSig))
      row.names(filt_gene_res) <- filt_genes
      colnames(filt_gene_res) <- colnames(fcHurdleSig)
      filt_gene_res <- data.frame(filt_gene_res)
      filt_gene_res[,1] <- 1
      filt_gene_res[,'Log2FC'] <- 0
      fcHurdleSig <- rbind( fcHurdleSig, filt_gene_res)
      fcHurdleSig$fdr <- p.adjust(fcHurdleSig[,1], 'BH')
      fcHurdleSig[is.na(fcHurdleSig$Log2FC),'Log2FC'] <- 0
    }
    res[[paste(cnd[1],'_v_',cnd[2], sep = "")]] <- list(DESig=fcHurdleSig, model=zlmCond, DEfull=fcHurdle, samples=colData(sca0)$cellType %in% cnd, data = sca0, summaryDt = summaryCond)
  }
  
  return(res)
}

deg_scde <- function(obj = NULL, ct = NULL, meta = NULL, control = 'isotonic', scde.model = T,min_cell_grp = 3, min_cell = 6, min.pct = 0.5){
  
  if(!is.null(obj)){
    mat <- round(obj$ct$bio)
    mat <- apply(mat,2,function(x) {storage.mode(x) <- 'integer'; x})
    meta <- obj$meta
  }else{
    mat <- round(ct)
    mat <- apply(mat,2,function(x) {storage.mode(x) <- 'integer'; x})
    meta <- meta
  }
  target <- setdiff(unique(meta$cellType), c(control))
  ctrl_mat <- mat[,row.names(subset(meta, cellType == control))]
  ctrl_mat <- ctrl_mat[rowSums(ctrl_mat != 0) > max(ncol(ctrl_mat)*min.pct, min_cell_grp),]
  trgt_mat <- mat[,row.names(subset(meta, cellType == target))]
  trgt_mat <- trgt_mat[rowSums(trgt_mat != 0) > max(ncol(trgt_mat)*min.pct, min_cell_grp),]
  genes <- union(row.names(ctrl_mat), row.names(trgt_mat))
  genes <- intersect(row.names(mat)[which(rowSums(mat > 0) > min_cell)], genes)
  mat <- mat[genes,c(row.names(subset(meta, cellType == control)), row.names(subset(meta, cellType == target)))]
  group <- as.character(meta[colnames(mat),]$cellType)
    #sg <-factor(meta[meta$cellType %in% compare,]$cellType)
    #names(sg) <- row.names(meta[meta$cellType %in% compare,])
  out <- tryCatch(
      {
        sg <- factor(group)
        sg <- relevel(sg, control)
        names(sg) <- colnames(mat)
        mat <- mat[,names(sg)]
        cd <- clean.counts(mat, min.lib.size = 250, min.reads = 1, min.detected = round(ncol(mat)*0))
        sg <- sg[colnames(cd)]
        cat("finished cleaning counts\n")
        if (scde.model){
          print(cd)
          print(sg)
          fim <- scde.error.models(counts=cd, groups = sg, n.cores=1, threshold.segmentation = T, verbose = 1, min.size.entries = 2000)
        }else{
          fim <- knn.error.models(cd, k = max(ncol(cd)/4, 5), n.cores = 1, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 0)
            } 
        cat("finished error models\n")
        valid.cells <- row.names(fim[fim$corr.a>0,])
        sg <- sg[valid.cells]
        cd <- cd[,valid.cells]
        fim <- fim[valid.cells,]
        prior <- scde.expression.prior(models = fim, counts=cd, length.out = 1000, show.plot = F)
        cat("finished expression priors\n")
        cat("begin calculating expression differences\n")
        ediff <- scde.expression.difference(fim, cd, prior, groups  =  sg, n.randomizations  =  1000, n.cores  =  1, verbose  =  1)
        p.values <- 2*pnorm(abs(ediff$Z),lower.tail=F) # 2-tailed p-value
        p.values.adj <- p.adjust(p.values, method="BH") # Adjusted to control for FDR
        res <- data.frame(row.names=rownames(ediff),lower_bound=ediff[,1], log2FC=ediff[,2], upper_bound=ediff[,3],p.val=p.values, padj=p.values.adj)
        res$diff <- res$padj < 0.05
        res
        }, error=function(con){
          message(con)
          return(fim)
        })
    out$Log2FC <- -out$log2FC
    out$gene_short_name <- row.names(out)
    return(out)
}


min.pct.filter <- function(mat, meta, control, min.pct = 0.3){
  target <- setdiff(unique(meta$cellType), c(control))
  ctrl_mat <- mat[,row.names(subset(meta, cellType == control))]
  ctrl_mat <- ctrl_mat[rowSums(ctrl_mat != 0) > ncol(ctrl_mat)*min.pct,]
  trgt_mat <- mat[,row.names(subset(meta, cellType == target))]
  trgt_mat <- trgt_mat[rowSums(trgt_mat != 0) > ncol(trgt_mat)*min.pct,]
  genes <- union(row.names(ctrl_mat), row.names(trgt_mat))
  mat <- mat[genes, c(row.names(subset(meta, cellType == control)), row.names(subset(meta, cellType == target)))]
  return(mat)
}


deg_monocle <- function(obj, control = 'isotonic', min.pct = 0.5, tpm2abs = T, min_cell_grp = 2, min_cell =5){
  require(monocle)
  if(tpm2abs){
  ct <- obj$tpm$bio
  }else{
    ct <- obj$ct$bio
  }
  
  meta <- obj$meta
  ct <- ct[rowSums(ct > 0) > ncol(ct)*0,]
  
  target <- setdiff(unique(meta$cellType), c(control))
  ctrl_mat <- ct[,row.names(subset(meta, cellType == control))]
  ctrl_mat <- ctrl_mat[rowSums(ctrl_mat > 0) > max(ncol(ctrl_mat)*min.pct, min_cell_grp),]
  trgt_mat <- ct[,row.names(subset(meta, cellType == target))]
  trgt_mat <- trgt_mat[rowSums(trgt_mat > 0) > max(ncol(trgt_mat)*min.pct, min_cell_grp),]
  genes <- union(row.names(ctrl_mat), row.names(trgt_mat))
  
  
  genes <- intersect(row.names(ct)[which(rowSums(ct > 0) > min_cell)], genes)
  filt_genes <- setdiff(row.names(ct), genes)
  
  mat <- ct[genes,]
  mat <- mat[genes, c(row.names(subset(obj$meta, cellType == control)), row.names(subset(obj$meta, cellType == target)))]
  meta <- obj$meta[colnames(mat),]
  meta$sensitivity <- colSums(mat[,row.names(meta)] != 0)
  feature <- obj$features
  
  #feature$external_gene_name  <- as.character(feature$gene)
  #feature <- unique(feature)
  #gene.id <- data.frame(row.names = row.names(feature), gene_short_name=feature$external_gene_name, stringsAsFactors = F)
  #genes <- intersect(row.names(mat)[which(rowSums(mat) > 0)], row.names(feature))
  #mat <- mat[genes,]
  feature <- feature[row.names(mat),,drop = F]
  fd <- new("AnnotatedDataFrame", data= feature)
  pd <- new("AnnotatedDataFrame", data= meta)
  if(tpm2abs){
  cds <- newCellDataSet(as.matrix(mat), phenoData = pd, featureData = fd, lowerDetectionLimit = 0.1, expressionFamily = tobit(Lower=0.1))
  ct  <- relative2abs(cds, method = "num_genes", modelFormulaStr = '~cellType')
  lower_detect = 0.5
  min_expr = 0.1
  }else{
    ct <- as.matrix(mat)
    lower_detect = 0.5
    min_expr = 0.5
  }
  cds <- newCellDataSet(as.matrix(ct), phenoData = pd, featureData = fd, lowerDetectionLimit = lower_detect, expressionFamily = negbinomial.size())
  cds <- estimateSizeFactors(cds)
  cds <- detectGenes(cds, min_expr = min_expr)
  #cds <- cds[which(!grepl('msp-', fData(cds)$gene_short_name)),]
  ce.abs.diff <- monocle_diff(cds, subtype = c(control, target), full_formula = '~cellType+num_genes_expressed', reduced_formula =  '~num_genes_expressed', min.pct = min.pct, num_cells = 0)
  #MeanExp <- tapply(colnames(as.matrix(ct)), as.character(meta[colnames(as.matrix(ct)),]$cellType), function(x){rowMeans(as.matrix(ct)[,x])})
  #ce.abs.diff$Log2FC <- log2((MeanExp[[target]]+1e-16)/(MeanExp[[control]]+1e-16))[row.names(ce.abs.diff)]
  ce.abs.diff$Log2FC <- ce.abs.diff$Log2FC_shrunk
  ce.abs.diff$diff <- ce.abs.diff$qval < 0.05
  return(ce.abs.diff)
}




performAllDEG <- function(deg_res, obj, res_list = NULL, 
                          control = 'mouseEgg', 
                          refilter_deg = F, 
                          METHODS = c('wilcox','bimod','roc', 't', 'LR', 'MAST', 'mast', 'monocle','bpsc','scdd','basics','limma', 'scde', 'edgeR'), 
                          alpha = 0.05, 
                          FCThresh = log2(1.5),
                          min_pct = 0.5){
  other <- setdiff(unique(obj$meta$cellType), control)[1]

  if(is.null(res_list)){
    res <- list()
    genes <- list()
  }else{
    res <- res_list
    genes <- res[['gene_names']]
    res$consensus <- NULL
    res$Upregulated <- NULL
    res$tpmLog2FC <- NULL
  }
  if(refilter_deg){
       for(n in names(res)){
         if(class(res[[n]]) == 'data.frame'){
           if(grepl('mast', n)){
             degs <- subset(res[[n]], abs(Log2FC) > FCThresh & fdr < 0.05)
             genes[[n]] <- row.names(degs)
           }else if(n %in% c('bpsc', 'scdd', 'limma', 'scde', 'edgeR', 'monocle', 'monocle_ct')){
             degs <- subset(res[[n]], diff == TRUE & abs(Log2FC) > FCThresh)
             genes[[n]] <- row.names(degs)
           }else if(n == 'basics'){
             degs <- subset(res[[n]], diff == TRUE & abs(Log2FC) > FCThresh & ProbDiffMean > 0.99)
             genes[[n]] <- row.names(degs)
           }else{
             if(n == 'roc'){
               degs <- subset(res[[n]], power > 0.7 & abs(Log2FC) > FCThresh)
             }else{
               degs <- subset(res[[n]], p_val_adj < alpha & abs(Log2FC) > FCThresh)
             }
             genes[[n]] <- row.names(degs)
           }
         }else{
           next
         }
       }
  }else{
    for(method in METHODS){
    cat('performing', method, '\n')
    if(method == 'mast'){
      res[['mast_tpm']] <- mast_diff(obj, control = control, tpm = T, nbins = 10, min_per_bin = 50, freq = min_pct, max_thres = 6, plot = T)[[1]]$DESig
      degs <- subset(res[['mast_tpm']], abs(Log2FC) > FCThresh & fdr < 0.05)
      genes[['mast_tpm']] <- row.names(degs)
      
      res[['mast_cpm']] <- mast_diff(obj, control = control, tpm = F, nbins = 10, min_per_bin = 50, freq = min_pct, max_thres = 6, plot = T)[[1]]$DESig
      degs <- subset(res[['mast_cpm']], abs(Log2FC) > FCThresh & fdr < 0.05)
      genes[['mast_cpm']] <- row.names(degs)
    }else if(method == 'monocle_ct'){
      #res[['monocle_ct']] <- subset(deg_monocle(obj, control = control, tpm2abs = F), qval < 0.05 & abs(Log2FC) > FCThresh)
      res[['monocle_ct']] <- deg_monocle(obj, control = control, tpm2abs = F, min.pct = min_pct)
      degs <- subset(res[['monocle_ct']], qval < 0.05 & abs(Log2FC) > FCThresh)
      genes[[method]] <- row.names(degs)
    }
    else if(method %in% c('bpsc', 'scdd', 'limma', 'scde', 'edgeR', 'monocle')){
      #res[[method]] <- subset(get(paste('deg_',method,sep=''))(obj, control = control), diff == TRUE & abs(Log2FC) > FCThresh)
      res[[method]] <- get(paste('deg_',method,sep=''))(obj, control = control, min.pct = min_pct)
      degs <- subset(res[[method]], diff == TRUE & abs(Log2FC) > FCThresh)
      genes[[method]] <- row.names(degs)
    }
    else if(method == 'basics'){
      res[['basics']] <- deg_basics(obj,  control = control, min.pct = min_pct, spike = T)
      degs <- subset(res[['basics']], diff == TRUE & abs(Log2FC) > FCThresh & ProbDiffMean > 0.99)
      genes[[method]] <- row.names(degs)
    }
    else{
      if(method == 'roc'){
        res[[method]] <- SeuratDEG(obj, DEmethod = method, reference = control, min.pct = min_pct)
        degs <- subset(res[[method]], power > 0.7 & abs(Log2FC) > FCThresh)

      }else{
        res[[method]] <- SeuratDEG(obj, DEmethod = method, reference = control, min.pct = min_pct)
        degs <- subset(res[[method]], p_val_adj < alpha & abs(Log2FC) > FCThresh)
      }
      genes[[method]] <- row.names(degs)
      }
    }
  }
  #res[['all_sig_genes']] <- consensus_genes
  res[['gene_names']] <- genes
  consensus_genes <- consensusDEG(res)
  res[['gene_names']][['consensus']] <- row.names(consensus_genes$consensus)
  res[['gene_names']][['innerQuantile']] <- Reduce(intersect, consensus_genes$inner_degs)
  res[['consensus']] <- consensus_genes
  pos <- table(do.call(c, lapply(res, FUN = function(x){if(class(x) == 'data.frame'){row.names(subset(x, Log2FC > 0))}})))
  neg <- -table(do.call(c, lapply(res, FUN = function(x){if(class(x) == 'data.frame'){row.names(subset(x, Log2FC < 0))}})))
  posneg<- intersect(names(pos), names(neg))
  pos[posneg] <- pos[posneg]+neg[posneg]
  res[['Upregulated']] <- c(pos,neg[!names(neg) %in% posneg]) > 0
  MeanExp <- tapply(colnames(obj$tpm$bio), as.character(obj$meta[colnames(obj$tpm$bio),]$cellType), function(x){rowMeans(obj$tpm$bio[,x])})
  logfc <- log2((MeanExp[[other]]+1e-16)/(MeanExp[[control]]+1e-16))
  res[['tpmLog2FC']] <- logfc
  deg_res[[paste(control,'_v_', other, sep = '')]] <- res
  return(deg_res)
}

consensusDEG <- function(degs, pct = 0.5){
  deg_summary <- summary(sapply(degs[['gene_names']], length))
  degs_inner <- lapply(names(degs[['gene_names']]), FUN = function(x){
      if(length(degs[['gene_names']][[x]]) > deg_summary['1st Qu.'] & length(degs[['gene_names']][[x]]) < deg_summary['3rd Qu.']){
        return(degs[['gene_names']][[x]])
      }
        }
      )
  names(degs_inner) <-names(degs[['gene_names']])
  for(n in names(degs_inner)){
    if(is.null(degs_inner[[n]])){
      degs_inner[[n]] <- NULL
    }
  }
  deg_ct <- do.call(rbind, lapply(names(degs), FUN = function(x){
    deg_df <- degs[[x]]
    if(class(deg_df) == 'data.frame'){
      genes <- degs[['gene_names']][[x]]
      return(data.frame(gene_names = genes, Log2FC = deg_df[genes,]$Log2FC, N = rep(1, length(genes))))
    }else{return(NULL)}
    }))
  #deg_ctfc <- data.frame(aggregate(Log2FC~gene_names, deg_ct, function(x){return(x[which(abs(x) == min(abs(x)))])}), row.names = 1)
  deg_ctfc <- data.frame(aggregate(Log2FC~gene_names, deg_ct, median), row.names = 1)
  deg_ctct <- data.frame(aggregate(N~gene_names, deg_ct, sum), row.names = 1)
  deg_inner_ct <- table(do.call(c, lapply(degs_inner, FUN = function(x){return(x)})))
  consensus_deg <- data.frame(row.names = row.names(deg_ctfc), Log2FC = deg_ctfc$Log2FC, N = deg_ctct[row.names(deg_ctfc),] )
  consensus_deg[abs(consensus_deg$Log2FC) > 10,]$Log2FC <- sign(consensus_deg[abs(consensus_deg$Log2FC) > 10,]$Log2FC)*median(abs(consensus_deg$Log2FC))
  consensus_deg_inner <- deg_inner_ct[deg_inner_ct >= (length(degs_inner)-1)]
  return(list(consensus = subset(consensus_deg, N > 5), consensus_all = consensus_deg, inner_quantile = consensus_deg_inner, inner_degs = degs_inner))
}


consensusPlot <- function(css, file='./deg/', title = ''){
  con <- css
  con$State <- con$N >= 5
  con[con$N >= 5,]$State <- 'DE'
  con[con$N < 5,]$State <- 'No-DE'
  con$N <- as.factor(con$N)
  p<-ggplot(con, aes(x=N, y=Log2FC, fill=State)) +
    geom_violin(trim=FALSE)+ylab(expression('-Log'['2']*'FC'))+xlab('N methods')+ggtitle(title)+theme(
      plot.title = element_text(hjust = 0.5, size = 14))
  p <- p + scale_fill_grey() + theme_classic()
  print(p)
  #ggsave(file, plot = p, width = 5, height = 3)
}




dtu_drimseq <- function(salmon_dir, meta, t2g){
  tx2gene <- read.csv(t2g)
  read_expression(salmon_dir, mode = 'salmon', tx2gene)
}


homologyCompare <- function(degMouse, degRat)
{
  require(clusterProfiler)
  require(org.Mm.eg.db)
  require(org.Rn.eg.db)
  mS2E <- bitr(row.names(degMouse), fromType = 'ALIAS', 'ENSEMBL', org.Mm.eg.db)
  mS2E <- data.frame(subset(mS2E[match(row.names(degMouse), mS2E$ALIAS),], !is.na(ALIAS)), row.names = 1)
  rS2E <- bitr(row.names(degRat), fromType = 'ALIAS', 'ENSEMBL', org.Rn.eg.db)
  rS2E <- data.frame(subset(rS2E[match(row.names(degRat), rS2E$ALIAS),], !is.na(ALIAS)), row.names = 1)
  degMouse <- degMouse[intersect(row.names(degMouse), row.names(mS2E)),]
  degRat <- degRat[intersect(row.names(degRat), row.names(rS2E)),]
  degM <- mS2E[row.names(degMouse),1 ]
  degMouseUp <- mS2E[row.names(subset(degMouse, Log2FC > 0)),1 ]
  degMouseDown <- mS2E[row.names(subset(degMouse, Log2FC < 0)),1 ]
  degR<- rS2E[row.names(degRat),1 ]
  degRatUp <- rS2E[row.names(subset(degRat, Log2FC > 0)), 1]
  degRatDown <- rS2E[row.names(subset(degRat, Log2FC < 0)),1 ]
  rOrtho <- intersect(degR, one2one$Rat.gene.stable.ID)
  mOrtho <- intersect(degM, one2one$Gene.stable.ID)
  #mOrthoSymbol <- data.frame(mS2E[match(mOrtho, mS2E$ENSEMBL),1] , row.names = mS2E$ENSEMBL[match(mOrtho, mS2E$ENSEMBL)])[mOrtho, 1, drop =F]
  #rOrthoSymbol <- data.frame(rS2E[match(rOrtho, rS2E$ENSEMBL),1], row.names = rS2E$ENSEMBL[match(rOrtho, rS2E$ENSEMBL)])[data.frame(one2one, row.names = 1)[mOrtho]$Rat.gene.stable.ID,1, drop =F]
  rOrtho <- subset(one2one, Rat.gene.stable.ID %in% rOrtho)$Gene.stable.ID
  
  rOrthoU <- intersect(degRatUp, one2one$Rat.gene.stable.ID)
  mOrthoU <- intersect(degMouseUp, one2one$Gene.stable.ID)
  rOrthoU <- subset(one2one, Rat.gene.stable.ID %in% rOrthoU)$Gene.stable.ID
  
  rOrthoD <- intersect(degRatDown, one2one$Rat.gene.stable.ID)
  mOrthoD <- intersect(degMouseDown, one2one$Gene.stable.ID)
  rOrthoD <- subset(one2one, Rat.gene.stable.ID %in% rOrthoD)$Gene.stable.ID
  
  
  return(list(mouseHomologDEG = mOrtho, mouseHomologDEGup = mOrthoU, mouseHomologDEGdown = mOrthoD, 
              ratHomologDEG = rOrtho, ratHomologDEGup = rOrthoU, ratHomologDEGdown = rOrthoD, mS2E =  mS2E, 
         rS2E =  rS2E))
}


deg_deseq <- function(obj = NULL, ct = NULL, meta = NULL, control = 'control', lfc_thresh = 1, deFactor = 'condition',min_cell_grp = 2, min_cell = 0, min.pct = 0, min.mean.expr = 5){
    
    if(!is.null(obj)){
      mat <- round(obj$ct$bio)
      mat <- apply(mat,2,function(x) {storage.mode(x) <- 'integer'; x})
      meta <- obj$meta
      meta$condition <- as.factor(meta[,deFactor])
    }else{
      mat <- round(ct)
      mat <- apply(mat,2,function(x) {storage.mode(x) <- 'integer'; x})
      meta <- meta
      meta$condition <- as.factor(meta[,deFactor])
    }
    target <- setdiff(unique(meta$condition), c(control))
    ctrl_mat <- mat[,row.names(subset(meta, condition == control))]
    ctrl_mat <- ctrl_mat[rowSums(ctrl_mat != 0) > max(ncol(ctrl_mat)*min.pct, min_cell_grp),]
    trgt_mat <- mat[,row.names(subset(meta, condition == target))]
    trgt_mat <- trgt_mat[rowSums(trgt_mat != 0) > max(ncol(trgt_mat)*min.pct, min_cell_grp),]
    genes <- union(row.names(ctrl_mat), row.names(trgt_mat))
    print(length(genes))
    genes <- intersect(row.names(mat)[which(rowSums(mat > 0) > min_cell & rowMeans(mat) > min.mean.expr)], genes)
    print(length(genes))
    mat <- mat[genes,c(row.names(subset(meta, condition == control)), row.names(subset(meta, condition == target)))]
    dds <- DESeqDataSetFromMatrix(countData = mat,
                                        colData = meta,
                                        design = ~ condition)
    dds$condition <- relevel(dds$condition, ref = control)
    dds <- DESeq(dds)
    res_ph <- DESeq2::results(dds, alpha = 0.05)
    res_thresh_1 <- DESeq2::results(dds, alpha=0.05,lfcThreshold=1, altHypothesis="greaterAbs")
    res_thresh_custom <- DESeq2::results(dds, alpha=0.05,lfcThreshold=lfc_thresh, altHypothesis="greaterAbs")
    reslfc <- lfcShrink(dds, coef=paste("condition",target,"vs",control, sep ='_'), type="apeglm")
    reslfc_thresh <- lfcShrink(dds, coef="condition_treated_vs_control", type="apeglm", lfcThreshold = lfc_thresh)
    
    return(list(dds = dds, ph = res_ph, thresh = res_thresh_custom, thresh_1 = res_thresh_1, lfc = reslfc, lfc_thresh = reslfc_thresh))
}


deg_utr <- function(file, ct, compare, meta, impute = F, method = 'fisher.test', combine_p = 'fisher'){
  ##read in the new dapars file that includes long. short and PDUI values
  dapars <- read.csv(file, sep = '\t', header = T, row.names = 1)
  print(dim(dapars))
  ## Filter first based on coverage of the each gene's entire body
  ## Only account for UTR coverage in genes that are assigned at have mean expression of at least 10 uniquely mapped reads across all samples
  genes <- row.names(ct)[rowMeans(ct) > 10]
  #genes <- intersect(genes, row.names(ct)[rowSums(ct[,row.names(subset(meta, cellType == compare[2]))] > 10) > 0])
  dapars$gene_short_names <- sapply(row.names(dapars), FUN = function(x){strsplit(x,"\\|")[[1]][2]})
  print(length(unique((dapars$gene_short_names))))
  
  dapars_orig <- dapars
  all_genes <- unique(dapars$gene_short_names)
  dapars <- dapars[dapars$gene_short_names %in% genes,]
  dapars <- subset(dapars, fit_value >= 10) # Maybe filter also based on regression fit value
  print(length(unique((dapars$gene_short_names))))
  # vector to store gene name and UTR region correspondence in case gene names is lost with imputation
  gene2region <- dapars$gene_short_names
  names(gene2region) <- sapply(row.names(dapars), FUN = function(x){strsplit(x,"\\|")[[1]][1]})
  
  #split file into long, short and pdui
  d_long <- dapars[,grepl('long_exp', colnames(dapars))]
  d_short <- dapars[,grepl('short_exp', colnames(dapars))]
  d_pdui <- dapars[,grepl('PDUI', colnames(dapars))]
  
  # change column names
  colnames(d_long) <- sapply(strsplit(colnames(d_long), '_'), FUN = function(x){strsplit(x[1], "\\.")[[1]][1]})
  colnames(d_short) <- sapply(strsplit(colnames(d_short), '_'), FUN = function(x){strsplit(x[1], "\\.")[[1]][1]})
  colnames(d_pdui) <- sapply(strsplit(colnames(d_pdui), '_'), FUN = function(x){strsplit(x[1], "\\.")[[1]][1]})
  # select filtering genes based on number of passes (Non NAs) and overall coverage (average > 2 either in long or short UTR in both conditions)
  grp = meta[colnames(d_pdui),'condition']
  cond1_ind = which(grp == compare[1])
  cond2_ind = which(grp == compare[2])
  #print(d_pdui[dapars$gene_short_name == 'Cdk1',])
  ## NA FILTER
  # First filter out all genes that have at more than 1 NAs in terms of coverage in both conditions
  na.filt.genes <- rowSums(!is.na(d_pdui[,cond1_ind])) >= 3 & rowSums(!is.na(d_pdui[,cond2_ind])) >= 3
  dapars <- dapars[na.filt.genes, ]
  
  print(length(na.filt.genes))
  print(sum(na.filt.genes))

  d_long <- d_long[na.filt.genes, ]
  d_short <- d_short[na.filt.genes, ]
  d_pdui <- d_pdui[na.filt.genes, ]
  
  print(dim(d_long))
  print(dim(d_short))
  print(dim(d_pdui))
  # Filter based on coverage of UTR regions
  c1.filt.genes <- rowMeans(d_long[,cond1_ind], na.rm = T) > 1 | rowMeans(d_short[,cond1_ind], na.rm = T) > 1
  print(sum(c1.filt.genes))
  c2.filt.genes <- rowMeans(d_long[,cond2_ind], na.rm = T) > 1 | rowMeans(d_short[,cond2_ind], na.rm = T) > 1
  print(sum(c2.filt.genes))
  print(sum(c2.filt.genes))
  final.filt.genes <- c1.filt.genes & c2.filt.genes 
  #print(sum(final.filt.genes))
  # filtering matrices with genes selected prior
  dapars <- dapars[final.filt.genes, ] 
  d_long <- d_long[final.filt.genes, ] %>%  mutate(
    across(everything(), replace_na, 0)
  )
  d_short <- d_short[final.filt.genes, ] %>%  mutate(
    across(everything(), replace_na, 0)
  )
  d_pdui <- d_pdui[final.filt.genes, ] 
  print(length(unique((dapars$gene_short_names))))
  #print('Cdk1' %in% dapars$gene_short_name)
  if(impute){
  dapars_out <- data.frame(Gene = row.names(dapars), dapars[1:3,], d_pdui)
  write.table(dapars_out, file = './temp.dp.tsv', sep = '\t', quote = F, row.names = F)
  d_pdui = scDaPars(raw_PDUI_file = './temp.dp.tsv',
                           out_dir = "apa/scDaPars_result",
                           filter_gene_thre = 0.2,
                           filter_cell_thre = 0.1, k= 8)
  method = 'ks.test'
  }
  #dapars_ratio <- dapars_ratio[rowSums(!is.na(dapars_ratio)) > 20 ,]
  
    #dapars_ratio[is.na(dapars_ratio)] <- 0

  #print(sum(rowSums(!is.na(dapars_ratio[,cond1_ind])) >= 10 & rowSums(!is.na(dapars_ratio[,cond2_ind])) >=10))
  #print(sum(rowSums(!is.na(dapars_ratio)) > 20))
  #dapars_ratio <- dapars_ratio[rowSums(!is.na(dapars_ratio[,cond1_ind])) >= 7 & rowSums(!is.na(dapars_ratio[,cond2_ind])) >=7,]
  cat('starting now\n')
  if(method == 'ks.test'){
    test <- apply(d_pdui, 1, FUN = function(x){
      if(sum(x[!is.na(x)]) == 0 | sum(!is.na(x[cond1_ind])) <= 3 | sum(!is.na(x[cond2_ind])) <= 3){
       c(1, 0)
      }else{
      c(ks.test(x[cond1_ind][!is.na(x[cond1_ind])], x[cond2_ind][!is.na(x[cond2_ind])])$p.value, mean(x[cond2_ind][!is.na(x[cond2_ind])]) - mean(x[cond1_ind][!is.na(x[cond1_ind])]))
      }
    })
    test <- t(test)
  }else{
    l1 = length(cond1_ind)
    l2 = length(cond2_ind)
    utrl1_mean <- rowMeans(d_long[, cond1_ind], na.rm = T)
    utrl2_mean <- rowMeans(d_long[, cond2_ind], na.rm = T)
    utrs1_mean <- rowMeans(d_short[, cond1_ind], na.rm = T)
    utrs2_mean <- rowMeans(d_short[, cond2_ind], na.rm = T)
      test <- do.call(rbind, pblapply(row.names(d_long), FUN = function(x){
      #utr_l1 <- d_long[x,cond1_ind][!is.na(d_long[x,cond1_ind])] # long utr coverage in condition 1
      #utr_l2 <- d_long[x,cond2_ind][!is.na(d_long[x,cond2_ind])] # short utr coverage in condition 2
      #utr_s1 <- d_short[x,cond1_ind][!is.na(d_short[x,cond1_ind])] # long utr coverage in condition 1
      #utr_s2 <- d_short[x,cond2_ind][!is.na(d_short[x,cond2_ind])] # short utr coverage in condition 2
      pdui_1 <- mean(d_pdui[x,cond1_ind][!is.na(d_pdui[x,cond1_ind])])
      pdui_2 <- mean(d_pdui[x,cond2_ind][!is.na(d_pdui[x,cond2_ind])])
      #c(fisher.test(x = rbind(c(mean(utr_l1), mean(utr_s1)), c(mean(utr_l2), mean(utr_s2))))$p.value, pdui_2 - pdui_1)
      twobytwo <- rbind(c(utrl1_mean[x], utrs1_mean[x]), c(utrl2_mean[x], utrs2_mean[x]))
      c(fisher.test(x = twobytwo)$p.value, pdui_2 - pdui_1)
      
      
      #if(x == 'XM_039113041.1|Pou2f2|NC_051336.1|-'){
       # print(rbind(c(sum(utr_l1)/l1, sum(utr_s1)/l1), c(sum(utr_l2)/l1, sum(utr_s2)/l2)))
      #}
            }))
  }
  cat('finished\n')
  #print(dim(test))
  #print(dim(d_long))
  row.names(test) <- row.names(d_long)
  colnames(test) <- c('pval', 'mean.diff')
  test <- data.frame(test)
  test$pval[test$pval > 1] = 1
  test$padj <- p.adjust(test$pval)
  test$fdr <- qvalue(test$pval)$qvalue
  
  if(impute){
    test$gene_short_names <- gene2region[row.names(test)]
  }else{
    test$gene_short_names <- sapply(row.names(test), FUN = function(x){strsplit(x,"\\|")[[1]][2]})
  }
  test$diff <- abs(test$mean.diff) > 0.2 & test$fdr < 0.05
  test$fit_value <- dapars[row.names(test),]$fit_value
  test$predicted_p_APA <- dapars_orig[row.names(test),]$Predicted_Proximal_APA
  test$loci <- dapars_orig[row.names(test),]$Loci
  test$strand = sapply(strsplit(row.names(test), '\\|'), FUN = function(x){x[4]})
  test$APA_dist = 0
  test[test$strand == '+',]$APA_dist <- abs(sapply(strsplit(test[test$strand == '+',]$loci, '-'), 
                                                  FUN = function(x){as.numeric(strsplit(x[1], ':')[[1]][2])}) - test[test$strand == '+',]$predicted_p_APA)-1
  test[test$strand == '-',]$APA_dist <- abs(sapply(strsplit(test[test$strand == '-',]$loci, '-'), 
                                                  FUN = function(x){as.numeric(x[2])}) - test[test$strand == '-',]$predicted_p_APA)-1
  print('done')
  #test$APA_dist <- abs(sapply(strsplit(test$loci, '-'), FUN = function(x){as.numeric(x[2])}) - test$predicted_p_APA)-1
  
  gene_res <- data.frame(do.call(rbind, tapply(row.names(test), test$gene_short_names, function(x){
    df <- test[x,]
    min_pval <- min(df[,'pval'])
    min_pval_ind = which(df[,'pval'] == min_pval)
    min_pval_ind <- min_pval_ind[which(abs(df[min_pval_ind, 'mean.diff']) == max(abs(df[min_pval_ind, 'mean.diff'])))][1]
    if(combine_p == 'fisher'){
      df[min_pval_ind,'pval'] <- metapod::combineParallelPValues(as.list(df[,'pval']), method = 'fisher')$p.value
    }
    return(cbind(df[min_pval_ind, ], dapars[x[min_pval_ind],c(1,2,3)]))
  })))
  gene_res$padj <- p.adjust(gene_res$pval)
  gene_res$fdr <- qvalue(gene_res$pval)$qvalue
  gene_res$diff <- abs(gene_res$mean.diff) > 0.2 & gene_res$fdr < 0.05
  print(length(unique((test$gene_short_names))))
  pdui <- dapars_orig[,grepl('PDUI', colnames(dapars_orig))]
  colnames(pdui) <- colnames(d_long)
  pdui <- pdui[row.names(d_long),]
  pdui_impute <- t(apply(pdui, 1, FUN = function(x){x[is.na(x)] = mean(x, na.rm =T); x}))
  
  return(list(deg= test, long = d_long, gene_res = gene_res, short = d_short, df = dapars_orig, pdui = pdui, pdui_imp = pdui_impute, gene_universe = all_genes))
}




test_ortho <- do.call(rbind, tapply(orthology[,c(1,3,4)], INDEX = orthology$DB.Class.Key, FUN = function(x){
  if(nrow(x) == 1 | length(unique(x[,2])) == 1 | length(unique(x[,2])) > 2){
    return(NULL)
  }else{
    m = unique(x[,2])[1]
    n = unique(x[,2])[2]
    db_class = unique(x[,1])
    new_ = cross_join(data.frame(x[x[,2] == m,3, drop = F]), data.frame(x[x[,2] == n,3, drop = F]))
    
    new_ <- cbind(rep(db_class, nrow(new_)), new_)
    colnames(new_) = c('db_class',m,n)
    new_
  }
}))
colnames(test_ortho) <- c('ortho_class', 'mouse', 'rat')

## Function to combine expression matrices between 2 species based on orthologous genes
## Orthologous genes must be in a table of 3 columns: ortho_class, spieces1, spieces2
## Default is use only one2one mappings of orthologous genes
combine2SpieciesOrthoExpression <- function(species1, species2, ortho_table, one2many = F, many2many = F){
  o2o_class <- names(table(ortho_table[,'ortho_class'])[table(ortho_table[,'ortho_class']) == 1])
  o2o_tbl <- ortho_table[match(o2o_class, ortho_table[,'ortho_class'] ),]
  combined <- do.call(rbind, lapply(o2o_tbl[,'ortho_class'], FUN =function(n){
    x = o2o_tbl[o2o_tbl[,'ortho_class'] == n,]
    if(x[,2] %in% row.names(species1)){
      sp1 = species1[x[,2],]
      
    }else{
      sp1 = rep(0, ncol(species1))
      #return(as.numeric(rep(0, ncol(species1)+ncol(species2))))
    }
    if(x[,3] %in% row.names(species2)){
      sp2 = species2[x[,3],]
      
    }else{
      sp2 = rep(0, ncol(species2))
      #return(as.numeric(rep(0, ncol(species1)+ncol(species2))))
    }
    new = c(as.numeric(sp1), as.numeric(sp2))
    return(new)
  }))
  colnames(combined) <- c(colnames(species1), colnames(species2))
  row.names(combined) <- o2o_tbl[,'ortho_class']
  combined <- combined[rowSums(combined) != 0,]
  combined
}


read_filter_splice <- function(splice_f, unsplic_f){
  
  spliced <- t(read.csv(splice_f, header = T, row.names = 1))
  spliced_raw <- t(read.csv(splice_f, header = T, row.names = 1))
  unsplic <- t(read.csv(unsplic_f, header = T, row.names = 1))
  unsplic_raw <- t(read.csv(unsplic_f, header = T, row.names = 1))
  grp1 <- grepl(pattern = 'F', colnames(spliced))
  grp2 <- grepl(pattern = 'U', colnames(spliced))
  grp1_filt_s <- rowMeans(spliced[,grp1]) > 5 & rowSums(spliced[,grp1] > 0) > sum(grp1)*0.1
  grp1_filt_u <- rowMeans(unsplic[,grp1]) > 1 & rowSums(unsplic[,grp1] > 0) > sum(grp1)*0.1
  
  grp2_filt_s <- rowMeans(spliced[,grp2]) > 5 & rowSums(spliced[,grp2] > 0) > sum(grp2)*0.1
  grp2_filt_u <- rowMeans(unsplic[,grp2]) > 1 & rowSums(unsplic[,grp2] > 0) > sum(grp2)*0.1
  
  spliced <- spliced[(grp1_filt_s &  grp1_filt_u) | (grp2_filt_s & grp2_filt_u),]
  unsplic <- unsplic[(grp1_filt_s &  grp1_filt_u) | (grp2_filt_s & grp2_filt_u),]
  
  
  test_cor <- do.call(rbind, lapply(1:nrow(spliced), FUN = function(x){
    c('pc' = cor(spliced[x,], unsplic[x,]),
      'sc' = cor(spliced[x,], unsplic[x,], method = 'spearman'),
      'kc' = cor(spliced[x,], unsplic[x,], method = 'kendall'),
      'rsq' = summary(lm(spliced[x,]~unsplic[x,]))$adj.r.squared,
      'err.rat' = sd(unsplic[x,])/sd(spliced[x,]),
      'expr.rat' = mean(unsplic[x,])/mean(spliced[x,])
    )}
  )
  )
  row.names(test_cor) <- row.names(spliced)
  test_cor <- data.frame(test_cor)
  genes <- row.names(subset(test_cor, rsq > 0.1 & kc > 0.1 & err.rat > 0.005 & err.rat < 5))
  return(list(spliced = spliced_raw, unspliced=unsplic_raw, metrics = test_cor, genes = genes))
}



fisher_proportion_test <- function(nascent, mature, groups, control = 'mouseEgg'){
  grp1 <- which(groups == control)
  grp2 = which(groups == setdiff(groups, control))
  nascent_c <- nascent[,grp1]
  nascent_t <- nascent[,grp2]
  mature_c <- mature[,grp1]
  mature_t <- mature[,grp2]
  nascent_p1 <- rowMeans(nascent_c/(mature_c+nascent_c), na.rm = T)
  nascent_p01 <- rowSums(nascent_c)/rowSums(mature_c)
  
  nascent_p2 <- rowMeans(nascent_t/(mature_t+nascent_t), na.rm = T)
  nascent_p02 <- rowSums(nascent_t)/rowSums(mature_t)
  res <- data.frame(do.call(rbind, lapply(1:nrow(nascent), FUN = function(i){
    #diff <- c(nascent_p2[i]- nascent_p1[i])/mean(nascent[i,]/c(mature[i,]+nascent[i,]), na.rm=T)
    diff <- log2(nascent_p02[i]/nascent_p01[i])#/(sum(nascent[i,])/sum(mature[i,]))
    c(fisher.test(x = rbind(c(mean(nascent_c[i,]), mean(mature_c[i,])), c(mean(nascent_t[i,]), mean(mature_t[i,]))))$p.value, diff)
  })))
  
  row.names(res) <- row.names(nascent)
  colnames(res) <- c('fisher.p', 'prop_diff')
  res[res$fisher.p >= 1,1] <- 1
  res$qvalue <- qvalue(res[,1])$qvalue
  res$features <- row.names(res)
  res
}
