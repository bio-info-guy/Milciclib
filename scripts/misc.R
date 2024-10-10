#Misc Functions
##GENERAL FUNCTIONS
# PCA
# HEATMAP
# SPIKEIN-BASED NORMALIZATION (a version of TPM conversion)


#PLOTTING SPIKEIN Linear Models
spikein_plot <- function(tpm, true, lm, xlab = "spikein tpm", ylab = "true spikein counts" , main = "spike-in plot")
{
  plot( NULL, xaxt="n", yaxt="n",
        log="xy", xlim = c( 0.01, 1.5*10^5 ), ylim = c( .05, 2*10^6 ),
        xlab = xlab, ylab = ylab , main = main)
  axis( 1, 2*10^(-1:5), c( "0.1", "1", "10", "100", "1000",
                           expression(10^4), expression(10^5) ) )
  axis( 2, 2*10^(-2:6), c( "0.01", "0.1", "1", "10" , "100", "1000",expression(10^4), expression(10^5), expression(10^6)), las=2 )
  abline( h=2*10^(-2:6), v=2*10^(-1:5), col="#D0D0D0", lwd=2 )
  points( rowMeans(tpm), true, pch=20, cex=1, col="#00207040" )
  xg <- 10^seq( -2, 6, length.out=1000 )
  lines(xg, 2^(coef(lm)[1]*log2(xg+1))-1)
}

#SPIKEIN NORMALIZATION WITH TPM VS TRUE COUNTS
spikein_normalization <- function(tpm, sptpm, meta, 
                                  sp.true="~/Desktop/mouse_rat/spikeins/noise_correction/ERCC_conc.tsv", 
                                  sizeFactor = "ERCC", 
                                  plot=TRUE, 
                                  deseq=TRUE, 
                                  sp_batch=FALSE, 
                                  ct_batch=sp_batch, 
                                  outlier="ERCC-00074", 
                                  return_ct=FALSE)
{
  if (sp_batch != FALSE)
  {
    print(sp_batch)
    meta <- meta[which(meta$batch == sp_batch),]
    spct <- spct[,row.names(meta)]
  }
  if (ct_batch != FALSE)
  {
    print(ct_batch)
    ct <- ct[, row.names(meta[which(meta$batch == ct_batch),])]
  }
  #read spike in length
  spikein.length <- read.table(splen, row.names=1)#take in spikein length
  spikein.length <- spikein.length[row.names(spct),]#reorder the spikein.length genes
  spikein.length <- data.frame(row.names = rownames(spikein.length),
                               length=spikein.length$V5, 
                               type=rep("spike", nrow(spikein.length)))#build data.frame of gene length
  #return spike in sizeFactors either from DESeq or library size
  sf <- tpm_norm(spct, spikein.length$length, 
                 deseq = deseq, return_sf = TRUE)#tpm normalization of spikein genes (deseq normalization factor)
  spikein.ct <- spct[,which(meta$spikein.rate > 0.001)]
  
  spike.norm <- tpm_norm(spikein.ct, spikein.length$length, deseq = deseq)
  
  spike.norm <- spike.norm[rowSums(spike.norm == 0)<=ncol(spike.norm)*0.6, ]
  spike.norm <- spike.norm[!rownames(spike.norm) %in% outlier,]
  spikein.true <- read.table(sp.true, header=TRUE, row.names = 1)
  spikein.true <- spikein.true[row.names(spike.norm),]
  sc.lm <- lm(log2(spikein.true+1)~rowMeans(log2(spike.norm+1))-1)
  print(summary(sc.lm))
  if ( plot == TRUE)
  {
    spikein_plot(spike.norm, spikein.true, sc.lm, main="Yeast ERCC Spike-in plot ")
  }
  if (return_ct == TRUE)
  {
    if(sizeFactor == "ERCC")
    {norm <- tpm_norm(ct, bio.len$effective.length, deseq = deseq, sf = sf)}
    else
    {norm <- tpm_norm(ct, bio.len$effective.length, deseq = deseq)}
    tc <- 2^(coef(sc.lm)[1]*log2(norm+1))-1
    return(tc)
  }
}

write_out_enrichment <- function(res.env, file="enrichment_result.tsv")
{
  env <- res.env
  for (item in ls(env))
  {
    if(!grepl('_v_',  item)){
      next
    }
    
    df <- env[[item]]$all
    df <- df[order(df$padj),]
    df$pval <- formatC(df$pval, format = "e", digits = 2)
    df$padj <- formatC(df$padj, format = "e", digits = 2)
    df$score <- round(df$score, 2)
    #View(df)
    #write(paste("\n",item, sep=""), file=file,  append = TRUE)
    write.table(df, file = paste(file,'_',item,'.tsv', sep=''),  quote = FALSE, sep = "\t",  row.names = F)
  }
}

##PCA and Heatmaps

#getting a dataframe of biocyc pathways to genes

biocyc <- yeast$biocyc_gene2path
kegg <- unique(data.frame(gene=as.character(yeast$features[yeast$features$kegg_enzyme != "",]$ensembl_gene_id), pathway=as.character(yeast$features[yeast$features$kegg_enzyme != "",]$kegg_enzyme)))
obsClustMatrix <- function(df, genes)
{
  sameCluster <- function(i,j,n, df)
  {
    return(length(intersect(df[df$gene == n[i],]$pathway, df[df$gene == n[j],]$pathway)) > 0)
  }
  y <- matrix(data=0, nrow=(length(genes)), ncol=(length(genes)))
  for (i in 1:length(genes))
  {
    for (j in i:length(genes))
    {
      if(sameCluster(i,j,genes, df))
      {
        y[i,j] <- 1
        y[j,i] <- 1
      }
    }
  }
  return(y)
}

bicycMat <- obsClustMatrix(biocyc)
row.names(biocycMat) <- unique(yeast$biocyc_gene2path$gene)
colnames(biocycMat) <- unique(yeast$biocyc_gene2path$gene)

keggMat <- obsClustMatrix(kegg, unique(kegg$gene))
row.names(keggMat) <- unique(kegg$gene)
colnames(keggMat) <- unique(kegg$gene)






###
#PWY-5079 aa_v_hyp
#PWY-7118 glu v aa
#PWY-7118 PWY-5076  glu v hyp
# ANAGLYCOLYSIS-PWY glu v iso
# PWY-821 iso v aa
#GLYCOLYSIS-YEAST-PWY iso v hyp

normalize_expr_data <- function(cds,
                                norm_method = c("log", "vstExprs", "none"),
                                pseudo_expr = 1,
                                relative_expr = TRUE){
  FM <- exprs(cds)
  use_for_ordering <- NULL
  # If the user has selected a subset of genes for use in ordering the cells
  # via setOrderingFilter(), subset the expression matrix.
  if (is.null(fData(cds)$use_for_ordering) == FALSE &&
      nrow(subset(fData(cds), use_for_ordering == TRUE)) > 0) {
    FM <- FM[fData(cds)$use_for_ordering, ]
  }
  
  norm_method <- match.arg(norm_method)
  if (cds@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
    
    # If we're going to be using log, and the user hasn't given us a pseudocount
    # set it to 1 by default.
    if (is.null(pseudo_expr)){
      if(norm_method == "log")
        pseudo_expr = 1
      else
        pseudo_expr = 0
    }
    
    #checkSizeFactors(cds)
    
    if (norm_method == "vstExprs") {
      if (relative_expr == FALSE)
        message("Warning: relative_expr is ignored when using norm_method == 'vstExprs'")
      
      if (is.null(fData(cds)$use_for_ordering) == FALSE &&
          nrow(subset(fData(cds), use_for_ordering == TRUE)) > 0) {
        VST_FM <- vstExprs(cds[fData(cds)$use_for_ordering,], round_vals = FALSE)
      }else{
        VST_FM <- vstExprs(cds, round_vals = FALSE)
      }
      
      if (is.null(VST_FM) == FALSE) {
        FM <- VST_FM
      }
      else {
        stop("Error: set the variance-stabilized value matrix with vstExprs(cds) <- computeVarianceStabilizedValues() before calling this function with use_vst=TRUE")
      }
    }else if (norm_method == "log") {
      # If we are using log, normalize by size factor before log-transforming
      
      if (relative_expr)
        FM <- Matrix::t(Matrix::t(FM)/sizeFactors(cds))
      
      if(is.null(pseudo_expr))
        pseudo_expr <- 1 
      FM <- FM + pseudo_expr
      FM <- log2(FM)
    }else if (norm_method == "none"){
      # If we are using log, normalize by size factor before log-transforming
      FM <- Matrix::t(Matrix::t(FM)/sizeFactors(cds))
      FM <- FM + pseudo_expr
    }
  }else if (cds@expressionFamily@vfamily == "binomialff") {
    if (norm_method == "none"){
      #If this is binomial data, transform expression values into TF-IDF scores.
      ncounts <- FM > 0
      ncounts[ncounts != 0] <- 1
      FM <- Matrix::t(Matrix::t(ncounts) * log(1 + ncol(ncounts)/rowSums(ncounts)))
    }else{
      stop("Error: the only normalization method supported with binomial data is 'none'")
    }
  }else if (cds@expressionFamily@vfamily == "Tobit") {
    FM <- FM + pseudo_expr
    if (norm_method == "none"){
      
    }else if (norm_method == "log"){
      FM <- log2(FM)
    }else{
      stop("Error: the only normalization methods supported with Tobit-distributed (e.g. FPKM/TPM) data are 'log' (recommended) or 'none'")
    }
  }else if (cds@expressionFamily@vfamily == "uninormal") {
    if (norm_method == "none"){
      FM <- FM + pseudo_expr
    }else{
      stop("Error: the only normalization method supported with gaussian data is 'none'")
    }
  }
  # if(norm_method != "none")
  #normalize_expr_data
  return (FM)
}


BrenneckeGetVariableGenes <- function(expr_mat, spikes=NA, suppress.plot=FALSE, fdr=0.1, minBiolDisp=0.5, fitMeanQuantile=0.8) {
  #require(statmod)
  
  rowVars <- function(x) { unlist(apply(x,1,var, na.rm=TRUE))}
  
  colGenes <- "black"
  colSp <- "blue"
  
  
  fullCountTable <- expr_mat;
  
  if (is.character(spikes)) {
    sp <- rownames(fullCountTable) %in% spikes;
    countsSp <- fullCountTable[sp,];
    countsGenes <- fullCountTable[!sp,];
  } else if (is.numeric(spikes)) {
    countsSp <- fullCountTable[spikes,];
    countsGenes <- fullCountTable[-spikes,];
  } else {
    countsSp <- fullCountTable;
    countsGenes <- fullCountTable;
  }
  
  meansSp <- rowMeans(countsSp, na.rm=TRUE)
  varsSp <- rowVars(countsSp)
  cv2Sp <- varsSp/meansSp^2
  View(data.frame(meansSp, varsSp))
  meansGenes <- rowMeans(countsGenes, na.rm=TRUE)
  varsGenes <- rowVars(countsGenes)
  cv2Genes <- varsGenes/meansGenes^2
  # Fit Model
  print(meansSp["YMR072W"])
  minMeanForFit <- unname( quantile( meansSp[ which( cv2Sp > 0.3 ) ], fitMeanQuantile))
  print(minMeanForFit)
  useForFit <- meansSp >= minMeanForFit
  if (sum(useForFit, na.rm=TRUE) < 20) {
    warning("Too few spike-ins exceed minMeanForFit, recomputing using all genes.")
    meansAll <- c(meansGenes, meansSp)
    cv2All <- c(cv2Genes,cv2Sp)
    minMeanForFit <- unname( quantile( meansAll[ which( cv2All > 0.3 ) ], 0.80))
    useForFit <- meansSp >= minMeanForFit
  }
  if (sum(useForFit, na.rm=TRUE) < 30) {warning(paste("Only", sum(useForFit), "spike-ins to be used in fitting, may result in poor fit."))}
  fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansSp[useForFit] ), cv2Sp[useForFit] )
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"])
  print(a0); print(a1)
  res <- cv2Genes - (a0 + a1/meansGenes);
  
  # Test
  psia1theta <- a1
  minBiolDisp <- minBiolDisp^2
  m <- ncol(countsSp);
  cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
  testDenom <- (meansGenes*psia1theta + meansGenes^2*cv2th)/(1+cv2th/m)
  p <- 1-pchisq(varsGenes * (m-1)/testDenom,m-1)
  padj <- p.adjust(p,"BH")
  sig <- padj < fdr
  sig[is.na(sig)] <- FALSE
  if (!suppress.plot) {
    plot( meansGenes,cv2Genes, xaxt="n", yaxt="n", log="xy",
          xlab = "average normalized read count",
          ylab = "squared coefficient of variation (CV^2)", col="white")
    axis( 1, 10^(-2:5), c( "0.01", "0.1", "1", "10", "100", "1000",
                           expression(10^4), expression(10^5) ) )
    axis( 2, 10^(-2:3), c( "0.01", "0.1", "1", "10", "100","1000" ), las=2 )
    abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
    # Plot the genes, use a different color if they are highly variable
    points( meansGenes, cv2Genes, pch=20, cex=.2,
            col = ifelse( padj < .1, "#C0007090", colGenes ) )
    # Plot/highlight the spike-ins if they are different from the genes
    if (length(meansSp) < length(meansGenes)) {
      points(meansSp, cv2Sp, pch=20, cex=.5, col=colSp)
    }
    # Add the technical noise fit
    xg <- 10^seq( -2, 6, length.out=1000 )
    lines( xg, (a1)/xg + a0, col="#FF000080", lwd=3 )
    # Add a curve showing the expectation for the chosen biological CV^2 thershold
    lines( xg, psia1theta/xg + a0 + minBiolDisp, lty="dashed", col="#C0007090", lwd=3)
  }
  TABLE <- data.frame(Gene = names(meansGenes)[sig], effect.size=res[sig], p.value = p[sig], q.value= padj[sig])
  TABLE <- TABLE[order(-TABLE[,2]),];
  return(TABLE)
}

setwd('../negative')
clusterProfilerPlots(cp.nr.lt.nonglm.all$negative)


closestFactor <- function(x){
  q = ceiling(sqrt(x))
  while(x%%q != 0){
    q = q+1
  }
  return(c(q, x/q))
}

makeUpSet <- function(deg1, deg2, name){
  deg1 <- deg1[deg2$primerid,]
  deg1[,name] <- deg2$coef
  deg1[deg1 > 0] <- 1
  deg1[deg1 < 0] <- 0
  deg1[is.na(deg1)] <- -1
  deg1_temp <- deg1; deg1_temp[deg1_temp < 0] <- 0
  comb <- make_comb_mat(deg1_temp,mode = 'intersect')
  a <- UpSet(comb[,colSums(comb) == 2 & comb[name,] == 1])
  deg1_temp <- 1-deg1; deg1_temp[deg1_temp == 2] <- 0
  comb <- make_comb_mat(deg1_temp,mode = 'intersect')
  b <- UpSet(comb[,colSums(comb) == 2 & comb[name,] == 1])
  return(list(pos = a, neg = b))
}

qvalue_truncp <- function(p, fdr.level = NULL, pfdr = FALSE, lfdr.out = TRUE, pi0 = NULL, ...) {
  # Argument checks
  p_in <- qvals_out <- lfdr_out <- p
  rm_na <- !is.na(p)
  p <- p[rm_na]
  if (min(p) < 0 || max(p) > 1) {
    stop("p-values not in valid range [0, 1].")
  } else if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level > 1)) {
    stop("'fdr.level' must be in (0, 1].")
  }
  p <- p / max(p)
  # Calculate pi0 estimate
  if (is.null(pi0)) {
    pi0s <- pi0est(p, ...)
  } else {
    if (pi0 > 0 && pi0 <= 1)  {
      pi0s = list()
      pi0s$pi0 = pi0
    } else {
      stop("pi0 is not (0,1]")
    }
  }
  
  # Calculate q-value estimates
  m <- length(p)
  i <- m:1L
  o <- order(p, decreasing = TRUE)
  ro <- order(o)
  if (pfdr) {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m / (i * (1 - (1 - p[o]) ^ m))))[ro]
  } else {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m /i ))[ro]
  }
  qvals_out[rm_na] <- qvals
  # Calculate local FDR estimates
  if (lfdr.out) {
    lfdr <- lfdr(p = p, pi0 = pi0s$pi0, ...)
    lfdr_out[rm_na] <- lfdr
  } else {
    lfdr_out <- NULL
  }
  
  # Return results
  if (!is.null(fdr.level)) {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, qvalues = qvals_out,
                   pvalues = p_in, lfdr = lfdr_out, fdr.level = fdr.level,
                   significant = (qvals <= fdr.level),
                   pi0.lambda = pi0s$pi0.lambda, lambda = pi0s$lambda,
                   pi0.smooth = pi0s$pi0.smooth)
  } else {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, qvalues = qvals_out,
                   pvalues = p_in, lfdr = lfdr_out, pi0.lambda = pi0s$pi0.lambda,
                   lambda = pi0s$lambda, pi0.smooth = pi0s$pi0.smooth)
  }
  class(retval) <- "qvalue"
  return(retval)
}
