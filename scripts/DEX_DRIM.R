read_rnasplice_dex_dtu <- function(dexseq_ds_rds, drim_seq_filt_rds){
  dexseq_ds <- readRDS(dexseq_ds_rds)
  drimseq_ds <- readRDS(drim_seq_filt_rds)
  dex_norm <- cbind(as.data.frame(stringr::str_split_fixed(rownames(counts(dexseq_ds)), ":", 2)), as.data.frame(counts(dexseq_ds, normalized = TRUE))[,1:(nrow(dexseq_ds@colData)/2)])
  colnames(dex_norm) = c("groupID", "featureID", as.character(colData(dexseq_ds)$sample.1)[1:(nrow(dexseq_ds@colData)/2)])
  row.names(dex_norm) = NULL
  obj <- list()
  obj$dexseq <- dexseq_ds
  obj$drimseq <- drimseq_ds
  obj$dex_norm <- dex_norm
  return(obj)
}

get_drim_prop <- function(cds){
  cts = as.matrix(subset(counts(cds), select = -c(gene_id, feature_id)))
  # Summarise count total per gene
  gene.cts = rowsum(cts, counts(cds)$gene_id)
  # Use total count per gene as count per transcript
  total.cts = gene.cts[match(counts(cds)$gene_id, rownames(gene.cts)),]
  # Calculate proportion of transcript in gene
  props = cts/total.cts
  
  colnames(props) = colnames(cts)
  props <- data.frame(props)
  
  props$gene_id <- counts(cds)$gene_id
  props$feature_id <- counts(cds)$feature_id
  rownames(props) = props$feature_id
  props
}

smallProportionSD <- function(d, filter = 0.1) {
  # Generate count table
  props = get_drim_prop(d)
  props = subset(props, select = -c(gene_id, feature_id))
  # Calculate standard deviation
  propSD = sqrt(rowVars(props))
  # Check if standard deviation of per-sample proportions is < 0.1
  propSD < filter
}

smallPropDiff <- function(d, group, filter= 0.2){
  props = get_drim_prop(d)
  props = subset(props, select = -c(gene_id, feature_id))
  lvl = unique(group)
  diff = abs(rowMeans(props[,group == lvl[1]]) - rowMeans(props[,group == lvl[2]]))
  diff < filter
}

drim_stageR <- function(obj, cond_names = c('MF' = 'mouseZygote', 'MU' = 'mouseEgg'), seed = 12345, otherCond = 'mouseZygote', correct_det = F, filter = T){
  no.na <- function(x) ifelse(is.na(x), 1, x) 
  drim <- obj$drimseq
  dex <- obj$dexseq
  samps = data.frame(sample_id = dex@colData$sample.1, group = cond_names[dex@colData$condition])[1:(nrow(dex@colData)/2),]
  design = model.matrix(~group, data = samps)
  set.seed(seed)
  if(correct_det){
    samps = data.frame( sample_id = dex@colData$sample.1, group = cond_names[dex@colData$condition], sensitivity = colSums(dex@assays@data$counts > 0))[1:(nrow(dex@colData)/2),]
    cts = as.matrix(subset(counts(drim), select = -c(gene_id, feature_id)))
    # Summarise count total per gene
    gene.cts = rowsum(cts, counts(drim)$gene_id)[,samps$sample_id]
    samps$sensitivity <- colSums(gene.cts > 0)/mean(colSums(gene.cts > 0))
    #samps$sensitivity <- rnorm(ncol(gene.cts), 1, sd = 0.00001)
    design = model.matrix(~group+sensitivity, data = samps)
  }
  
  print(samps)
  #print(design)
  
  if(filter){
    n = 10#nrow(samps)
    n.small = 5#min(table(samps$group))
    expr_min = 10 #10
    drim <- dmFilter(drim, min_samps_feature_expr = n.small, min_feature_expr = expr_min,
                     min_samps_feature_prop = n.small, min_feature_prop = 0.1, 
                     min_samps_gene_expr = n, min_gene_expr = expr_min)
  }
  print(drim)
  system.time({
    drim <- dmPrecision(drim, design = design, add_uniform=F, prec_subset=0.1)    # Estimate the precision (Higher dispersion is associated with lower precision)
    drim <- dmFit(drim, design = design, add_uniform=F)          # Fit regression coefficients
    drim <- dmTest(drim, coef = paste("group", otherCond, sep = ''))     # Perform null hypothesis testing on the coefficient of interest
  })
  
  res.g = DRIMSeq::results(drim)
  res.g$pvalue <- no.na(res.g$pvalue) 
  #res.t = DRIMSeq::results(drim, level = "feature")
  #res.t$pvalue <- no.na(res.t$pvalue)
  filt = smallProportionSD(drim)
  res.t.filt = DRIMSeq::results(drim, level = "feature")
  res.t.filt$pvalue <- no.na(res.t.filt$pvalue)
  res.t.filt$pvalue[filt] = 1
  res.t.filt$adj_pvalue[filt] = 1
  pScreen = res.g$pvalue
  names(pScreen) = res.g$gene_id
  pConfirmation = matrix(res.t.filt$pvalue, ncol = 1)
  dimnames(pConfirmation) = list(res.t.filt$feature_id, "transcript")
  tx2gene = data.frame(res.t.filt[,c("feature_id", "gene_id")], res.t.filt[,c("feature_id", "gene_id")])
  
  stageRObj = stageRTx(pScreen = pScreen, 
                       pConfirmation = pConfirmation, 
                       pScreenAdjusted = FALSE, 
                       tx2gene = tx2gene[,1:2])
  
  stageRObj = stageWiseAdjustment(stageRObj, method = "dtu", alpha = 0.05)
  drim.padj = getAdjustedPValues(stageRObj, order = FALSE, onlySignificantGenes = T)
  drim.padj = merge(tx2gene, drim.padj, by.x = c("gene_id","feature_id"), by.y = c("geneID","txID"))
  
  res_obj <- list()
  res_obj$drim <- drim
  res_obj$filt <- filt
  res_obj$res.t.filt <- res.t.filt
  res_obj$res.g <- res.g
  res_obj$pScreen <- pScreen
  res_obj$pConfirmation <- pConfirmation
  res_obj$stageRObj <- stageRObj
  res_obj$stageR_res <- drim.padj
  res_obj$samps = samps
  res_obj
}

stageR_dexseqRes <- function(dex){
  #dex = DEXSeqResults(dex_obj, independentFiltering = FALSE)
  qval = perGeneQValue(dex)
  gene_pvals = data.frame(groupID= names(qval), padj=qval)
  #gene_pvals <- read.csv(per_gene, header = T,  sep = '\t')
  t_pvals <- dex$pvalue
  t_pvals[is.na(t_pvals)] <- 1
  #res.t = DRIMSeq::results(drim, level = "feature")
  #res.t$pvalue <- no.na(res.t$pvalue)
  pScreen = gene_pvals$padj
  names(pScreen) = gene_pvals$groupID
  pConfirmation = matrix(t_pvals, ncol = 1)
  dimnames(pConfirmation) = list(dex$featureID, "transcript")
  #View(pConfirmation)
  tx2gene = data.frame(dex[,c("featureID", "groupID" )], dex[,c("featureID", "groupID")])
  #View(tx2gene)
  stageRObj = stageRTx(pScreen = pScreen, 
                       pConfirmation = pConfirmation, 
                       pScreenAdjusted = T, 
                       tx2gene = tx2gene[,1:2])
  
  stageRObj = stageWiseAdjustment(stageRObj, method = "dtu", alpha = 0.05)
  drim.padj = getAdjustedPValues(stageRObj, order = FALSE, onlySignificantGenes = F)
}


refilt_drim_res <- function(drim, group = NULL, prop_diff = 0.2, prop_sd = 0.1){
  res.g = DRIMSeq::results(drim)
  res.g$pvalue <- no.na(res.g$pvalue) 
  #res.t = DRIMSeq::results(drim, level = "feature")
  #res.t$pvalue <- no.na(res.t$pvalue)
  res.t.filt = DRIMSeq::results(drim, level = "feature")
  res.t.filt$pvalue <- no.na(res.t.filt$pvalue)
  filt <- which(rep(F, nrow(res.t.filt)))
  print(length(filt))
  if(prop_diff){
    filt = union(which(smallPropDiff(drim, group = group, filter = prop_diff)), filt)
    print(length(filt))
  }
  if(prop_sd){
    filt = union(which(smallProportionSD(drim, filter = prop_sd)),filt)
    print(length(filt))
  }
  
  res.t.filt$pvalue[filt] = 1
  res.t.filt$adj_pvalue[filt] = 1
  pScreen = res.g$pvalue
  names(pScreen) = res.g$gene_id
  pConfirmation = matrix(res.t.filt$pvalue, ncol = 1)
  dimnames(pConfirmation) = list(res.t.filt$feature_id, "transcript")
  tx2gene = data.frame(res.t.filt[,c("feature_id", "gene_id")], res.t.filt[,c("feature_id", "gene_id")])
  
  stageRObj = stageRTx(pScreen = pScreen, 
                       pConfirmation = pConfirmation, 
                       pScreenAdjusted = FALSE, 
                       tx2gene = tx2gene[,1:2])
  
  stageRObj = stageWiseAdjustment(stageRObj, method = "dtu", alpha = 0.05)
  drim.padj = getAdjustedPValues(stageRObj, order = FALSE, onlySignificantGenes = T)
  drim.padj = merge(tx2gene, drim.padj, by.x = c("gene_id","feature_id"), by.y = c("geneID","txID"))
}


stageR_wrapper <- function(screen, screen_names, confirm, confirm_names, s2c){
  pScreen = screen
  names(pScreen) = screen_names
  pConfirmation = matrix(confirm, ncol = 1)
  dimnames(pConfirmation) = list(confirm_names, "transcript")
  #View(pConfirmation)
  tx2gene = s2c
  #View(tx2gene)
  stageRObj = stageRTx(pScreen = pScreen, 
                       pConfirmation = pConfirmation, 
                       pScreenAdjusted = F, 
                       tx2gene = tx2gene[,1:2])
  
  stageRObj = stageWiseAdjustment(stageRObj, method = "dte", alpha = 0.05)
  drim.padj = getAdjustedPValues(stageRObj, order = FALSE, onlySignificantGenes = F)
}

prop_diff <- function(props, group){
  lvl = unique(group)
  diff = abs(rowMeans(props[,group == lvl[1]], na.rm = T) - rowMeans(props[,group == lvl[2]], na.rm = T))
}


get_dexseq_prop <- function(dex_obj){
  drim.prop = reshape2::melt(dex_obj$dex_norm, id = c("groupID", "featureID"))
  drim.prop = drim.prop[order(drim.prop$groupID, drim.prop$variable, drim.prop$featureID),]

  # Calculate proportions from counts
  system.time({
    drim.prop = drim.prop %>%
     group_by(groupID, variable) %>%
     mutate(total = sum(value)) %>%
     group_by(variable, add=TRUE) %>%
     mutate(prop = value/total)
  })

  # Convert the data.frame to wide-format data with reshape2::dcast
  drim.prop = reshape2::dcast(drim.prop[,c(1,2,3,6)], groupID + featureID ~ variable)

}



plotDEXSeqDTU <- function(expData = NULL, geneID = NULL, samps = NULL, isProportion = FALSE) {
  colnames(expData)[1:2] = c("gid","tid")
  sub = subset(expData, gid == geneID)
  colnames(samps) <- c('sample_id', 'group')
  sub = reshape2::melt(sub, id = c("gid", "tid"))
  sub = merge(samps, sub, by.x = "sample_id", by.y = "variable")
  if(!isProportion) {
    sub$value = log2(sub$value+1)
  }
  
  #clrs = c("dodgerblue3", "maroon2",  "forestgreen", "darkorange1", "blueviolet", "firebrick2",
  # "deepskyblue", "orchid2", "chartreuse3", "gold", "slateblue1", "tomato" , "blue", "magenta", "green3",
  # "yellow", "purple3", "red" ,"darkslategray1", "lightpink1", "lightgreen", "khaki1", "plum3", "salmon")
  clrs = DOT_COLOR
  print(clrs)
  p = ggplot(sub, aes(tid, value, color = group, fill = group)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.8, lwd = 0.5) +
    stat_summary(fun = mean, geom = "point", color = "black", shape = 5, size = 3, position=position_dodge(width = 0.8)) +
    scale_color_manual(values = clrs) + scale_fill_manual(values = clrs) +
    geom_quasirandom(color = "black", size = 1, dodge.width = 0.8) + theme_bw() +
    ggtitle(geneID) + xlab("Transcripts")+theme(axis.text.x = element_text(angle = 15, vjust = 0.5, hjust=0.5, size = 8), axis.title = element_text(size = 15), title = element_text(size = 16))
  
  if(!isProportion) {
    p = p + ylab("log(Expression)")
  } else {
    p = p + ylab("Proportions")
  }
  p
}
