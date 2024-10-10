if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("PADOG", "GSVA", "AnnotationDbi", "topGO",
                       "pathview", "gage", "globaltest", "limma", "edgeR", "safe",
                       "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db"))


egsea_limma_obj <- function(ct,  meta, control){
  obj <- unique_entrez_mat(ct)
  mat <- obj$ct
  genes <- obj$genes
  meta <- meta[colnames(mat),]
  group <- as.character(meta$condition)
  group <- relevel(as.factor(group), control)
  dge <- DGEList(counts=mat)
  design <- model.matrix(~0+group)
    dge <- calcNormFactors(dge)
    v <- voom(dge, design, plot=TRUE)
    #fit <- lmFit(v, design)
    #fit <- eBayes(fit)
    v$genes <- genes
    v
}

test_limma <- egsea_limma_obj(counts(salmon_res5$dds), colData(salmon_res5$dds), control='control')
gs.annots = buildIdx(entrezIDs=row.names(test_limma$E), species="human", 
                     msigdb.gsets=c("h"), go.part = TRUE)


unique_entrez_mat <- function(df){
  df <- data.frame(df)
  df$genes <- ENSEMBL2ALIAS[row.names(df)]
  bitr_res <- bitr(df$genes, 'ALIAS', 'ENTREZID', org.Hs.eg.db)
  df <- df[df$genes %in% bitr_res$ALIAS,]
  ALIAS2ENTREZ <- tapply(bitr_res$ENTREZID, bitr_res$ALIAS, function(x) x[1])
  df$entrez <- ALIAS2ENTREZ[df$genes]
  temp_df <- df[,!grepl('entrez|genes', colnames(df), perl = T)]
  
  new_df <- do.call(rbind, tapply(temp_df, df$entrez, function(x){
    x_m <- as.matrix(x)
    cv = rowVars(x_m)/rowMeans(x_m)
    x[which(cv == max(cv)),]
  }))
  genes <- data.frame(do.call(rbind, lapply(row.names(new_df), function(x){
    x_m <- as.matrix(temp_df[df$entrez == x,])
    cv = rowVars(x_m)/rowMeans(x_m)
    g = c(x, df[df$entrez == x,'genes'][which(cv == max(cv))])
  })))
  #genes <- do.call(rbind(lapply(row.names(new_df), function(x){c(x, ALIAS2ENTREZ[ALIAS2ENTREZ$ENTREZID == x, 1][1])})))
  list(ct=new_df, genes = genes)
}


msiggs <- buildMSigDBIdx(
  row.names(test_limma$E),
  species = "Homo sapiens",
  geneSets = "h",
  go.part = FALSE,
  min.size = 10
)

gsets <- buildIdx(
  row.names(test_limma$E),
  species = "human", msigdb.gsets = c("h"), kegg.exclude = c("Metabolism"))


contrast.matrix = makeContrasts('grouptreated-groupcontrol',
                                levels=test_limma$design)
gsa = egsea(voom.results = test_limma, contrasts = contrast.matrix, gs.annots = msiggs,
            symbolsMap = test_limma$genes, baseGSEAs =  egsea.base(), report.dir = "./egsea_res",
            sort.by = "p.adj", num.threads = 4, report = T, interactive = T, 
            keep.base=T, keep.set.scores = T)




