require(GenomicFeatures)
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(PADOG)
  library(GSVA)
  library(org.Hs.eg.db)
  library(ggrepel)
  library(apeglm)
  library(heatmap3)
  library(ggplot2)
  library(ggfortify)
  library(stringr)
  library(RColorBrewer)
  library(MKmisc)
  library(DESeq2)
  library(reticulate)
  library(edgeR)
  library(GSEABase)
  library(limma)
  library(aplot)
  library(reshape2)
  library(data.table)
  library(knitr)
  library(stringr)
  library(rsvd)
  library(pcaMethods)
  library(robust)
  library(MASS)
  library(umap)
  library(tximport)
  require(ReactomePA)
  require(clusterProfiler)
  require(meshes)
  library(msigdbr)
  library(doParallel)
  library(grid)
  library(gridExtra)
  library(pathview)
  library(igraph)
  library(ggraph)
  library(qvalue)
  library(ggrepel)
  library(DRIMSeq)
  library(ggbeeswarm)
  library(stageR)
  library(DEXSeq)
  library(GenomicFeatures)
  library(EGSEA)
})


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


prepareCount <- function(cds, dirc)
{
  cells <- row.names(cds$meta)
  # tpm expression file
  #genes_table <- read.table(paste(FILES,dirc, "genes.csv", sep=""), row.names = 1, header=T, sep = ',')
  
  #if(ribo_filter){
  # genes_table <- genes_table[!grepl('ribosom', genes_table$Gene.description) & !grepl('rRNA', genes_table$Gene.description),]
  #}
  #genes <- row.names(subset(genes_table, Gene.featureType == 'ORF'))
  
  salmon = read_expression(paste(FILES, dirc, 'salmon', sep='/'), mode = 'salmon', tx2gene = paste(FILES, dirc, 'tx2gene.tsv', sep = '/'))
  star = read.csv(paste(FILES, dirc, 'gene_featureCount.txt', sep='/'), header = T, row.names = 1, sep = '\t')
  cds[['salmon']] <- salmon
  tpm <- salmon$gene$abundance[,cells]
  cds[['tpm']] <- list()
  cds[['tpm']][['bio']] <- tpm[which(!grepl("ERCC", row.names(tpm))),]
  cds[['tpm']][['bio']] <- t(t(cds[['tpm']][['bio']])*1e6/colSums(cds[['tpm']][['bio']]))
  cds[['tpm']][['spike']] <- tpm[which(grepl("ERCC", row.names(tpm))),]
  cds[['ct']] <- list(); 
  cds[['ct']][['bio']] <- star[which(!grepl("ERCC", row.names(star))),cells]
  cds[['ct']][['spike']] <- star[which(grepl("ERCC", row.names(star))),cells]
  return(cds)
}

prepareGeneFeatures <- function(cds, gtf_file = NULL, dirc = 'mouse_data',org_type = 'mouse')
{
  if(org_type == 'mouse'){
    org_type = 'mmusculus_gene_ensembl'
    kegg_org = 'mmu'
    
  }else if(org_type == 'rat'){
    org_type = 'rnorvegicus_gene_ensembl'
    kegg_org = 'rno'
  }
  genes <- union(rownames(cds$tpm$bio), rownames(cds$ct$bio))
  gtf = rtracklayer::import(gtf_file)
  gtf = data.frame(gtf[gtf$type == 'gene',])
  gtf <- gtf[,colSums(is.na(gtf)) != nrow(gtf)]
  gtf <- data.frame(row.names = gtf$gene_id, gtf)
  gtf$gene_short_name <- gtf$gene
  gtf <- gtf[genes,]
  cds$features <- gtf
  cds$features_tpm <- cds$features[row.names(cds$tpm$bio),]
  cds$features_ct <- cds$features[row.names(cds$ct$bio),]
  


  return(cds)
}


get_organism_items <- function(organisms){
  if(organisms == 'mouse'){
    require(org.Mm.eg.db)
    orgdb = org.Mm.eg.db
    orgabv = 'mmu'
    orgname = "Mus musculus"
  }else if(organisms == 'rat'){
    require(org.Rn.eg.db)
    orgdb = org.Rn.eg.db
    orgabv = 'rno'
    orgname = "Rattus norvegicus"
  }
  else if(organisms == 'celegans'){
    require(org.Ce.eg.db)
    orgdb = org.Ce.eg.db
    orgabv = 'cel'
    orgname = 'Caenorhabditis elegans'
  }
  else if(organisms == 'human'){
    require(org.Hs.eg.db)
    orgdb = org.Hs.eg.db
    orgabv = 'hsa'
    orgname = 'Homo sapiens'
  }
  return(list(orgdb = orgdb, orgkegg = orgabv, orgname = orgname))
}


enrich_CP <- function(ora_genes, organisms, n_type = 'ENSEMBL',universe = NULL, classic = T, GO_BP_only = F, enrich_all = T, Msig = NULL, alpha = 0.5, combine = T, simple_combine = T, full_combine = T){
  require(ReactomePA)
  require(clusterProfiler)
  require(org.Sc.sgd.db)
  require(meshes)
  require(dplyr)
  items = get_organism_items(organisms = organisms)
  orgabv = items$orgkegg
  orgname = items$orgname
  orgdb = items$orgdb
  
  if(n_type == 'ENSEMBL'){
    ora_genes <- ENSEMBL2ALIAS[ora_genes]
    universe <- ENSEMBL2ALIAS[universe]
    n_type = 'ALIAS'
  }
  
  oraL <- tryCatch({
    egdf <- bitr(ora_genes, n_type, 'ENTREZID', orgdb) %>% distinct(eval(as.name('ALIAS')), .keep_all = T) %>% data.frame(row.names = 1)
    egdf$ENTREZID
  },error=function(cond){return(NULL)})
  universe <-  tryCatch({
    egdf <- bitr(universe, n_type, 'ENTREZID', orgdb) %>% distinct(eval(as.name('ALIAS')), .keep_all = T) %>% data.frame(row.names = 1)
    egdf$ENTREZID
  },error=function(cond){return(NULL)})
  
  if(is.null(oraL) || length(oraL) == 0)
  {return(NULL)}
  
  if(organisms == 'celegans'){
    oraL_kegg <- paste('CELE_', WORM.GENES[ora_genes,]$Sequence.Name, sep ='')
  }else{
    oraL_kegg <- oraL
  }
  #return(list(g = oraL, u=universe))
  
  GSE_results <- list()
  # GO Enrichment
  if(GO_BP_only | classic){
    GSE_results[['GO_BP_ora']] <- tryCatch({setReadable(enrichGO(gene = oraL,
                                                                 universe      = universe,
                                                                 OrgDb         = orgdb,
                                                                 ont           = "BP",
                                                                 pAdjustMethod = "BH",
                                                                 pvalueCutoff  = alpha, maxGSSize = 500,minGSSize = 10, 
                                                                 qvalueCutoff = alpha), OrgDb = orgdb)},error=function(cond){return(NULL)})
  }
  if(!GO_BP_only & classic){
    GSE_results[["WKP_ora"]] <- tryCatch({setReadable(enrichWP(oraL, organism = orgname, maxGSSize = 500, minGSSize = 10, universe = universe, pvalueCutoff  = alpha, qvalueCutoff = alpha ), OrgDb = orgdb)},error=function(cond){return(NULL)})
    
    GSE_results[['GO_CC_ora']] <- tryCatch({setReadable(enrichGO(gene = oraL,
                                                                 universe      = universe,
                                                                 OrgDb         = orgdb,
                                                                 ont           = "CC",
                                                                 pAdjustMethod = "BH",
                                                                 pvalueCutoff  = alpha,
                                                                 qvalueCutoff = alpha,maxGSSize = 500,minGSSize = 10, 
    ), OrgDb = orgdb)},error=function(cond){return(NULL)})
    
    GSE_results[['GO_MF_ora']] <- tryCatch({setReadable(enrichGO(gene = oraL,
                                                                 universe      = universe,
                                                                 OrgDb         = orgdb,
                                                                 ont           = "MF",
                                                                 pAdjustMethod = "BH",
                                                                 pvalueCutoff  = alpha, maxGSSize = 500,minGSSize = 10, 
                                                                 qvalueCutoff = alpha), OrgDb = orgdb)},error=function(cond){return(NULL)})
    
    
    #KEGG Enrichment
    
    GSE_results[['KEGG_ora']] <- tryCatch({setReadable(enrichKEGG(gene = oraL, universe = universe,
                                                                  organism     = orgabv, maxGSSize = 500,minGSSize = 10, 
                                                                  pvalueCutoff  = alpha,qvalueCutoff = alpha), OrgDb = orgdb, keyType = 'ENTREZID')},error=function(cond){return(NULL)})
    
    
    GSE_results[['MKEGG_ora']] <- tryCatch({setReadable(enrichMKEGG(gene = oraL, universe = universe, maxGSSize = 500,minGSSize = 10, 
                                                                    organism = orgabv, pvalueCutoff  = alpha,qvalueCutoff = alpha), OrgDb = orgdb, keyType = 'ENTREZID')},error=function(cond){return(NULL)})
    
    # Reactome Enrichment
    GSE_results[['REACT_ora']] <- tryCatch({setReadable(enrichPathway(gene=oraL,organism = organisms, maxGSSize = 500,minGSSize = 10,  universe = universe, pvalueCutoff  = alpha,qvalueCutoff = alpha), OrgDb = orgdb, keyType = 'ENTREZID')},error=function(cond){return(NULL)})
  }
  
  Msig_res <- NULL
  if(!is.null(Msig)){
    Msig_res <- lapply(Msig, function(x){
      cat = strsplit(x, '-')[[1]][1]
      sub_cat = ''
      if(length(strsplit(x, '-')[[1]]) > 1){
        sub_cat = strsplit(x, '-')[[1]][2]
      }
      df = get_msig(organisms, cat = cat, sub_cat = sub_cat, gmt_dir = './dataset/Msigdb/') 
      tryCatch({setReadable(enricher(oraL, TERM2GENE=df, maxGSSize = 500,minGSSize = 10,  universe = universe, pvalueCutoff = alpha, qvalueCutoff = alpha), OrgDb = orgdb, keyType = 'ENTREZID')},error=function(cond){return(NULL)})
    })
  }
  
  if(!is.null(Msig_res)){
    names(Msig_res) <- Msig
    for(n in names(Msig_res)){
      GSE_results[[n]] <- Msig_res[[n]]
    }
  }
  if(combine){
    if(simple_combine == T & classic){
      combined <- GSE_results[['GO_BP_ora']]
      res_ <- c(c("WKP_ora", 'GO_BP_ora','KEGG_ora','MKEGG_ora','REACT_ora'), Msig)
      combined@result <- do.call(rbind, lapply(res_, FUN = function(x){
        r=GSE_results[[x]]@result
        if(ncol(r) == 11){
          r <- r[,3:11]
        }
        r
      }))
      row.names(combined@result) <- combined@result$ID
      combined@geneSets <- do.call(c, lapply(res_, FUN = function(x){GSE_results[[x]]@geneSets}))
      names(combined@geneSets) <- do.call(c, lapply(res_, FUN = function(x){names(GSE_results[[x]]@geneSets)}))
      combined@result$old_qvalue <- combined@result$qvalue
      combined@result$qvalue <- qvalue(combined@result$pvalue)$qvalue
      GSE_results[['combined']] <- combined
    }
    if(full_combine == T){
      all_sets <- NULL
      all_sets_n <- NULL
      if(classic){
        wiki <- data.frame(clusterProfiler:::get_wp_data(orgname))
        wiki_t2g <- wiki[, c('wpid', 'gene')]
        colnames(wiki_t2g) <- c('term', 'gene')
        wiki_t2n <- unique(wiki[, c('wpid', 'name')])
        colnames(wiki_t2n) <- c('term', 'name')
        row.names(wiki_t2n) <- wiki_t2n$wpid
        react <- as.list(ReactomePA:::get_Reactome_DATA(organisms))
        go_bp <- as.list(clusterProfiler:::get_GO_data(orgdb, 'BP', "ENTREZID"))
        kegg <- as.list(clusterProfiler:::prepare_KEGG(orgabv, "KEGG", "ncbi-geneid"))
        go_kegg_react_list <- c(react$PATHID2EXTID, go_bp$PATHID2EXTID, kegg$PATHID2EXTID)
        go_kegg_react_names <- c(react$PATHID2NAME, go_bp$PATHID2NAME, kegg$PATHID2NAME)
        go_kegg_react_P2G <- data.frame(do.call(rbind, lapply(names(go_kegg_react_list), FUN = function(x){
          cbind(rep(x, length(go_kegg_react_list[[x]])), go_kegg_react_list[[x]])
        })))
        colnames(go_kegg_react_P2G) <- c('term', 'gene')
        
        go_kegg_react_P2N <- data.frame(term = names(go_kegg_react_names), name = go_kegg_react_names)
        
        all_sets <- rbind(go_kegg_react_P2G, wiki_t2g)
        all_sets_n <- rbind(go_kegg_react_P2N, wiki_t2n)
      }
      
      Msig_df <- NULL
      if(!is.null(Msig)){
        Msig_df <- do.call(rbind, lapply(Msig, function(x){
          cat = strsplit(x, '-')[[1]][1]
          sub_cat = ''
          if(length(strsplit(x, '-')[[1]]) > 1){
            sub_cat = strsplit(x, '-')[[1]][2]
          }
          df = get_msig(organisms, cat = cat, sub_cat = sub_cat, gmt_dir = './dataset/Msigdb/') 
          colnames(df) <- c('term', 'gene')
          df
        }))
        Msig_df <- data.frame(Msig_df)
        View(Msig_df)
        colnames(Msig_df) <- c('term', 'gene')
        all_sets <- rbind(all_sets, Msig_df)
        all_sets_n <- rbind(all_sets_n, data.frame(term=unique(Msig_df[,1]), name=unique(Msig_df[,1])))
      }
      if(!is.null(all_sets)){
        GSE_results[['combined_full']] <- tryCatch({setReadable(enricher(oraL, TERM2GENE = all_sets,
                                                                         TERM2NAME = all_sets_n, 
                                                                         maxGSSize = 500,minGSSize = 10,  
                                                                         universe = universe, pvalueCutoff = alpha, 
                                                                         qvalueCutoff = alpha), OrgDb = orgdb, keyType = 'ENTREZID')},error=function(cond){return(NULL)})}
    }
  }
  
  return(GSE_results)
}


gse_CP <- function( organisms, logFC=NULL, n_type = 'ENSEMBL', classic = T, simplify_go = T, combine = F, simple_combine = F, full_combine = F, Msig = NULL, alpha = 1, disease= F){
  require(ReactomePA)
  require(clusterProfiler)
  require(org.Sc.sgd.db)
  require(meshes)
  require(dplyr)
  if(simple_combine){
    alpha = 1
  }
  if(!classic){
    simple_combine = F
  }
  
  items = get_organism_items(organisms = organisms)
  orgabv = items$orgkegg
  orgname = items$orgname
  orgdb = items$orgdb
  GSE_results <- list()
  gse_list0 <- NULL
  gse_list <- NULL
  if(!is.null(logFC)){ 
    gse_list <- logFC
    if(n_type == 'ENSEMBL'){
      gse_list <- tapply(gse_list, ENSEMBL2ALIAS[names(gse_list)], FUN = function(x){
        x[which(abs(x) == max(abs(x)))]
      })
      n_type ='ALIAS'
      gnames <- names(gse_list)
      gse_list <- as.vector(gse_list)
      names(gse_list) <- gnames
    }
    gse_list <- sort(gse_list, T)
    gse_list0 <- sort(gse_list0, T)
    egdf <- bitr(names(gse_list), n_type, 'ENTREZID', orgdb) %>% distinct(eval(as.name('ALIAS')), .keep_all = T) %>% data.frame(row.names = 1)
    gse_list <- gse_list[intersect(names(gse_list), row.names(egdf))]
    names(gse_list) <- egdf[names(gse_list),]$ENTREZID
    #return(gse_list)
    if(classic){
      GSE_results[["WKP_gse"]] <- setReadable(gseWP(gse_list, eps = 0, organism = orgname, minGSSize = 10, nPermSimple = 100000, maxGSSize = 500, pvalueCutoff=alpha), OrgDb = orgdb, keyType = 'ENTREZID')
      #GSE_results[["BIO_gse"]] <- setReadable0(GSEA(gse_list, TERM2GENE = BIOCYC[,c(1,2)], TERM2NAME = BIOCYC[,c(1,3)], nPerm = 1000, minGSSize = 5, maxGSSize = 500), gene2symbol = entrez2symbol, keyType = 'ENTREZ')
      GSE_results[['GO_BP_gse']] <- setReadable(gseGO(geneList = gse_list,
                                                      OrgDb        = orgdb, 
                                                      keyType = "ENTREZID", nPermSimple = 10000,
                                                      ont          = "BP",eps = 0,
                                                      minGSSize    = 10,
                                                      maxGSSize    = 500,
                                                      pvalueCutoff = alpha,
                                                      verbose      = FALSE, by = 'fgsea'), OrgDb = orgdb, keyType = 'ENTREZID')
      
      GSE_results[['GO_CC_gse']] <- setReadable(gseGO(geneList = gse_list,
                                                      OrgDb        = orgdb,
                                                      keyType = "ENTREZID",
                                                      ont          = "CC",nPermSimple = 10000,
                                                      minGSSize    = 10,eps = 0,
                                                      maxGSSize    = 500,
                                                      pvalueCutoff = alpha,
                                                      verbose      = FALSE), OrgDb = orgdb, keyType = 'ENTREZID')
      
      GSE_results[['GO_MF_gse']] <- setReadable(gseGO(geneList     = gse_list,
                                                      OrgDb        = orgdb,keyType = "ENTREZID",
                                                      ont          = "MF",
                                                      minGSSize    = 10,eps = 0,
                                                      maxGSSize    = 500,nPermSimple = 10000,
                                                      pvalueCutoff = alpha,
                                                      verbose      = FALSE), OrgDb = orgdb, keyType = 'ENTREZID')
      GSE_results[['KEGG_gse']]<- setReadable(gseKEGG(geneList = gse_list, 
                                                      organism     = orgabv,
                                                      minGSSize = 10, eps = 0,
                                                      maxGSSize = 500,nPermSimple = 10000,
                                                      pvalueCutoff = alpha,
                                                      verbose      = FALSE, keyType = 'ncbi-geneid'),  OrgDb = orgdb, keyType = 'ENTREZID')
      GSE_results[['MKEGG_gse']] <- setReadable(gseMKEGG(gene = gse_list,  minGSSize = 10, eps = 0, nPermSimple = 10000, maxGSSize = 500, organism = orgabv, pvalueCutoff=alpha),  OrgDb = orgdb, keyType = 'ENTREZID')
      GSE_results[['REACT_gse']] <- setReadable(gsePathway(geneList =  gse_list, organism = organisms, nPermSimple = 10000, minGSSize = 10, maxGSSize = 500, 
                                                           pvalueCutoff=alpha,eps = 0,
                                                           pAdjustMethod="BH", verbose=FALSE), OrgDb = orgdb, keyType = 'ENTREZID')
    }
    Msig_res <- NULL
    if(!is.null(Msig)){
      Msig_res <- lapply(Msig, function(x){
        cat = strsplit(x, '-')[[1]][1]
        sub_cat = ''
        if(length(strsplit(x, '-')[[1]]) > 1){
          sub_cat = strsplit(x, '-')[[1]][2]
        }
        df = get_msig(organisms, cat = cat, sub_cat = sub_cat, gmt_dir = './dataset/Msigdb/') 
        tryCatch({setReadable(GSEA(
          gse_list,
          exponent = 1,
          minGSSize = 10,
          maxGSSize = 500,
          eps = 0,nPermSimple = 10000,
          pvalueCutoff = alpha,
          pAdjustMethod = "BH",
          TERM2GENE = df,
          verbose = TRUE,
          seed = FALSE,
          by = "fgsea"), OrgDb = orgdb, keyType = 'ENTREZID')},error=function(cond){return(NULL)})
      })
    }
    
    if(!is.null(Msig_res)){
      names(Msig_res) <- Msig
      for(n in names(Msig_res)){
        GSE_results[[n]] <- Msig_res[[n]]
      }
    }
    #for(r in names(GSE_results)){
    #GSE_results[[r]]@result <- subset(GSE_results[[r]]@result, qvalue < 0.05)
    #}
    if(combine){
      if(simple_combine == T){
        combined <- GSE_results[['GO_BP_gse']]
        res_ <- c(c("WKP_gse", 'GO_BP_gse','KEGG_gse','REACT_gse'), Msig)
        combined@result <- do.call(rbind, lapply(res_, FUN = function(x){GSE_results[[x]]@result}))
        row.names(combined@result) <- combined@result$ID
        combined@geneSets <- do.call(c, lapply(res_, FUN = function(x){GSE_results[[x]]@geneSets}))
        names(combined@geneSets) <- do.call(c, lapply(res_, FUN = function(x){names(GSE_results[[x]]@geneSets)}))
        combined@result$old_qvalue <- combined@result$qvalue
        combined@result$qvalue <- qvalue(combined@result$pvalue)$qvalue
        GSE_results[['combined']] <- combined
        GSE_results[['combined_up']] <- combined
        GSE_results[['combined_up']]@result <- subset(GSE_results[['combined_up']]@result, NES > 0)
        GSE_results[['combined_down']] <- combined
        GSE_results[['combined_down']]@result <- subset(GSE_results[['combined_down']]@result, NES < 0)
      }
      if(full_combine == T){
        all_sets <- NULL
        all_sets_n <- NULL
        if(classic){
          wiki <- data.frame(clusterProfiler:::get_wp_data(orgname))
          wiki_t2g <- wiki[, c('wpid', 'gene')]
          colnames(wiki_t2g) <- c('term', 'gene')
          wiki_t2n <- unique(wiki[, c('wpid', 'name')])
          colnames(wiki_t2n) <- c('term', 'name')
          row.names(wiki_t2n) <- wiki_t2n$wpid
          react <- as.list(ReactomePA:::get_Reactome_DATA(organisms))
          go_bp <- as.list(clusterProfiler:::get_GO_data(orgdb, 'BP', "ENTREZID"))
          kegg <- as.list(clusterProfiler:::prepare_KEGG(orgabv, "KEGG", "ncbi-geneid"))
          go_kegg_react_list <- c(react$PATHID2EXTID, go_bp$PATHID2EXTID, kegg$PATHID2EXTID)
          go_kegg_react_names <- c(react$PATHID2NAME, go_bp$PATHID2NAME, kegg$PATHID2NAME)
          go_kegg_react_P2G <- data.frame(do.call(rbind, lapply(names(go_kegg_react_list), FUN = function(x){
            cbind(rep(x, length(go_kegg_react_list[[x]])), go_kegg_react_list[[x]])
          })))
          colnames(go_kegg_react_P2G) <- c('term', 'gene')
          
          go_kegg_react_P2N <- data.frame(term = names(go_kegg_react_names), name = go_kegg_react_names)
          
          all_sets <- rbind(go_kegg_react_P2G, wiki_t2g)
          all_sets_n <- rbind(go_kegg_react_P2N, wiki_t2n)
        }
        Msig_df <- NULL
        if(!is.null(Msig)){
          Msig_df <- do.call(rbind, lapply(Msig, function(x){
            cat = strsplit(x, '-')[[1]][1]
            sub_cat = ''
            if(length(strsplit(x, '-')[[1]]) > 1){
              sub_cat = strsplit(x, '-')[[1]][2]
            }
            df = get_msig(organisms, cat = cat, sub_cat = sub_cat, gmt_dir = './dataset/Msigdb/') 
            colnames(df) <- c('term', 'gene')
            df
          }))
          Msig_df <- data.frame(Msig_df)
          colnames(Msig_df) <- c('term', 'gene')
          all_sets <- rbind(all_sets, Msig_df)
          all_sets_n <- rbind(all_sets_n, data.frame(term=unique(Msig_df[,1]), name=unique(Msig_df[,1])))
        }
        GSE_results[['combined_full']] <- setReadable(GSEA(
          gse_list,
          exponent = 1,
          minGSSize = 10,
          maxGSSize = 500,
          eps = 0,nPermSimple = 10000,
          pvalueCutoff = alpha,
          pAdjustMethod = "BH",
          gson = NULL,
          TERM2GENE = all_sets,
          TERM2NAME = all_sets_n,
          verbose = TRUE,
          seed = FALSE,
          by = "fgsea",
        ), OrgDb = orgdb, keyType = 'ENTREZID')
      }
    }
    
  }
  
  
  
  GSE_results[['gfc0']] <- gse_list0
  GSE_results[['gfc']] <- gse_list
  
  return(GSE_results)
}

addSmallLegend <- function(myPlot, pointSize = 1.5, textSize = 7, spaceLegend = 0.3) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

custom_cnet_plot <- function(cp_res, top_n_cat = 10, category = NULL, gene_color = NULL, gene_color2 = NULL, layout = 'fr', color_cat_pval = F){
  cp_res1 <- cp_res
  cp_res1@result <- cp_res1@result[category,]
  cp_res1@pvalueCutoff <- 1
  cp_res1@qvalueCutoff <- 1
  pal <- c(rev(c(colorRampPalette(brewer.pal(9, 'Blues'))(32)[1:16],colorRampPalette(brewer.pal(9, 'Blues'))(128)[64:128])), 
           c(colorRampPalette(brewer.pal(9, 'Reds'))(32)[1:16],colorRampPalette(brewer.pal(9, 'Reds'))(128)[64:128]))
  plot1 <- cnetplot( cp_res1, showCategory = top_n_cat, foldChange = gene_color, layout = layout, node_label = 'none')+
    scale_colour_gradientn(name = expression('Log'['2']*'FC'), limits= c(-4, 4), 
                           colours = pal, 
                           na.value = 'black')
  cat_df <- plot1$data[c(1:max(top_n_cat, length(category))),]
  cat_df$qval <- cp_res@result[match(cat_df$name, cp_res@result$Description),]$qvalue
  if(!color_cat_pval){
    plot1 <- plot1+ggnewscale::new_scale_color() +
      ggraph::geom_node_point(aes_(size=~size), data = cat_df, color = 'black') + 
      scale_size(limits = c(1,150), breaks = c(10,20,40,80), range = c(1,5))
  }
  plot1$data$name <- stringr::str_wrap(plot1$data$name, 25)
  plot1 <- plot1+ggraph::geom_node_text(aes_(label=~name), data = plot1$data[c(1:max(top_n_cat, length(category))),], size = 3, bg.color = "white", repel=TRUE)
  
  
  if(!is.null(gene_color2)){
    plot2 <- plot1
    plot2$data$color[-c(1:max(top_n_cat, length(category)))] <- gene_color2[plot1$data$name[-c(1:max(top_n_cat, length(category)))]]
    
  }
  if(color_cat_pval){
    plot1 <- plot1+ggnewscale::new_scale_color() +
      ggraph::geom_node_point(aes_(color=~qval, size=~size), data = cat_df) + 
      scale_size(limits = c(1,150), breaks = c(10,20,40,80), range = c(1,5)) +
      scale_colour_gradientn(name = "qvalue", colours = colorRampPalette(rev(brewer.pal(9,  'Purples')))(255)[0:200], limits= c(5e-4, 0.5), breaks = c( 5e-4, 5e-3, 5e-2, 5e-1), trans = 'log10', oob = scales::oob_squish)
    #theme(plot.margin=unit(c(0,0,0,0),"mm"), aspect.ratio = 1)
    if(!is.null(gene_color2)){
      plot2 <- plot2+ggnewscale::new_scale_color() +
        ggraph::geom_node_point(aes_(color=~qval, size=~size), data = cat_df) + 
        scale_size(limits = c(1,150), breaks = c(10,20,40,80), range = c(1,5)) +
        scale_colour_gradientn(name = "qvalue", colours = colorRampPalette(rev(brewer.pal(9,  'Purples')))(255)[0:200], limits= c(5e-4, 0.5), breaks = c(5e-4, 5e-3, 5e-2, 5e-1), trans = 'log10', oob = scales::oob_squish)
      
    }
  }
  return(list(plot1 = plot1, plot2 = plot2))
}

cp_tree_ridge_plot <- function(res, n_cat = 50, nclust = 5, alpha = 0.1, by_keys = NULL){
  require(ggtree)
  require(enrichplot)
  res@result <- subset(res@result, qvalue < alpha)
  if(sum(grepl('GO', res@result$ID))){
    res@ontology <- 'BP'
    #print(res@setType)
    res <- clusterProfiler::simplify(res)
  }
  
  gcolors <- paletteer_d("ggsci::springfield_simpsons")[1:nclust]
  if(!is.null(by_keys)){
    sets <- find_key_gs(res@result, keys = c('(cancer|carcinoma)', 
                                             '(cycle|mitoti|mitosi)',
                                             'HALLMARK', 'repair', 'ribos', 'death', 'RNA',
                                             'histon', 'chromat', 'methy', 'DNA',
                                             'colorec'), key_length = c(8, 7,5,2,2,2, 2, 2, 2, 2, 2, 2))
    res@result <- res@result[sets,]
  }
  
  res@result$Description[nchar(res@result$Description) > 60] <- res@result$ID[nchar(res@result$Description) > 60]
  
  
  res@result <- res@result[tapply(res@result$ID, res@result$Description, function(x) x[1]),]
  res <- enrichplot::pairwise_termsim(res)
  #res@result <- res@result[tapply(res@result$ID, res@result$Description, function(x) x[order(res@result[x, 'qvalue'])][1]),]
  View(res@result)
  
  res_tree <- addSmallLegend(enrichplot::treeplot(res, showCategory = min(n_cat, nrow(res@result)), 
                                                  geneClusterPanel = 'pie', 
                                                  cluster.params = list(method = "ward.D2", n = nclust, color = gcolors, label_words_n = 4,label_format = 25), 
                                                  color = NULL, offset_tiplab = 0.8, fontsize = 3)+
                               geom_tiplab(offset = 0.8, hjust = 0, 
                                           show.legend = FALSE, 
                                           align=TRUE, size = 2.3)+xlim(c(0,20))+scale_size(name = "number of genes",
                                                                                            range = c(1, 4)))
  tree <- res_tree
  
  res_tree$layers[c(7,8)] <- NULL
  res_tree$layers[c(3,4)] <- NULL

  res_ridge <- addSmallLegend(ridgeplot(res, showCategory = min(n_cat, nrow(res@result)))+scale_fill_viridis_c()+
                                theme(axis.title.y=element_blank(),axis.text.y=element_blank())+xlim(c(-4, 4))+geom_vline(xintercept=0, linetype="dashed", color = "red")+xlab('Log2FoldChange'), pointSize = 3, textSize = 10, spaceLegend = 0.9)
  
  res_ridge$data$label <- res_ridge$data$category
  #return(list(ridge = res_ridge, tree=res_tree))
  ridge_tree <- res_ridge  %>% insert_left(res_tree, width = 3)
  return(list(tree = res_tree, ridge = test_ridge, both = ridge_tree))
  #res_ridge  %>% insert_left(res_tree, width = 2)
  
  
  
}

tf_cnet_plot <- function(res, DE){
  new_res <- res
  sig_gs <- subset(res@result, qvalue < 0.05)
  gs <- sig_gs$ID
  gs_genes <- res@geneSets[gs]
  sig_gs$core_enrichment <- do.call(rbind, lapply(gs_genes, function(x){
    sub_DE <- DE[match(x, DE$entrez),]
    sig_g <- paste(subset(sub_DE, padj < 0.05)$gene_short_name, collapse = '/')
  }))
  new_res@result <- sig_gs
  new_res
}




volcanoplot <- function(res, pval_col = 'fdr', fc_col = 'Log2FC', lfcthresh = F, top_genes = NULL, alpha = 0.05, fc = log2(2), main = NULL){
  res <- subset(res, !is.na(res[[fc_col]]) & !is.na(res[[pval_col]]))
  res$Log2FC <- res[[fc_col]]
  res$padjust <- res[[pval_col]]
  if(is.null(top_genes)){
    top_up_genes <- row.names(subset(res, padjust < alpha)[order(subset(res, padjust < 0.05)$Log2FC, decreasing = T),])[1:20]
    top_down_genes <- row.names(subset(res, padjust < alpha)[order(subset(res, padjust < 0.05)$Log2FC),])[1:20]
    top_genes <- c(top_up_genes, top_down_genes)
  }
  res$deg <- "NoDE"
  if(lfcthresh){
    fc <- 0
  }
  res[res$padjust < alpha & res$Log2FC > fc,]$deg <- "up-regulated"
  res[res$padjust < alpha & res$Log2FC < -fc,]$deg <- "down-regulated"
  res[!row.names(res) %in% top_genes,]$gene_short_name <- NA
  res[res$deg == 'NoDE',]$gene_short_name <- NA
  #res[grepl('LOC',row.names(res)) | grepl('Rik', row.names(res)) | grepl('Gm', row.names(res)),]$gene_short_name <- NA
  p <- ggplot(data=res, aes(x=Log2FC, y=-log10(padjust), color = deg, label=gene_short_name)) + geom_point(size = 0.2)+geom_label_repel(size=3, max.overlaps = 300) 
  theme <- theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 16), 
    axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 16), 
    legend.position = "none",
  )
  p <- p+theme+xlim(min(res$Log2FC),max(res$Log2FC))+ ylim(0, 13)+ggtitle(main)+ylab(expression('-Log'['10']*'Pvalue'))+xlab(expression('-Log'['2']*'FC'))+scale_color_manual(values = c( "#6f95e6", "grey","#a4302a"))
  return(p)
}


geneHeatmap <- function(cds, exprs_mat=NULL, tpm = F, scale = F,
                        genes = NULL, clusterGenes = T,
                        dist_row = 'pearson', dist_col = 'euclidean',
                        main = "Gene Expression Heatmap",
                        fontsize = 5, legend = T, row_color = NULL,
                        annotations, clust_n = 1, mmethod = 'ward.D2') {
  
  if(is.null(exprs_mat)){
    exprs_mat <- as.matrix(log2(t(t(cds)/estimateSizeFactorsForMatrix(as.matrix(round(cds))))+1))
  }
  if(scale){
    exprs_mat <- t(scale(t(exprs_mat)))
  }else{
    exprs_mat <- log2(exprs_mat+1)
  }
  if(tpm){
    cluster_mat <- log(cds+1)
  }else{
    cluster_mat <- DESeq2::varianceStabilizingTransformation(as.matrix(round(cds)))
  }
  if(!is.null(genes)){
    cluster_mat = cluster_mat[genes,]
    exprs_mat <- exprs_mat[genes,]
  }
  
  
  DISTMETHOD_row = dist_row
  DISTMETHOD_col = dist_col
  dist_func_col = function(df){
    #df <- 2**df-1
    #print(dim(df))
    #mat <- DESeq2::varianceStabilizingTransformation(as.matrix(round(t(df))))
    if (DISTMETHOD_col %in% c("pearson", "spearman", "kandall", "cosine", "mcd", "ogk"))
    {
      
      #mat <- df
      require(MKmisc)
      
      cdist <- corDist(t(cluster_mat)[rownames(df),], method = DISTMETHOD_col)
      
      #print("calculating correlation distance")
    }
    else
    {
      cdist <- dist(t(cluster_mat)[rownames(df),], method = DISTMETHOD_col)
      #print("calculating non-correlation distance")
    }
    print(DISTMETHOD_col)
    print(dim(as.matrix(cdist)))
    cdist
  }
  dist_func_row = function(df){
    #df <- 2**df-1
    #mat <- DESeq2::varianceStabilizingTransformation(as.matrix(round(df)))
    if (DISTMETHOD_row %in% c("pearson", "spearman", "kandall", "cosine", "mcd", "ogk"))
    {
      #mat <- DESeq2::varianceStabilizingTransformation(as.matrix(df))
      #mat <- df
      require(MKmisc)
      cdist <- corDist(cluster_mat, method = DISTMETHOD_row)
      
      #print("calculating correlation distance")
    }
    else
    {
      cdist <- dist(cluster_mat, method = DISTMETHOD_row)
      #print("calculating non-correlation distance")
    }
    print(DISTMETHOD_row)
    print(dim(as.matrix(cdist)))
    cdist
  }
  #View(as.matrix(cdist))
  #print(color)
  anno_color <-
    list(
      condition= c('lightblue', 'maroon'),
      logfc = c("negative" = 'red', 'positive' = 'green')
      
    )
  if (nrow(cds) >= 0) {
    show_rown <- F
  } else{
    show_rown <- T
  }
  
  annotations = annotations[,'condition', drop = F]
  if(!is.null(row_color)){
    #logfc <- c("FALSE" = 'negative', 'TRUE' = 'positive')[as.character(row_color > 0)]
    #names(logfc) <- names(row_color)
    logfc = data.frame(row.names = names(row_color), pval = row_color)
  }else{
    logfc = NULL
  }
  if(clust_n == 'auto'){
    clust_n = NbClust(data = cds, diss = cdist, distance = NULL, method = mmethod, index = 'all', max.nc = 50)$Best.partition
  }
  breaksList = seq(0, 13, by = 1)
  color_pal <- colorRampPalette(c('black','red','gold'))(12)
  color_pal <- colorRampPalette(c('blue','white','red'))(12)
  heat_obj <- Heatmap(matrix = exprs_mat, col = color_pal, name = 'Gene Heatmap',
                      na_col = "grey",
                      color_space = "LAB",
                      rect_gp = gpar(col = NA),
                      border = T,
                      border_gp = gpar(col = "black"),
                      cell_fun = NULL,
                      layer_fun = NULL,
                      row_title = character(0),
                      row_title_side = c("left"),
                      row_title_gp = gpar(fontsize = 8),
                      row_title_rot = 0,
                      column_title = character(0),
                      column_title_side = c("top"),
                      column_title_gp = gpar(fontsize = 13.2),
                      column_title_rot = 0,
                      
                      cluster_rows = TRUE,
                      cluster_row_slices = F,
                      clustering_distance_rows = dist_func_row,
                      clustering_method_rows = "ward.D2",
                      row_dend_side = c("left"),
                      row_dend_width = unit(50, "mm"),
                      show_row_dend = F,
                      row_dend_gp = gpar(),
                      
                      cluster_columns = F,
                      cluster_column_slices =  F,
                      clustering_distance_columns = dist_func_col,
                      clustering_method_columns = "ward.D2",
                      column_dend_side = c("top"),
                      column_dend_height = unit(10, "mm"),
                      show_column_dend = F,
                      column_dend_gp = gpar(),
                      
                      row_order = NULL,
                      column_order = NULL,
                      
                      row_names_side = c( "left"),
                      show_row_names = F,
                      row_names_max_width = unit(6, "cm"),
                      row_names_gp = gpar(fontsize = 8),
                      row_names_rot = 0,
                      row_names_centered = FALSE,
                      column_labels = colnames(cds),
                      column_names_side = c("bottom"),
                      show_column_names = F,
                      column_names_max_height = unit(8, "cm"),
                      column_names_gp = gpar(fontsize = 7),
                      column_names_rot = 270,
                      column_names_centered = FALSE,
                      
                      
                      bottom_annotation = HeatmapAnnotation( condition = anno_text(colnames(cds), rot = 270, 
                                                                                   gp = gpar(fontsize = 4), location = 0.5, 
                                                                                   just = "center"), which = c('column')),
                      left_annotation = NULL,
                      right_annotation = NULL,
                      
                      km = 1,
                      split = NULL,
                      row_km = 1,
                      row_km_repeats = 1,
                      row_split = 2,
                      column_km = 1,
                      column_km_repeats = 1,
                      column_split = factor(annotations$condition),
                      gap = unit(1, "mm"),
                      row_gap = unit(3, "mm"),
                      column_gap = unit(1, "mm"),
                      
                      heatmap_width = unit(1, "npc"),
                      width = NULL,
                      heatmap_height = unit(1, "npc"),
                      height = NULL,
                      
                      show_heatmap_legend = TRUE,
                      heatmap_legend_param = list(title = 'Z score')
                      
  )
  return(heat_obj)
}


subset_rmats_by_coverage <- function(rmats, min_cov = 20){
  is1 <- do.call(rbind, lapply(strsplit(rmats$IJC_SAMPLE_1, ','), as.numeric))
  is2 <- do.call(rbind, lapply(strsplit(rmats$IJC_SAMPLE_2, ','), as.numeric))
  ss1 <- do.call(rbind, lapply(strsplit(rmats$SJC_SAMPLE_1, ','), as.numeric))
  ss2 <- do.call(rbind, lapply(strsplit(rmats$SJC_SAMPLE_2, ','), as.numeric))
  new_rmats <- rmats[(rowMeans(is1) > min_cov | rowMeans(ss1) > min_cov) & (rowMeans(is2) > min_cov | rowMeans(ss2) > min_cov),]
  new_rmats$FDR <- p.adjust(new_rmats$PValue, 'BH')
  new_rmats
}



rmats_read_enrich <- function(res_dir,  org = 'mouse', gs_diff = 0.2, filter_low = T, filter_thresh = 20){
  rmats <- list()
  rmats[['a3ss']] <- read.csv(paste(res_dir, 'A3SS.MATS.JC.txt', sep = '/'), sep = '\t', header = T, row.names = 1)
  rmats[['a5ss']] <- read.csv(paste(res_dir, 'A5SS.MATS.JC.txt', sep = '/'), sep = '\t', header = T, row.names = 1)
  rmats[['se']] <- read.csv(paste(res_dir, 'SE.MATS.JC.txt', sep = '/'), sep = '\t', header = T, row.names = 1)
  rmats[['ri']] <- read.csv(paste(res_dir, 'RI.MATS.JC.txt', sep = '/'), sep = '\t', header = T, row.names = 1)
  rmats[['mxe']] <- read.csv(paste(res_dir, 'MXE.MATS.JC.txt', sep = '/'), sep = '\t', header = T, row.names = 1)
  rmats_gs <- list()
  rmats_f <- list()
  for(n in names(rmats)){
    rmats[[n]]$splice_type <- toupper(n)
    df <- rmats[[n]]
    
    if(filter_low){
      df <- subset_rmats_by_coverage(rmats[[n]], min_cov = filter_thresh)
      
    }
    df$diff <- df$FDR < 0.05 & abs(df$IncLevelDifference) > gs_diff
    rmats_f[[n]] <- df
    
    rmats_gs[[n]] <- enrich_CP(subset(df, FDR < 0.05 & abs(IncLevelDifference) > gs_diff)$GeneID, universe = df$GeneID , organisms = org, classic = T, Msig = c('H','C3-TFT:GTRD'), simple_combine = F, combine = T, full_combine = T)
  }
  return(list(rmats = rmats, rmats_gs = rmats_gs, rmats_f = rmats_f))
}

cp_tree_bar_plot <- function(res, splice_res,n_cat = 50, nclust = 5, alpha = 0.1, by_keys = NULL){
  require(ggtree)
  require(enrichplot)
  res@result <- subset(res@result, qvalue < alpha)
  if(sum(grepl('GO', res@result$ID))){
    res@ontology <- 'BP'
    #print(res@setType)
    res <- clusterProfiler::simplify(res)
  }
  gcolors <- paletteer_d("ggsci::springfield_simpsons")[1:nclust]
  if(!is.null(by_keys)){
    sets <- find_key_gs(res@result, keys = c('(cancer|carcinoma)', 
                                             '(cycle|mitoti|mitosi)',
                                             'HALLMARK', 'repair', 'ribos', 'death', 'RNA',
                                             'histon', 'chromat', 'methy', 'DNA',
                                             'colorec'), key_length = c(8, 7,5,2,2,2, 2, 2, 2, 2, 2, 2))
    res@result <- res@result[sets,]
  }
  
  res@result$Description[nchar(res@result$Description) > 60] <- res@result$ID[nchar(res@result$Description) > 60]
  res@result$Description<- str_wrap(res@result$Description, width = 30)
  
  res@result <- res@result[tapply(res@result$ID, res@result$Description, function(x) x[1]),]
  res <- enrichplot::pairwise_termsim(res)
  #res@result <- res@result[tapply(res@result$ID, res@result$Description, function(x) x[order(res@result[x, 'qvalue'])][1]),]
  res_tree <- enrichplot::treeplot(res, showCategory = min(n_cat, nrow(res@result)), 
                                   geneClusterPanel = 'pie', 
                                   cluster.params = list(method = "ward.D2", n = nclust, color = gcolors, label_words_n = 4,label_format = 25), 
                                   color = NULL, offset_tiplab = 0.8, fontsize = 3)+
    geom_tiplab( offset = 0.8, hjust = 0, 
                 show.legend = FALSE, lineheight = 0.8,
                 align=TRUE, size = 3)+xlim(c(0,20))+scale_size(name = "number of genes",
                                                                range = c(1, 4))+
    guides(size = guide_legend(title='Number of Genes'))
  
  res_tree$data$FDR <- as.numeric(res@result[match(res_tree$data$label, res@result$Description),'qvalue'])
  View(res_tree$data)
  res_tree <- addSmallLegend( res_tree+new_scale_color()+
                                geom_point(aes(x=x,y=y, color = FDR, size = count))+
                                scale_color_viridis()+new_scale_color() , 
                              pointSize = 3, textSize = 10, spaceLegend = 0.9)
  sub_tree_data <- subset(res_tree$data, !is.na(label))
  sub_tree_data$geneID <- res@result[match(sub_tree_data$label, res@result$Description),'geneID']
  res_tree$layers[c(7,8)] <- NULL
  res_tree$layers[c(3,4)] <- NULL
  change_df <- data.frame(do.call(rbind,  lapply(strsplit(sub_tree_data$geneID, '/'), function(x) c(length(intersect(x, splice_res$a3ss$geneSymbol)), 
                                                                                                    length(intersect(x, splice_res$a5ss$geneSymbol)),
                                                                                                    length(intersect(x, splice_res$se$geneSymbol)),
                                                                                                    length(intersect(x, splice_res$ri$geneSymbol)),
                                                                                                    length(intersect(x, splice_res$mxe$geneSymbol))))))
  colnames(change_df) <- c('A3SS', 'A5SS', 'SE', 'IR', 'MXE')
  change_df$Description <- sub_tree_data$label
  #return(change_df)
  #res_tree$data <- cbind(res_tree$data, change_df)
  
  
  bar <- addSmallLegend(change_df %>% melt(id.vars = c('Description'), variable.name = 'SpliceType') %>% 
                          ggplot( aes(fill=SpliceType, y=value, x=Description)) +scale_y_continuous(labels = scales::percent_format(scale = 100))+ coord_flip()+
                          geom_bar(position="fill", stat="identity") +
                          ggsci::scale_fill_lancet() +
                          ggtitle("") +
                          theme_classic() + 
                          theme(axis.title.y = element_blank(), 
                                axis.text.y=element_blank(),
                                axis.ticks.length.y = unit(0.2, "cm"))+
                          ylab('Splice Type Percentage')+
                          guides(fill=guide_legend(title="Splice type")), pointSize = 3, textSize = 10, spaceLegend = 0.9)
  
  bar$data$label <- bar$data$Description
  #res_tree$layers[c(3,4)] <- NULL
  #return(list(ridge = res_ridge, tree=res_tree))
  bar_tree <- bar  %>% insert_left(res_tree, width = 2)
  return(list(tree = res_tree, bar = bar, both = bar_tree))
  #res_ridge  %>% insert_left(res_tree, width = 2)
  
  
}

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

