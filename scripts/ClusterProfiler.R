#### Slim down GO BP terms prior to enrichment with information content and similarity using GoSemSim
library(msigdbr)
reduce_go_by_level <- function(organism, cutoff = 0.9, minSize=10, maxSize = 500, ont = 'BP', keep_by = min, level_fun = max, measure = 'Rel'){
  if(organism == 'mouse'){
    require(org.Mm.eg.db)
    orgdb = org.Mm.eg.db
  }else if(organism == 'rat'){
    require(org.Rn.eg.db)
    orgdb = org.Rn.eg.db
  }
  go_bp <- clusterProfiler:::get_GO_data(orgdb, ont,'ENTREZID')
  semData = GOSemSim::godata(OrgDb = orgdb, ont=ont)
  go_ids = names(sapply(go_bp$PATHID2EXTID, length)[sapply(go_bp$PATHID2EXTID, length) >= minSize & sapply(go_bp$PATHID2EXTID, length) <= maxSize])
  lvls <- do.call(rbind, lapply(1:19, FUN = function(x){goterms = clusterProfiler:::getGOLevel(ont, x); data.frame(ID=goterms, level=rep(x, length(goterms)))}))
  lvls_max <- do.call(rbind, tapply(data.frame(lvls), lvls[,1], FUN = function(x){
      data.frame(ID=x[1,1], level = level_fun(x[,2]))
    })
    )
    df <- data.frame(row.names = go_ids, ID = go_ids, level = lvls_max[go_ids,2])
  print(dim(df))
  reduce_go <- clusterProfiler:::simplify_internal(df, cutoff = cutoff, measure = measure, ontology = ont, by = 'level', select_fun = max, semData = semData)
  #return(list(reduced=reduce_go))
  go_list <- go_bp$PATHID2EXTID[reduce_go$ID]
  #print(length(go_kegg_react_list))
  go_names <- go_bp$PATHID2NAME[reduce_go$ID]
  #print(print(length(go_kegg_react_names)))
  go_P2G <- data.frame(do.call(rbind, lapply(names(go_list), FUN = function(x){
    cbind(rep(x, length(go_list[[x]])), go_list[[x]])
  })))
  return(list(reduced=reduce_go, T2G = go_P2G, T2N = go_names))
  #print(dim(go_kegg_react_P2G))
}
# GO
# Using clusterProfiler GO

# Kegg
# Using clusterProfiler KEGG

#Mesh
# Using cluster Profiler Mesh


#Misgdb gene sets with combination of Msigdbr and gmt files
get_msig <- function(organism, cat, sub_cat='', gmt_dir=NULL){
  require(clusterProfiler)
  require(msigdbr)
  if(cat == 'C4' & sub_cat == '3CA'){
    return(clusterProfiler::read.gmt(paste(gmt_dir, '3CA.gmt', sep = '/')))
  }
  else if(cat == 'C2' & sub_cat == 'CP_KEGG_MEDICUS'){
    return(clusterProfiler::read.gmt(paste(gmt_dir, 'CP_KEGG_MEDICUS.gmt', sep = '/')))
  }
  else if(cat == 'NCG'){
    ncg <- read.csv('~/../Downloads/NCG_cancerdrivers_annotation_supporting_evidence.tsv', sep ='\t', header = T)
    ncg <- ncg[ncg$cancer_type != '',]
    ncg$cancer_type <- paste('NCG_', ncg$cancer_type, sep = '')
    ncg[,c('cancer_type', 'entrez')]
  }
  else{
    msigdbr(species = organism, category = cat, subcategory = sub_cat) %>% dplyr::select(gs_name, entrez_gene)
  }
}

# get organism stuff for clusterprofiler
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

#double axis plot for clusterprofiler results
ratio_pval_plot <- function(res, top_n = 10, alpha = 0.05, cut_by = 'qvalue'){
  res <- res[res[,'qvalue'] < alpha,]
  res <- res[1:min(nrow(res), top_n),]
  NES = res$NES
  #fsetSize = sapply(strsplit(res[,'BgRatio'], '/'), function(x){as.numeric(x[1])})
  pval <- res$qvalue
  df <- data.frame(row.names=row.names(res), description=res[,'Description'], NES = NES,  logp = -log10(pval), pvalue = pval)
  first <- df$NES
  second <- df$logp
  
  max_first  <- max(first)   # Specify max of first y axis
  max_second <- max(second) # Specify max of second y axis
  min_first  <- min(first)   # Specify min of first y axis
  min_second <- min(second) # Specify min of second y axis
  
  # scale and shift variables calculated based on desired mins and maxes
  scale = (max_second - min_second)/(max_first - min_first)
  scale = median(second)
  print(second)
  print(scale)
  shift =  median(first) - median(second*scale)
  
  # Function to scale secondary axis
  scale_function <- function(x, scale, shift){
    return ((x)*scale)
  }
  
  # Function to scale secondary variable values
  inv_scale_function <- function(x, scale, shift){
    return ((x)/scale)
  }
  print(inv_scale_function(df$logp, scale, shift))
  df <- df[order(abs(df$NES), decreasing = T),]
  df$description <- factor(df$description, levels = df$description)
  pkpd <- ggplot(df, aes(x = description)) +
    geom_bar( aes(y=NES), stat="identity", color = 'blue', fill = 'blue') +
    geom_line(aes(y = inv_scale_function(logp, scale, shift), group = 1, color = c('-log10 FDR')), stat = 'identity') +
    scale_y_continuous(sec.axis = sec_axis(~scale_function(., scale, shift), name="-log10 FDR")) +
    labs(x = "Gene Set", y = "Ratio of GeneSet", color = "") +
    scale_color_manual(values = c("black", "blue"))+theme_classic()
  return(pkpd)
}

ClusterProfilerORAResult <- function(res, cds, mode = c('all', 'positive', 'negative'), method = 'ORA',organism = 'mouse', thresh = log2(2), gsea_rank = c('logfc', 'pval_sign', 'pval_logfc')){
  all.res <- list()
  sigTerms <- list()
  pval_col <- c('padj', 'qval', 'qval', 'fdr', 'fdr', 'adj.P.Val', 'FDR', 'ProbDiffMean')
  names(pval_col) <- c(  'scde', 'monocle', 'monocle_ct', 'mast_tpm', 'mast_cpm', 'limma', 'edgeR', 'basics')
  if(method == 'GSE'){
    mode = 'all'
  }
  for(comp in sort(names(res))){
    if(grepl('mouse', comp)){
      require(org.Mm.eg.db)
      organism = 'mouse'}
    else if(grepl('rat', comp)){
      require(org.Rn.eg.db)
        organism = 'rat'}
    
    test <- res[[comp]]
    geneFC <- test$Upregulated
    res0 <- list()
    
    for( x in names(test[['gene_names']])){
      cat('\n', x)
      genes.ORA <- test[['gene_names']][[x]]
      if(x %in% c('consensus', 'innerQuantile')){
        universe <- unique(do.call(c, lapply(test, FUN = function(x){
          if(class(x) == 'data.frame'){
            return(row.names(x))
          }else{
            return(c())
          }})))
      }else{
        universe <- row.names(test[[x]])
      }
      if(length(genes.ORA) < 2){next}
      if(mode == 'all'){
        genes.ORA <- test[['gene_names']][[x]]
        
        
        if(!(x %in% c('consensus', 'innerQuantile', 'wilcox'))){
          
          foldChange <- test[[x]]$Log2FC
          names(foldChange) <- row.names(test[[x]])
        }else{
          
          foldChange <- test[['tpmLog2FC']]}
      }else if(mode == 'positive'){
          if(!x %in% c('consensus', 'innerQuantile', 'wilcox')){
            genes.ORA <- row.names(subset(test[[x]][genes.ORA,], Log2FC > abs(thresh)))
            foldChange <- subset(test[[x]][genes.ORA,], Log2FC > abs(thresh))$Log2FC
            names(foldChange) <- genes.ORA
          }else{
            foldChange <- test[['tpmLog2FC']]
            deg_fc <- foldChange[genes.ORA]
            genes.ORA <- names(deg_fc[deg_fc > abs(thresh)])
            
          }
      }else{
          if(!x %in% c('consensus', 'innerQuantile', 'wilcox')){
            genes.ORA <- row.names(subset(test[[x]][genes.ORA,], Log2FC < -abs(thresh)))
            foldChange <- subset(test[[x]][genes.ORA,], Log2FC < -abs(thresh))$Log2FC
            names(foldChange) <- genes.ORA
          }else{
            foldChange <- test[['tpmLog2FC']]
            deg_fc <- foldChange[genes.ORA]
            genes.ORA <- names(deg_fc[deg_fc < -abs(thresh)])
          }
      }
      if(method == 'ORA'){
          #print(length(genes.ORA))
          #print(genes.ORA[1:10])
          #print(length(foldChange))
          r <- enrich_CP(genes.ORA, logFC = foldChange, universe = universe, organisms = organism)
        
      }else{
        
        if(x %in% c(  'scde', 'monocle', 'monocle_ct', 'mast_tpm', 'mast_cpm', 'limma', 'edgeR', 'basics'))
        {
          if(gsea_rank == 'pval_logfc'){
            foldChange <- -log10(test[[x]][,pval_col[x]]+1e-32)*foldChange
          }else if(gsea_rank == 'pval_sign'){
            foldChange <- -log10(test[[x]][,pval_col[x]]+1e-32)*sign(foldChange)
          }
          r <- gse_CP(genes.ORA, logFC = foldChange, GSE = T, organisms = organism)

        }else{
          r <- NULL
        }
      }
      if(!is.null(r)){
      for(n in names(r)){
        cat('\t',n)
        if(n == "gfc" | n == "gfc0"){next}
        else{
          if(nrow(r[[n]]@result) == 0){
            r[[n]] <- NULL
          }else{
            if(length(subset(r[[n]]@result, qvalue < 0.05)$Description) > 0){
              
              sigTerms[[comp]][[n]][[x]] <- subset(r[[n]]@result, qvalue < 0.05)[,'Description',drop=F]
              sigTerms[[comp]][[n]][[x]] <- sigTerms[[comp]][[n]][[x]][,1]
            }
          }
        }}
      all.res[[x]][[comp]] <- r
      }
    }
    #names(all.res) <- names(test[['gene_names']])
  }
  all.res[['sigTerms']] <- sigTerms
  return(all.res)
}

custom_ridgeplot <- function(gse, terms = NULL, top_n = 5){
  test <- gse
  result <- NULL
  if(top_n > 0){
    result <- rbind(subset(test@result, qvalue < 0.05 & NES > 0)[1:top_n,], subset(test@result, qvalue < 0.05 & NES < 0)[1:top_n,])
  }
  if(!is.null(terms)){
    test@result <- rbind(result, test@result[terms,])
  }
  test@result$p.adjust <- test@result$qvalue
  ridgeplot(test)
}

clusterProfilerPlots <- function(deg_ora, alpha = 0.05, dir = './results/Enrichment/', title = '', color_by = 'cellType'){
  theme0 <-  theme(plot.title = element_text(hjust = 0.5, size=20, face = "bold"), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), axis.text.y = element_text(size=12, face = 'bold'),
                   axis.text.x = element_text(size=15, face = 'bold'), axis.title.x  = element_text(size=20, face = 'bold'))
    RESULT_N <- c("Wikipath ORA", "Wikipath GSEA", "Biocyc ORA", "Biocyc GSEA", "BP GO ORA", "CC GO ORA", "MF GO ORA", "BP GO GSEA", "CC GO GSEA", "MF GO GSEA", 'KEGG ORA','KEGG GSEA','MKEGG ORA','MKEGG GSEA', 'REACT ORA', 'REACT GSEA','MESH ORA', 'MESH GSEA' )
  names(RESULT_N) <- c("WKP_ora", "WKP_gse", "BIO_ora", "BIO_gse", "GO_BP_ora", "GO_CC_ora", "GO_MF_ora", "GO_BP_gse", "GO_CC_gse", "GO_MF_gse", 'KEGG_ora','KEGG_gse','MKEGG_ora','MKEGG_gse', 'REACT_ora', 'REACT_gse','MESH_ora', 'MESH_gse' )
  cnetplot_list <- list()
  barplot_list <- list()
  upsetplot_list <- list()
  condition_colors <- list(mouseEgg = c(colorRampPalette(brewer.pal(9, 'Blues'))(32)[1:16],colorRampPalette(brewer.pal(9, 'Blues'))(128)[64:128]) , 
                           mouseZygote = c(colorRampPalette(brewer.pal(9, 'Reds'))(32)[1:16],colorRampPalette(brewer.pal(9, 'Reds'))(128)[64:128]), 
                           ratEgg = c(colorRampPalette(brewer.pal(9, 'Blues'))(32)[1:16],colorRampPalette(brewer.pal(9, 'Blues'))(128)[64:128]), 
                           ratZygote = c(colorRampPalette(brewer.pal(9, 'Reds'))(32)[1:16],colorRampPalette(brewer.pal(9, 'Reds'))(128)[64:128]))
  for(meth in names(deg_ora)){
    if(meth == 'sigTerms'){next}
    res <- deg_ora[[meth]]
    for(comp in names(res)){
      res_Name <- paste(strsplit(comp, split = '_')[[1]], collapse = " ")
      r <- res[[comp]]
      color2use <- strsplit(comp, split = '_v_')[[1]]
      if(color_by == 'cellType')
      {
      pal <- c(rev(condition_colors[[color2use[1]]]), condition_colors[[color2use[2]]])
      }else{
        pal <- c(colorRampPalette(brewer.pal(9, 'Blues'))(32)[1:16],colorRampPalette(brewer.pal(9, 'Blues'))(128)[64:128])
      }
      for(n in names(r)){
        cat(meth, comp, n, '\n')
        if(n %in% c('gfc', 'gfc0')){next}
         gfc <- ifelse(grepl('KEGG', n), list(res[[comp]][['gfc']]), list(res[[comp]][['gfc']]))[[1]]
         enrich_res <- r[[n]]
         enrich_res@result$Description[is.na(enrich_res@result$Description)] <- enrich_res@result$ID
         if(grepl('REACT', n)){
          enrich_res@result$Description <- sapply(enrich_res@result$Description, function(x){gsub("Nonsense Mediated Decay", 'NMD', gsub("\\s*\\([^\\)]+\\)","",x))})
          enrich_res@result$Description <- sapply(enrich_res@result$Description, function(x){str_wrap(gsub("Exon Junction Complex", 'EJC', gsub("\\s*\\([^\\)]+\\)","",x)), width = 50)})
         }
         enrich_res@result <- subset(enrich_res@result, !grepl('disease', Description, ignore.case = T) & !grepl('pathy', Description))
         enrich_res@result <- enrich_res@result[which(as.numeric(sapply(enrich_res@result$BgRatio, function(x){strsplit(x, '/')[[1]][1]})) < 250),]
         enrich_res@result$qvalue <- qvalue_truncp(enrich_res@result$pvalue)$qvalue
         temp <- subset(enrich_res@result, qvalue < 0.05)
         if(nrow(temp) < 1){
           enrich_res@result <- subset(enrich_res@result, qvalue < 0.1)
           enrich_res@qvalueCutoff <- 0.1
         }else{
         enrich_res@result <- temp
         }
         enrich_res@pvalueCutoff <- 1
         if(grepl('GO', n)){
           temp <- clusterProfiler::simplify(enrich_res, by = 'p.adjust')
           if(nrow(temp@result) < 1){
             enrich_res@result <- subset(enrich_res@result, qvalue < 0.1)
             enrich_res@qvalueCutoff <- 0.1
           }else{
             enrich_res <- temp
           }
         }
         qvals <- enrich_res@result$qvalue
         padjs <- enrich_res@result$p.adjust
         if(nrow(enrich_res@result) == 0 | sum(qvals[!is.na(qvals)] < 0.1) < 1){
            next
         }else{
           node_label <- ifelse(sum(enrich_res@result$Count[1:min(10, nrow(enrich_res@result))]) > 100, 'category', 'all')
           cex_label_category = ifelse(node_label == 'category', 2, 1)
            pcnet <- suppressMessages(suppressWarnings(cnetplot(enrich_res, showCategory = 10, foldChange = gfc, layout = 'graphopt', node_label = 'none', cex_label_category = 2)+
                                                       ggtitle(paste(title, RESULT_N[n], sep = ': ' ))+
                                                       theme(plot.title = element_text(hjust = 0.5, size=20, face = "bold"))+
                                                       scale_colour_gradientn(name = expression('Log'['2']*'FC'), limits= c(-13, 13), colours = pal, na.value = 'khaki')))
            if(length(unique(pcnet$data[!is.na(pcnet$data$color),]$name)) < 2){
              print(unique(pcnet$data[!is.na(pcnet$data$color),]$name))
              next
            }
            pcnet$data[!is.na(pcnet$data$color),]$size <- 1
            cat_df <- subset(pcnet$data, is.na(color))
            cat_df$qval <- enrich_res@result[match(cat_df$name, enrich_res@result$Description),]$qvalue
            #View(cat_df)
            pcnet <- pcnet + 
                   ggnewscale::new_scale_color() +
                   ggraph::geom_node_point(aes_(color=~qval, size=~size), data = cat_df) + 
                   scale_size(limits = c(1,150), breaks = c(10,20,40,80,140), range = c(3,10) ) +
                   scale_colour_gradientn(name = "qvalue", colours = colorRampPalette(rev(brewer.pal(9,  'Purples')))(255)[0:200], limits= c(1e-6, 0.1), breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1), trans = 'log10', oob = scales::oob_squish)+
                   theme(plot.margin=unit(c(0,0,0,0),"mm"), aspect.ratio = 1)
            
          if(node_label == 'category'){
            pcnet <- pcnet+ ggraph::geom_node_text(aes_(label=~name), data = pcnet$data[!grepl('^RP[ABCFLNOPST]', pcnet$data$name) & !grepl('MRP', pcnet$data$name) & !grepl('RLP', pcnet$data$name) & !is.na(pcnet$data$color),],
                                                   size = 2, bg.color = "white", repel=TRUE)
           pcnet <- pcnet+ ggraph::geom_node_text(aes_(label=~name), data = pcnet$data[is.na(pcnet$data$color),],
                                         size = 5, bg.color = "white", repel=TRUE)
          }else if(node_label == 'all'){
            if(nrow(pcnet$data[!is.na(pcnet$data$color),]) > 50){
              pcnet <- pcnet+ ggraph::geom_node_text(aes_(label=~name), data = pcnet$data[!grepl('^RP[ABCFLNOPST]', pcnet$data$name) & !grepl('MRP', pcnet$data$name) & !grepl('RLP', pcnet$data$name) & !is.na(pcnet$data$color),],
                                                     size = 2, bg.color = "white", repel=TRUE)
            }else{
           pcnet <- pcnet+ 
                    ggraph::geom_node_text(aes_(label=~name), data = pcnet$data[!is.na(pcnet$data$color),],
                                                   size = 2.5, bg.color = "white", repel=TRUE)
           }
            pcnet <- pcnet+ 
                     ggraph::geom_node_text(aes_(label=~name), data = pcnet$data[is.na(pcnet$data$color),],
                                                   size = 5, bg.color = "white", repel=TRUE)
          }else{
            pcnet <- pcnet
          }
          
            ggsave(paste(comp,n,meth ,'cnet.png', sep = '.'), path = dir, plot = pcnet, width = 6, height = 6, dpi = 300)
          if(grepl('gse', n)){
            pbar <- barplot(enrich_res, showCategory=10, label_format = 5, legend.text = F)+
              scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) + 
              ggtitle(paste(meth, RESULT_N[n], sep = ': ' ))+
              theme0+
              ggplot2::ylab(res_Name)
            pbar <- suppressMessages(suppressWarnings(pbar+scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(9,  'Purples')))(255)[0:200], name = 'p.adjust',
                                                 guide=guide_colorbar(reverse=F ))))
            ggsave(paste(comp, n,meth ,'bar.png', sep = '.'), path = dir, plot = pbar,width = 6, height = 6, dpi = 300)
            
            barplot_list[[comp]][[n]][[meth]] <- pbar
            #cowplot::plot_grid(pcnet , pbar, ncol=2, nrow = 1, labels=NULL, rel_width = c(1.3, 1.5))
          }#else{
            #cowplot::plot_grid(pcnet, ncol=1, labels=NULL, rel_widths=c(.9, 0.9))
          #}
          cnetplot_list[[comp]][[n]][[meth]] <- pcnet
          #ggsave(paste('./results/Enrichment/', comp, '.',n,'.',meth ,'.png', sep = ''), width = 20, height = 8, dpi = 300)
        }
      }
    }
  }
  
  #multiplotClusterProfiler(cnetplot_list)
  #multiplotClusterProfiler(barplot_list, type = 'bar')
  return(list(pcnet = cnetplot_list, pbar = barplot_list))
}







multiplotClusterProfiler <- function(res_plot_list, type = 'cnet', dir = './results/Enrichment/'){
for(comp in names(res_plot_list)){
  r <- res_plot_list[[comp]]
  for(R in names(r)){
    plots0 <- r[[R]]
    for(i in 1:(length(plots0)-1)){plots0[[i]] <- plots0[[i]]+theme(legend.position = "none", title = element_blank())}
    plot_g <- cowplot::plot_grid(plotlist = plots0, ncol=min(3, ceiling(sqrt(length(plots0)))), nrow = ceiling(length(plots0)/min(3, ceiling(sqrt(length(plots0))))), labels=NULL, label_size = 22)
    title <- cowplot::ggdraw() + 
      cowplot::draw_label(
        paste(gsub('_', ' ', comp), ': ', RESULT_N[R]),
        fontface = 'bold', size = 30, 
        x = 0,
        hjust = 0
      ) 
      
    plot1 <- cowplot::plot_grid(
      title, plot_g,
      ncol = 1,
      # rel_heights values control vertical title margins
      rel_heights = c(0.05, 1)
    )+theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
    ggsave(paste(comp, R , type, 'png', sep = '.'), plot = plot1, path = dir, width = min(3, ceiling(sqrt(length(plots0))))*6, height = ceiling(length(plots0)/min(3, ceiling(sqrt(length(plots0)))))*6, dpi = 300)
    }
  }
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

 #'graphopt' 'fr' 'kk' 'drl'  'lgl'



extract_geneSets <- function(x, n) {
  n <- update_n(x, n)
  
  if (inherits(x, 'list')) {
    geneSets <- x
  } else {
    geneSets <- geneInCategory(x) ## use core gene for gsea result
    y <- as.data.frame(x@result)
    geneSets <- geneSets[y$ID]
    names(geneSets) <- y$Description        
  }
  #print(geneSets)
  if (is.numeric(n)) {
    print(n)
    return(geneSets)
    return(geneSets[1:n])
  }
  return(geneSets[n]) ## if n is a vector of Description
}

update_n <- function(x, showCategory) {
  if (!is.numeric(showCategory)) {
    if (inherits(x, 'list')) {
      showCategory <- showCategory[showCategory %in% names(x)]
    }
    return(showCategory)
  }
  
  ## geneSets <- geneInCategory(x) ## use core gene for gsea result
  n <- showCategory
  if (inherits(x, 'list')) {
    nn <- length(x)
  } else {
    nn <- nrow(x@result)
  }
  if (nn < n) {
    n <- nn
  }
  
  return(n)
}

setReadable0 <- function (x, gene2symbol, keyType = "auto") {
  if(is.null(x)){return(x)}
  if (!(is(x, "enrichResult") || is(x, "groupGOResult") || 
        is(x, "gseaResult"))) 
    stop("input should be an 'enrichResult' or 'gseaResult' object...")
  isGSEA <- FALSE
  if (is(x, "gseaResult")) 
    isGSEA <- TRUE
  if (keyType == "auto") {
    keyType <- x@keytype
    if (keyType == "UNKNOWN") {
      stop("can't determine keyType automatically; need to set 'keyType' explicitly...")
    }
  }
  if (x@readable) 
    return(x)
  gc <- geneInCategory(x)
  if (isGSEA) {
    genes <- names(x@geneList)
  }
  else {
    genes <- x@gene
  }
  gn <- gene2symbol
  gc <- lapply(gc, function(i) gn[i])
  res <- x@result
  gc <- gc[as.character(res$ID)]
  geneID <- sapply(gc, paste0, collapse = "/")
  if (isGSEA) {
    res$core_enrichment <- unlist(geneID)
  }
  else {
    res$geneID <- unlist(geneID)
  }
  x@gene2Symbol <- gn
  x@result <- res
  x@keytype <- keyType
  x@readable <- TRUE
  return(x)
}

list2graph <- function(inputList) {
  x <- list2df(inputList)
  g <- graph.data.frame(x, directed=FALSE)
  return(g)
}








