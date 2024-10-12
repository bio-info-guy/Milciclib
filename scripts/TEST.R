if(!"rWikiPathways" %in% installed.packages()){
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
  BiocManager::install("rWikiPathways")
}
library(rWikiPathways)

if(!"RCy3" %in% installed.packages()){
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
  BiocManager::install("RCy3")
}
library(RCy3)

# Use install.packages() for the following, if necessary:
library(magrittr)
cytoscapePing()
gbm.pathways <- findPathwaysByText('colorectal') # many pathways returned
human.gbm.pathways <- gbm.pathways %>% 
  dplyr::filter(species == "Homo sapiens")
human.gbm.wpids <- human.gbm.pathways$id
commandsRun(paste0('wikipathways import-as-pathway id=', 'WP4172'))



for(p in RCy3::getNetworkList()){exportImage(p, type='PNG', network = p, resolution = 300, units = 'in', height = 7, width = 8) }



assignSetType <- function(res){}


ct_gse$combined_full@result$SetType <- c( 'Msigdb_C3_TFT_Legacy')
ct_gse$combined_full@result[grepl('hsa', ct_gse$combined_full@result$ID), 'SetType'] <- 'KEGG'
ct_gse$combined_full@result[grepl('HSA', ct_gse$combined_full@result$ID), 'SetType'] <- 'Reactome'
ct_gse$combined_full@result[grepl('WP', ct_gse$combined_full@result$ID), 'SetType'] <- 'WIkipathways'
ct_gse$combined_full@result[grepl('GO:', ct_gse$combined_full@result$ID), 'SetType'] <- 'GO_Biological_Process'
ct_gse$combined_full@result[grepl('HALLMARK', ct_gse$combined_full@result$ID), 'SetType'] <- 'Msigdb_HALLMARK'
ct_gse$combined_full@result[grepl('TARGET_', ct_gse$combined_full@result$ID), 'SetType'] <- 'Msigdb_C3_TFT_GTRD'
ct_gse$combined_full@result[grepl('NCG', ct_gse$combined_full@result$ID), 'SetType'] <- 'Network_Cancer_Genes'

ct_gse$combined_full@result$ID[grepl(regex("_[U0-9]"), ct_gse$combined_full@result$ID)]

ct_gse_nonTF <- ct_gse$combined_full
ct_gse_nonTF@result <- subset(ct_gse_nonTF@result,  !grepl('_UNKNOWN', ID) & !grepl('TARGET_', ID) & !grepl(regex("*_[0-9]*"), ID) | grepl('HALLMARK', ID) | grepl('NCG', ID))


ct_gse_TF <- ct_gse$combined_full
ct_gse_TF@result <- subset(ct_gse_TF@result, grepl('TARGET_', ID) | grepl(regex("_[0-9]"), ID))

write.csv(subset(ct_gse_nonTF@result, qvalue < 0.1)[,c(1,2,3,4,5,6,8,11,12)], './GSEA_enriched_GeneSets.csv', quote = F)
write.csv(subset(ct_gse_TF@result, qvalue < 0.1)[,c(1,2,3,4,5,6,8,11,12)], './GSEA_enriched_MSIGDB_C3_Regulators.csv', quote = F)




find_key_gs <- function(res, keys = NULL, key_length = 5, alpha = 0.1){
  if(is.null(keys)){
    return(NULL)
  }
  if(length(key_length) != length(keys)){
    key_length = rep(key_length[1], length(keys))
  }
  res <- subset(res, qvalue < alpha)
  gs <- unique(do.call(c, lapply(1:length(keys), function(x){
    ids = res$ID[grepl(regex(keys[x]), res$Description, ignore.case = T)]
    ids <- ids[1:min(length(ids), key_length[x])]
  })))
  return(gs)
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
  #res_tree$layers[c(3,4)] <- NULL
  return(res_tree)
  res_ridge <- addSmallLegend(ridgeplot(res, showCategory = min(n_cat, nrow(res@result)))+scale_fill_viridis_c()+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank())+xlim(c(-4, 4))+geom_vline(xintercept=0, linetype="dashed", color = "red")+xlab('Log2FoldChange'), pointSize = 3, textSize = 10, spaceLegend = 0.9)
  
  res_ridge$data$label <- res_ridge$data$category
  #return(list(ridge = res_ridge, tree=res_tree))
  ridge_tree <- res_ridge  %>% insert_left(res_tree, width = 3)
  return(list(tree = res_tree, ridge = test_ridge, both = ridge_tree))
  #res_ridge  %>% insert_left(res_tree, width = 2)
  

  
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



View(crc_rmats_f$rmats_gs$se$combined_full@result)


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

test_se_ora <- crc_rmats_f$rmats_gs$se$combined_full
test_se_ora@result <- subset(test_se_ora@result,  !grepl('_UNKNOWN', ID) & !grepl(regex("*_[0-9]*"), ID) | grepl('TARGET', ID) | grepl('HALLMARK', ID) | grepl('NCG', ID))


all_TFs <- na.omit(crc$features[match(sapply(strsplit(names(ct_gse_TF_target@geneSets)[grepl('TARGET_', names(ct_gse_TF_target@geneSets))], '_'), function(x) x[[1]][1]), crc$features$gene_name),'gene_id'])

tf_de_res <- final_res[intersect(row.names(final_res), all_TFs),]

View(ct_gse_TF@result[paste(subset(tf_de_res, padj < 0.05)$gene_short_name, '_TARGET_GENES', sep = ''),])

addSmallLegend <- function(myPlot, pointSize = 1.5, textSize = 7, spaceLegend = 0.3) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
