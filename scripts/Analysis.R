


FILES <- 'D:/gDrive/Colon_cancer/Milciclib/'
human_gtf <- makeTxDbFromGFF('./dataset/Homo_sapiens.GRCh38.111.gtf', 'gtf')



options(ucscChromosomeNames=FALSE)

# Seperate Conditions
crc <- list()
crc$meta <- data.frame(row.names = c('C1', 'C2', 'C3', 'M1', 'M2', 'M3'), condition = c(rep('control', 3), rep('treated', 3)))
crc <- prepareCount(crc, dirc = './dataset/star_salmon/')
crc <- prepareGeneFeatures(crc, './dataset/Homo_sapiens.GRCh38.111.gtf')
crc$old_ct <- read.csv('./A-vs-B/counts/raw_counts.csv', row.names = 1)uop
crc$meta$mito_rate <- colSums(crc$tpm$bio[row.names(subset(crc$features_tpm, seqnames == 'MT')),row.names(crc$meta)])/colSums(crc$tpm$bio[,row.names(crc$meta)])
crc$meta$rRNA_rate <- colSums(crc$tpm$bio[row.names(subset(crc$features_tpm, gene_biotype == 'rRNA')),row.names(crc$meta)])/colSums(crc$tpm$bio[,row.names(crc$meta)])
crc$biotype_tpm <- do.call(rbind, tapply(row.names(crc$features_tpm), crc$features_tpm$gene_biotype, FUN = function(x){colSums(crc$tpm$bio[intersect(row.names(crc$tpm$bio), x),, drop = F])/colSums(crc$tpm$bio)}))
crc$biotype_ct <- do.call(rbind, tapply(row.names(crc$features), crc$features$gene_biotype, FUN = function(x){colSums(crc$ct$bio[intersect(row.names(crc$ct$bio), x),, drop = F])/colSums(crc$ct$bio)}))
crc$fcts <- crc$ct$bio
crc$scts <- crc$salmon$gene$counts
crc$stpm <- crc$salmon$gene$abundance
#crc$ct$bio <- crc$ct$bio[crc$features[row.names(crc$ct$bio),'gene_biotype'] %in% c('protein_coding', 'Mt_rRNA', 'lncrna', 'processed_pseudogene'),]
#crc$salmon$gene$counts <- crc$salmon$gene$counts[crc$features[row.names(crc$salmon$gene$counts),'gene_biotype'] %in% c('protein_coding', 'Mt_rRNA', 'lncrna', 'processed_pseudogene'),]
#crc$old_ct <- crc$old_ct[crc$features[row.names(crc$old_ct),'gene_biotype'] %in% c('protein_coding', 'Mt_rRNA', 'lncRNA', 'processed_pseudogene'),]
ENSEMBL2ALIAS <- crc$features$gene_name
names(ENSEMBL2ALIAS) <- row.names(crc$features)

# Plots

prev_res <- read.csv('./A-vs-B/Differential_expression_analysis_table.csv', row.names = 1, header = T)
prev_norm <- read.csv('./A-vs-B/counts/normalized_counts.csv', row.names = 1, header = T)

# it seems that the previous analysis used a filter of mean counts > 10 for genes, resulting in 15847 genes
# after which size factors were calculated
old_res <- deg_deseq(ct=crc$old_ct[rowMeans(crc$old_ct) > 10,], meta = crc$meta, lfc_thresh = log2(1.5), min_cell_grp = 0, min_cell = 0,min.pct = 0)


ct_res <- deg_deseq(crc, lfc_thresh = log2(1.5), min.mean.expr = 5, min_cell_grp = 2, min_cell = 0, min.pct = 0)
salmon_res <- deg_deseq(ct = crc$salmon$gene$counts, meta = crc$meta, lfc_thresh = log2(1.5),  min.mean.expr = 2, min_cell_grp = 2, min_cell = 0, min.pct = 0)


ct_res5 <- deg_deseq(crc, lfc_thresh = log2(1.5), min.mean.expr = 5, min_cell_grp = 2, min_cell = 0, min.pct = 0)
salmon_res5 <- deg_deseq(ct = crc$salmon$gene$counts, meta = crc$meta, lfc_thresh = log2(1.5),  min.mean.expr = 5, min_cell_grp = 2, min_cell = 0, min.pct = 0)



prev_ora_up <- enrich_CP(row.names(subset(prev_res, padj < 0.05 & log2FoldChange > 0.01)), 'human', universe = row.names(crc$features), Msig = c('H', 'NCG', 'C4-3CA','C2-CP:BIOCARTA', 'C3-TFT:GTRD',"C3-MIR:MIRDB", 'C3-TFT:TFT_Legacy','C4-CGN', 'C4-CM', 'C2-CP_KEGG_MEDICUS'), combine = T, simple_combine = T, full_combine = T)
prev_ora_down <- enrich_CP(row.names(subset(prev_res, padj < 0.05 & log2FoldChange < -0.01)), 'human',universe = row.names(crc$features), Msig = c('H', 'NCG', 'C4-3CA','C2-CP:BIOCARTA', 'C3-TFT:GTRD',"C3-MIR:MIRDB", 'C3-TFT:TFT_Legacy','C4-CGN', 'C4-CM', 'C2-CP_KEGG_MEDICUS'), combine = T, simple_combine = T, full_combine = T)


prev_gse <- gse_CP(logFC = sapply(row.names(prev_res), function(x){prev_res[x,'log2FoldChange']}), organisms = 'human', 
                 Msig = c('H', 'C4-3CA','C2-CP:BIOCARTA', 'C3-TFT:GTRD',"C3-MIR:MIRDB", 'C3-TFT:TFT_Legacy','C4-CGN', 'C4-CM', 'C2-CP_KEGG_MEDICUS', 'NCG'), combine = T, simple_combine = T, classic = T, full_combine = T)

### ORA analysis was not performed, using GSEA instead
#ct_ora_up <- enrich_CP(row.names(subset(ct_res$lfc, padj < 0.05 & log2FoldChange > 1)), 'human', universe = row.names(ct_res$lfc), Msig = c('H'), combine = T, simple_combine = T, full_combine = T)
#ct_ora_down <- enrich_CP(row.names(subset(ct_res$lfc, padj < 0.05 & log2FoldChange < -1)), 'human', universe = row.names(ct_res$lfc),Msig = c('H'), combine = T, simple_combine = T, full_combine = T)


#ct_thresh_ora_up <- enrich_CP(row.names(subset(ct_res$thresh, padj < 0.05 & log2FoldChange > 0)), 'human', 
#                              universe = row.names(subset(ct_res$thresh, !is.na(padj))), Msig = c('H',"C3-TFT:GTRD", 'C3-TFT:TFT_Legacy'), 
#                              combine = T, simple_combine = T, full_combine = T)

#ct_thresh_ora_down <- enrich_CP(row.names(subset(ct_res$thresh, padj < 0.05 & log2FoldChange < 0)), 'human', 
#                                universe = row.names(subset(ct_res$thresh, !is.na(padj))),Msig = c('H',"C3-TFT:GTRD", 'C3-TFT:TFT_Legacy'), 
#                                combine = T, simple_combine = T, full_combine = T)


#ct_thresh_tf_ora_up <- enrich_CP(row.names(subset(ct_res$thresh, padj < 0.05 & log2FoldChange > 0)), 'human', 
 #                                classic = F, universe = row.names(subset(ct_res$thresh, !is.na(padj))), 
 #                                Msig = c("C3-TFT:GTRD", 'C3-TFT:TFT_Legacy'), combine = T, simple_combine = T, 
  #                               full_combine = F)

#ct_thresh_tf_ora_down <- enrich_CP(row.names(subset(ct_res$thresh, padj < 0.05 & log2FoldChange < 0)), 'human',
 #                                  classic = F, universe = row.names(subset(ct_res$thresh, !is.na(padj))),
  #                                 Msig = c("C3-TFT:GTRD", 'C3-TFT:TFT_Legacy'), combine = T, simple_combine = T, 
  #                                 full_combine = F)



#salmon_ora_up <- enrich_CP(row.names(subset(salmon_res$lfc, padj < 0.05 & log2FoldChange > 1)), 'human', universe = row.names(salmon_res$lfc), Msig = c('H'), combine = T, simple_combine = T, full_combine = T)
#salmon_ora_down <- enrich_CP(row.names(subset(salmon_res$lfc, padj < 0.05 & log2FoldChange < -1)), 'human', universe = row.names(salmon_res$lfc),Msig = c('H'), combine = T, simple_combine = T, full_combine = T)


#salmon_thresh_ora_up <- enrich_CP(row.names(subset(salmon_res$thresh, padj < 0.05 & log2FoldChange > 0)), 'human', 
 #                                 universe = row.names(subset(salmon_res$thresh, !is.na(padj))), 
  #                                Msig = c('H'), combine = T, simple_combine = T, full_combine = T)

#salmon_thresh_ora_down <- enrich_CP(row.names(subset(salmon_res$thresh, padj < 0.05 & log2FoldChange < 0)), 'human', 
   #                                 universe = row.names(subset(salmon_res$thresh, !is.na(padj))),Msig = c('H'), 
      #                              combine = T, simple_combine = T, full_combine = T)

#salmon_thresh_tf_ora_up <- enrich_CP(row.names(subset(salmon_res$thresh, padj < 0.05 & log2FoldChange > 0)), 'human',
   #                                  classic = F, universe = row.names(subset(salmon_res$thresh, !is.na(padj))),
       #                              Msig = c("C3-TFT:GTRD", 'C3-TFT:TFT_Legacy'), combine = T, simple_combine = T, 
         #                            full_combine = T)

#salmon_thresh_tf_ora_down <- enrich_CP(row.names(subset(salmon_res$thresh, padj < 0.05 & log2FoldChange < 0)), 'human',
      #                                 classic = F,universe = row.names(subset(salmon_res$thresh, !is.na(padj))),
       #                                Msig = c("C3-TFT:GTRD", 'C3-TFT:TFT_Legacy'), combine = T, simple_combine = T, 
        #                               full_combine = T)



### Using combined sets of Reactome, KEGG, GO:BP, Wikipathways and Msigdb Hallmark, TFT:GTRD and TFT_Legacy
ct_gse <- gse_CP(logFC = subset(ct_res$lfc, !is.na(padj))$log2FoldChange, organisms = 'human', 
                 Msig = c( 'H', 'NCG' ,'C3-TFT:GTRD', 'C3-TFT:TFT_Legacy'), combine = T, simple_combine = T, classic = T, full_combine = T)

ct_tf_gse <- gse_CP(logFC = subset(ct_res$lfc, !is.na(padj))$log2FoldChange, organisms = 'human',
                 Msig = c('C3-TFT:GTRD', 'C3-TFT:TFT_Legacy'), combine = T, simple_combine = F, full_combine = T, 
                 classic = F)



salmon_gse <- gse_CP(logFC = subset(salmon_res$lfc, !is.na(padj))$log2FoldChange, organisms = 'human', 
                 Msig = c('H', 'C3-TFT:GTRD', 'C3-TFT:TFT_Legacy'), combine = T, simple_combine = T, full_combine = T)

salmon_tf_gse <- gse_CP(logFC = subset(salmon_res$lfc, !is.na(padj))$log2FoldChange, organisms = 'human', 
                    Msig = c('C3-TFT:GTRD', 'C3-TFT:TFT_Legacy'), combine = T, simple_combine = F, full_combine = T, 
                    classic = F)



ct_gse5 <- gse_CP(logFC = subset(ct_res5$lfc, !is.na(padj))$log2FoldChange, organisms = 'human', 
                 Msig = c('H', 'C3-TFT:GTRD', 'C3-TFT:TFT_Legacy'), combine = T, simple_combine = T, full_combine = T)



salmon_gse5 <- gse_CP(logFC = subset(salmon_res5$lfc, !is.na(padj))$log2FoldChange, organisms = 'human', 
                     Msig = c('H', 'C3-TFT:GTRD', 'C3-TFT:TFT_Legacy'), combine = T, simple_combine = T, full_combine = T)


salmon_gse_list <- gse_CP(logFC = subset(salmon_res$lfc, !is.na(padj))$log2FoldChange, organisms = 'human', 
                          Msig = c('C3-TFT:GTRD', 'C3-TFT:TFT_Legacy'), combine = T, simple_combine = F, full_combine = T, 
                          classic = F)
s_dgn <- gseDGN(salmon_gse_list)
s_do <- gseDO(salmon_gse_list)
s_ncg <- gseNCG(salmon_gse_list)

ct_ora_down_list<- enrich_CP(row.names(subset(ct_res$thresh, padj < 0.05 & log2FoldChange < 0.05)), 'human', universe = row.names(ct_res$thresh),Msig = c('H'), combine = T, simple_combine = T, full_combine = T)


ct_gse_list <- gse_CP(logFC = subset(ct_res$lfc, !is.na(padj))$log2FoldChange, organisms = 'human', 
                          Msig = c('C3-TFT:GTRD', 'C3-TFT:TFT_Legacy'), combine = T, simple_combine = F, full_combine = T, 
                          classic = F)
c_ncg <- gseNCG(ct_gse_list, pvalueCutoff = 1)

c_ncg <- enrichNCG(row.names(subset(ct_res$thresh, padj < 0.05)), universe = row.names(ct_res$thresh), minGSSize = 10, maxGSSize = 500, readable = T)


ncg <- read.csv('~/../Downloads/NCG_cancerdrivers_annotation_supporting_evidence.tsv', sep ='\t', header = T)

ggsave('./deg/mouse_gsea_ridge.png', height = 8, width =4, units = 'in')
                                                                                                                                                                                                                                              
custom_ridgeplot(mast_gse_rat$combined,  terms = c('GO:0016581', 'GO:0090545', 'R-RNO-73728', 'R-RNO-212300', 'GO:0030527', 
     'rno03010', 'GO:0090092', 'rno00020','rno00280','R-RNO-927802',
    'GO:0004930','rno01200','R-RNO-204005'), top_n = 0)+theme(plot.title = element_text(hjust = 1, size = 12), axis.text.y = element_text(face="bold", color="black", size=8))+ggtitle('GSEA of Rat DEGs')
                                                                                                                                                                                                                                          
ggsave('./deg/rat_gsea_ridge.png', height = 8, width =4, units = 'in')


### FINAL deg comes from threshold at mean expression > 2 reads across all samples and using featureCount reads

final_res <- data.frame(salmon_res5$lfc)
final_res$gene_short_name <- crc$features[row.names(final_res), 'gene_name']
final_res$entrez <- tapply(bitr(final_res$gene_short_name, 'ALIAS', 'ENTREZID', org.Hs.eg.db)$ENTREZID, bitr(final_res$gene_short_name, 'ALIAS', 'ENTREZID', org.Hs.eg.db)$ALIAS, function(x){x[1]})[final_res$gene_short_name]
final_res$lfcShrink <- salmon_res5$lfc$log2FoldChange
final_res$geneID <- row.names(final_res)
volcanoplot(final_res, fc = log2(2), lfcthresh = F, pval_col = 'padj', fc_col = 'lfcShrink', alpha = 0.05)+xlim(c(-10, 10))+
  geom_hline(yintercept=-log10(0.05), linetype="dashed", 
             color = "red", size=0.5)+geom_vline(xintercept=-log2(2), linetype="dashed", color = "black", size=0.5)+
              geom_vline(xintercept=log2(2), linetype="dashed", color = "black", size=0.5)+
  annotate(geom="text", x=-8, y=0.5, label=paste('Downregulated:\n', nrow(subset(final_res, log2FoldChange < -1 & padj < 0.05)), 'Genes'),color="black")+
  annotate(geom="text", x=8, y=0.5, label=paste('Upregulated:\n', nrow(subset(final_res, log2FoldChange > 1 & padj < 0.05)), 'Genes'),color="black")
ggsave('./volcano_plot.png', width = 5, height = 4)

hmap <- geneHeatmap(counts(ct_res$dds), exprs_mat = counts(ct_res$dds), genes = row.names(subset(ct_res$thresh, padj < 0.05)), annotations = crc$meta, scale = T, label_genes = T)
png(file='./hmap.png', res = 300, width = 4, height = 7, units = 'in')
draw(hmap)
dev.off()

enriched_tfs <- row.names(crc$features)[match(sapply(subset(ct_gse$combined@result, grepl('TARGET_', ID) & qvalue < 0.05)$ID, function(x){strsplit(x, '_')[[1]][1]}), crc$features$gene_name)]
enriched_tfs <- na.omit(c(enriched_tfs, row.names(crc$features)[match(unique(sapply(subset(ct_gse$combined@result, 
                                                                            !grepl('TARGET_', ID) & 
                                                                              !grepl('HALLMARK', ID) & 
                                                                              !grepl('UNKNOWN', ID) & 
                                                                              grepl('_', ID) & qvalue < 0.05)$ID,
                                                                     function(x){strsplit(x, '_')[[1]][1]})), crc$features$gene_name)]))





View(data.frame(ct_res$thresh[intersect(row.names(ct_res$thresh), enriched_tfs),]))

#### As of Jan 16th, we use DEGs from MAST, input is CPM produced by HISAT2
#### Only protein-coding, transcribed_pseudogene and lncRNA are used
#### genes must be expressed in 10% of either zygote or oocyte samples
#### No other filtering is done through MAST
#### Genes are recorded in the vectors below




ggplot(percentage_summary, aes(Species, Percentage, fill = cellType)) + geom_bar(stat="identity", color="black",position=position_dodge()) + 
  geom_errorbar(aes(ymin = Percentage - sd_Perc, ymax = Percentage + sd_Perc), color = '#ffa500',width = 0.2, position=position_dodge(.9))+
  theme_classic()+scale_fill_manual(values=DOT_COLOR)+
  annotate("text", x = 1, y = 0.023, label = "", size = 5) +annotate("text", x = 2, y = 0.025, label = "**", size = 5)+ggtitle('Unspliced Read %')+theme(plot.title = element_text(hjust = 0.5, size = 12))
ggsave('percent_gene_unspliced.png', width = 3,height = 2.5)

ggplot(num_summary, aes(Species, Num_of_Genes, fill = cellType)) + geom_bar(stat="identity", color="black",position=position_dodge()) + 
  geom_errorbar(aes(ymin = Num_of_Genes - sd_num, ymax = Num_of_Genes + sd_num), color = '#ffa500',width = 0.2, position=position_dodge(.9))+
  theme_classic()+scale_fill_manual(values=DOT_COLOR)+
  annotate("text", x = 1, y = 2900, label = "*", size = 5) +annotate("text", x = 2, y = 2900, label = "**", size = 5)+ggtitle('# Genes w Nascent Reads')+theme(plot.title = element_text(hjust = 0.5, size = 12))
ggsave('number_gene_unspliced.png', width = 3,height = 2.5)


#mouse_unsplic_full_mast <- mast_diff(ct = mouse_intron$unspliced[mouse_intron$genes,], meta = mouse$meta, control = 'mouseEgg', tpm = F, nbins = 10, min_per_bin = 50, freq = 0.1, max_thres = 6, plot = T, min_cell_grp = 5, min_cell = 8)[[1]]$DESig

#rat_unsplic_full_mast <- mast_diff(ct = rat_intron$unspliced[rat_intron$genes,], meta = rat$meta, control = 'ratEgg', tpm = F, nbins = 10, min_per_bin = 50, freq = 0.1, max_thres = 6, plot = T)[[1]]$DESig











### SUPPA2

## SUPPA2 results

## diffsplice parameters are abs(dPSI) > 0.2 and SUPPA2 pval < 0.05



suppa_crc_iso_res <- read.csv('./dataset/salmon_splice/salmon/suppa/diffsplice/per_isoform/M-C_transcript_diffsplice.dpsi', header = T, row.names = 1, sep = '\t')
suppa_crc_iso_res$gene_short_name <- sapply(row.names(suppa_crc_iso_res), FUN = function(x){strsplit(x,"\\;")[[1]][1]})
suppa_crc_iso_res$tname <- sapply(row.names(suppa_crc_iso_res), FUN = function(x){strsplit(x,"\\;")[[1]][2]})

suppa_crc_iso_res <- suppa_crc_iso_res[suppa_crc_iso_res[,1] != 'NaN',]
suppa_crc_iso_res$diff <- abs(suppa_crc_iso_res[,1]) > 0.2 & suppa_crc_iso_res[,2] < 0.05

suppa_crc_iso_ora <- enrich_CP(subset(suppa_crc_iso_res, diff)$gene_short_name, universe = suppa_crc_iso_res$gene_short_name, organisms = 'human', classic = T,
                                 Msig = c('H','C3-TFT:GTRD', 'C3-TFT:TFT_Legacy', 'NCG'), simple_combine = F, combine = T, full_combine = T)

suppa_mouse_eve_ora <- lapply(unique(suppa_mouse_eve_res$event_type), FUN = function(x){
  sub_df <- subset(suppa_mouse_eve_res, event_type == x)
  diff_g = subset(sub_df, diff)$gene_short_name
  if(length(diff_g) <= 2){
    return(NULL)
  }else{
    enrich_CP(diff_g, universe = sub_df$gene_short_name, organisms = 'mouse')
  }
})

names(suppa_rat_eve_ora) <- unique(suppa_rat_eve_res$event_type)
suppa_rat_eve_ora$all <- enrich_CP(do.call(c, lapply(unique(subset(suppa_rat_eve_res, diff)$gene_short_name), FUN =function(x){strsplit(x, '_and_')[[1]]}))
                                   , universe = do.call(c, lapply(unique(suppa_mouse_eve_res$gene_short_name), FUN =function(x){strsplit(x, '_and_')[[1]]})), 
                                   organisms = 'rat')





### DRIMSEQ/STAGER



crc_dtu <- read_rnasplice_dex_dtu('./dataset/salmon_splice/salmon/dexseq_dtu/results/dexseq/DEXSeqDataSet.M-C.rds', './dataset/salmon_splice/salmon/dexseq_dtu/filter/drimseq/dmDSdata.rds')

crc_dtu$prop <- get_dexseq_prop(crc_dtu)



dexseq_crc <- DEXSeqResults(crc_dtu$dexseq, independentFiltering = T)
dexseq_crc_stageR <- stageR_dexseqRes(dexseq_crc)
#read.csv('./splicing/mouse/suppa/mouse_suppa_drimseq/dexseq_dtu/results/dexseq/DEXSeqResults.MF-MU.tsv', header = T,  sep = '\t')

#Number of transcripts DE causing DTU 
crc_det_bygene <- tapply(subset(dexseq_crc_stageR, gene < 0.05), subset(dexseq_crc_stageR, gene < 0.05)$geneID, FUN = function(x){sum(x$transcript < 0.05)})

stager_crc_ora <- enrich_CP(subset(dexseq_crc_stageR, gene < 0.05)$geneID, universe = dexseq_crc_stageR$geneID, organisms = 'human', classic = T, Msig = c('H', 'C3-TFT:GTRD', 'C3-TFT:TFT_Legacy', "NCG"), alpha = 1, combine = T, simple_combine = T, full_combine = T)

crc_prop_diff <- data.frame(gene = crc_dtu$prop$groupID, transcript=crc_dtu$prop$featureID, diff = prop_diff(crc_dtu$prop[,-c(1,2)], crc_dtu$drimseq@samples$condition))
crc_prop_gene <- tapply(crc_prop_diff, crc_prop_diff$gene, FUN = function(x){sum(x$diff) > 0.5})
crc_prop_gene <- names(crc_prop_gene[crc_prop_gene])
print(length(intersect(crc_prop_gene, subset(dexseq_crc_stageR, gene < 0.05)$geneID)))



# plot dexseq proportions

plotDEXSeqDTU(crc_dtu$prop, 'Foxm1', rat_dtu$drimseq@samples, isProportion = T)





##Rmats results
crc_rmats <- rmats_read_enrich('./dataset/gbam_splice/star_salmon/rmats/M-C/rmats_post/', 'human')
write.csv(crc_rmats$rmats_f$se, './putative_results/rmats_filtered_skipped_exons.csv', quote = F)
write.csv(crc_rmats$rmats_f$a3ss, './putative_results/rmats_filtered_A3SS.csv', quote = F)
write.csv(crc_rmats$rmats_f$a5ss, './putative_results/rmats_filtered_A5SS.csv', quote = F)
write.csv(crc_rmats$rmats_f$ri, './putative_results/rmats_filtered_retained_intron.csv', quote = F)
write.csv(crc_rmats$rmats_f$mxe, './putative_results/rmats_filtered_mixed_exons.csv', quote = F)


full_rmats_res <- do.call(rbind.data.frame, lapply(crc_rmats$rmats_f, FUN = function(x){
 x[,c("geneSymbol", "GeneID","splice_type", "PValue", "FDR","IncLevel1", 'IncLevel2', "IncLevelDifference", "diff",
                            "chr", "strand",  "SJC_SAMPLE_1", "IJC_SAMPLE_1",  "SJC_SAMPLE_2", "IJC_SAMPLE_2")]
}))
colnames(full_rmats_res)[colnames(full_rmats_res) == 'diff'] <- 'Significant'

write.csv(full_rmats_res, './putative_results/rmats_filtered_results.csv', quote = F, row.names = F)

write.csv(subset(full_rmats_res, Significant), './putative_results/rmats_filtered_significant_results.csv', quote = F, row.names = F)


rmats_summary <- do.call(rbind, lapply(crc_rmats$rmats_f, function(x) c(length(unique(subset(x, FDR < 0.05 & abs(IncLevelDifference) > 0.2)$GeneID)), 
                                                       length(unique(x$GeneID))-length(unique(subset(x, FDR < 0.05 & abs(IncLevelDifference) > 0.2)$GeneID)))))
colnames(rmats_summary) <- c('Significant', 'Non-significant')
row.names(rmats_summary) <- c('A3SS', 'A5SS', 'SE', 'IR', 'MXE')
rmats_summary <- data.frame(rmats_summary)
rmats_summary$Splice_Type <- row.names(rmats_summary)
rmats_summary %>% melt(id.vars = c('Splice_Type'), variable.name = 'Significant') %>% 
  group_by(Splice_Type) %>% 
  mutate(perc = round(100*(value/sum(value)), 2)) %>% 
  ggplot( aes(fill=Significant, y=value, x=Splice_Type)) +scale_y_continuous()+ 
  geom_bar(position="stack", stat="identity") +
  paletteer::scale_fill_paletteer_d("wesanderson::Darjeeling1") +
  ggtitle("") +
  theme_classic() + 
  theme(
    legend.position = c(0.1, 1),
    legend.justification = c("left", "top"),
    legend.box.just = "right",
    legend.text=element_text(size=10),
    legend.title = element_text(size = 12),
    legend.margin = margin(6, 10, 6, 6), 
    axis.text.x = element_text(size=12, vjust=-0, color="black"),
    axis.text.y = element_text(size=12, vjust=-0, color="black"),
    axis.title.x = element_text(size=14, vjust=-0, color="black"), 
    axis.title.y = element_text(size=14, vjust=2, color="black"))+
  ylab('Number of Genes')+xlab('Splice Type')+
  guides(fill=guide_legend(title=""))+geom_text(aes(label = paste0(value, ' (', perc,"%)"), y = value), position = position_stack(vjust = 0.5), size = 2.5, color = 'black')
ggsave('./rmats_splicing_results.png', height = 3, width = 5, dpi = 300)


crc_rmats_ora <- enrich_CP(do.call(c, lapply(crc_rmats$rmats_f, FUN = function(x){subset(x, FDR < 0.05 & abs(IncLevelDifference) > 0.2)$GeneID})),
                           universe = do.call(c, lapply(crc_rmats$rmats_f, FUN = function(x){x$GeneID})) , organisms = 'human', classic = T, Msig = c('H','C3-TFT:GTRD'), 
                           simple_combine = F, combine = T, full_combine = T)

write.csv(subset(crc_rmats_ora$combined_full@result, qvalue < 0.1), file = './rmats_ora_significant.csv', quote = F, row.names = F)
tree_plots <- cp_tree_bar_plot(crc_rmats_ora$combined_full, crc_rmats$rmats_f)
tree_plots$both
ggsave('tree_plot.png', tree_plots$both, height = 5, width =6)
#### DAPARS APA
## The 

crc_dapars_sf_st25 <- deg_utr('./combined.txt', crc$ct$bio, c('control', 'treatment'), crc$meta)

crc_dapars_sf_st25$gene_res[crc_dapars_sf_st25$gene_res$APA_dist <= 75,]$diff = F

utr_crc_sf_st25_up <- enrich_CP(row.names(subset(crc_dapars_sf_st25$gene_res, diff & mean.diff > 0.1)), 
                                universe = row.names(crc_dapars_sf_st25$gene_res), organisms = 'human', 
                                combine = T, classic = T, full_combine = T,
                                Msig = c('H', 'C3-TFT:GTRD', 'C3-TFT:TFT_Legacy', 'NCG'))

utr_crc_sf_st25_down <- enrich_CP(row.names(subset(crc_dapars_sf_st25$gene_res, diff & mean.diff < -0.1)),
                                  universe = row.names(crc_dapars_sf_st25$gene_res), organisms = 'human', 
                                  combine = T, classic = T, full_combine = T, 
                                  Msig = c('H', 'C3-TFT:GTRD', 'C3-TFT:TFT_Legacy', 'NCG'))

utr_crc_sf_st25_gs <- enrich_CP(row.names(subset(crc_dapars_sf_st25$gene_res, diff)), 
                                universe = row.names(crc_dapars_sf_st25$gene_res), organisms = 'human', 
                                combine = T, classic = T, full_combine = T,
                                Msig = c('H', 'C3-TFT:GTRD', 'C3-TFT:TFT_Legacy', 'NCG'))





#test_cols[unique(do.call(c, lapply(utr_mouse_sf_st25_down$GO_BP_ora@result$geneID[1:40], FUN = function(x){strsplit(x,'/')[[1]]})))] <- 'red'
#test_alpha[unique(do.call(c, lapply(utr_mouse_sf_st25_down$GO_BP_ora@result$geneID[1:40], FUN = function(x){strsplit(x,'/')[[1]]})))] <- 1

test_cols[intersect(row.names(subset(mouse_dapars_sf_st25$gene_res,mean.diff < -0.2 & fdr < 0.05)), row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > log2(2))))] <- 'red'
test_alpha[intersect(row.names(subset(mouse_dapars_sf_st25$gene_res,mean.diff < -0.2 & fdr < 0.05)), row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > log2(2))))] <- 1

png('short_UTR_vs_Expr_mouse.png', width = 3, height = 3.6, units = 'in', res = 300)
plot(subset(mouse_dapars_sf_st25$gene_res, mean.diff < -0.2 & fdr < 0.05)$mean.diff, 
     mast_res$mouseEgg_v_mouseZygote$DESig[row.names(subset(mouse_dapars_sf_st25$gene_res, mean.diff < -0.2 & fdr < 0.05)),]$Log2FC, 
     pch = 20, col = alpha(test_cols,  test_alpha), ylab='', xlab='')
title(line =0.8, main = expression('Shortened UTR vs Log'[2]*'FC'))
title( line=1.8, cex.lab=1.2, family="Arial", xlab = 'UTR PDUI Diff', ylab = expression('Expr Log'[2]*'FC'))
dev.off()
plot(mouse_dapars_sf_st25$gene_res[intersect(mouse_dapars_sf_st25$gene_res$gene_short_names, row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, abs(Log2FC) > 1 & fdr < 0.05))),]$mean.diff, 
     mast_res$mouseEgg_v_mouseZygote$DESig[intersect(mouse_dapars_sf_st25$gene_res$gene_short_names, row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, abs(Log2FC) > 1 & fdr < 0.05))),]$Log2FC, 
     pch = 20, col = alpha(test_cols,  test_alpha))



dotplot(clusterProfiler::simplify(utr_mouse_sf_st25_gs_stricter$GO_BP_ora, 0.5), showCategory=10, color='pvalue') +theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size = 8, colour =rep('red', 10)), 
                                                                            axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, size = 8), 
                                                                            legend.text = element_text(size=8), legend.key.size = unit(1, 'cm'),  
                                                                            legend.key.height = unit(0.3, 'cm'),legend.key.width = unit(0.3, 'cm'),
                                                                            legend.title = element_text(size=8))+ggtitle('Mouse DAP GO: BP')+
  scale_fill_continuous(low="red", high="blue", name = 'pvalue', guide=guide_colorbar(reverse=TRUE), limits=c(0,0.002))+
  scale_size(range=c(3,8), limits = c(5,220))
ggsave('Mouse_dap_GOBP.png', width = 4, height = 4.5)

cnetplot(clusterProfiler::simplify(utr_mouse_sf_st25_gs$GO_BP_ora, 0.6), 10, foldChange = sapply(row.names(mouse_dapars_sf_st25$gene_res), FUN = function(x){
  mouse_dapars_sf_st25$gene_res[x,]$mean.diff
}), cex_label_gene = 0.3, cex_label_category = 1)




### Plotting Genes Regions From R

test_vp <- plot_range_coverage('NC_000082.7:87,251,316-87,281,194', strand = '+', mouse_gtf, c('./dataset/igv/MF.sf.average.sorted.bw', './dataset/igv/MU.sf.average.sorted.bw'), log = F, cols = DOT_COLOR[c('mouseZygote', 'mouseEgg')], cell_names = c('mouseZygote', 'mouseEgg'), file_name_suffix = 'mouse', gene_name = 'Usp16')
test_vp <- plot_range_coverage('NC_051349.1:35,194,767-35,197,558', strand = '-', rat_gtf, c('./dataset/igv/RF.sf.average.sorted.bw', './dataset/igv/RU.sf.average.sorted.bw'), log = T, cols = DOT_COLOR[c('ratZygote', 'ratEgg')], cell_names = c('ratZygote', 'ratEgg'), file_name_suffix = 'rat', 'Oog1')
test_vp <- plot_range_coverage('NC_000073.7:15,130,624-15,132,520', strand = '+', mouse_gtf, c('./dataset/igv/MF.sf.average.sorted.bw', './dataset/igv/MU.sf.average.sorted.bw'), log = T, cols = DOT_COLOR[c('mouseZygote', 'mouseEgg')], cell_names = c('mouseZygote', 'mouseEgg'), file_name_suffix = 'mouse')


test_vp <- plot_utr_coverage('Cdc25a', utr_res = rat_dapars_sf_st25$gene_res, rat_gtf, c('./dataset/igv/RF.sf.average.sorted.bw', './dataset/igv/RU.sf.average.sorted.bw'), cols= c(DOT_COLOR['ratZygote'], DOT_COLOR['ratEgg']),cell_names = c('Rat Zygote', 'Rat Egg'), file_name_suffix = 'rat')
test_vp <- plot_utr_coverage('Cnot6l', utr_res = rat_dapars_sf_st25$gene_res, rat_gtf, c('./dataset/igv/RF.sf.average.sorted.bw', './dataset/igv/RU.sf.average.sorted.bw'), cols= c(DOT_COLOR['ratZygote'], DOT_COLOR['ratEgg']),cell_names = c('Rat Zygote', 'Rat Egg'), file_name_suffix = 'rat', loci = 'NC_051349.1:13495326-13502511')
test_vp <- plot_utr_coverage('Nek7', utr_res = mouse_dapars_sf_st25$gene_res, mouse_gtf, c('./dataset/igv/MF.sf.average.sorted.bw', './dataset/igv/MU.sf.average.sorted.bw'), cols= c(DOT_COLOR['mouseZygote'], DOT_COLOR['mouseEgg']),cell_names = c('Mouse Zygote', 'Mouse Egg'), file_name_suffix = 'mouse')
test_vp <- plot_utr_coverage('Btg4', utr_res = mouse_dapars_sf_st25$gene_res, mouse_gtf, c('./dataset/igv/MF.sf.average.sorted.bw', './dataset/igv/MU.sf.average.sorted.bw'), cols= c(DOT_COLOR['mouseZygote'], DOT_COLOR['mouseEgg']),cell_names = c('Mouse Zygote', 'Mouse Egg'), file_name_suffix = 'mouse')

