require(GenomicFeatures)
suppressPackageStartupMessages({
  library(ComplexHeatmap)
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
})


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



dnp_deg_mouse <- list('ct Up DEG'= row.names(subset(ct_res$thresh, padj < 0.05 & log2FoldChange > 0)),
                      'ct Down DEG'= row.names(subset(ct_res$thresh, padj < 0.05 & log2FoldChange < 0)),
                      'salmon up'= row.names(subset(salmon_res$thresh, padj < 0.05 & log2FoldChange > 0)),
                      'salmon down' =row.names(subset(salmon_res$thresh, padj < 0.05 & log2FoldChange < 0)))

make_comb_mat(dnp_deg_mouse)
UpSet(make_comb_mat(dnp_deg_mouse)[1:4], comb_col = c('black'))




volcanoplot(mast_res$mouseEgg_v_mouseZygote$DESig, pval_col = 'Pr..Chisq.',fc = log2(2), top_genes = c('Nup37', 'Nup54', 'Obox1', 'Obox2', 'Obox5', 'Obox7', 'Nup35', 'Rpl9', 
                                                                                                       'Rpl3', 'Rpl26', 'Rpl12', 'Rpl19', 'Rpl18', 'Rpl3','Mastl','Cdc20', 
                                                                                                       'Cnot7','Mad2l1', 'Ube2e1','Taf6', 'Taf9','Gtf2b', 'Gtf2a2', 'Ndufc1',
                                                                                                       'Ndufa2', 'Ndufb7','Atp5h','Atp1a1','Ndufa8','Atp5j','Ndufa3', 'Dlat', ))
ggsave('./deg/mouse_volcano_mast.png', width = 4, height = 4)
volcanoplot(mast_rat$ratEgg_v_ratZygote$DESig, pval_col = 'Pr..Chisq.',fc = log2(2), top_genes = c('Rps20', 'Rps2', 'Rps15', 'Rplp1','Hist2h3c2','H2ac1', 'H2ax', 'Gata3',
                                                                                                   'Sdhc','Aco2','Idh1','Sdhd','Pdha1','Dlst','Ndufab1','Ndufc1'))
ggsave('./deg/rat_volcano_mast.png', width = 4, height = 4)

#makeDEGmaps(mouse, res = DEG$mouseEgg_v_mouseZygote, main = 'Mouse Egg vs Zygote', output_dir = './deg/mouse_deg')
#makeDEGmaps(rat, res = DEG$ratEgg_v_ratZygote, main = 'Rat Egg vs Zygote', output_dir = './deg/rat_deg')

#makeDEGmaps(mouse, res = DEG$mouseEgg_v_mouseZygote, main = 'Mouse Egg vs Zygote', output_dir = './deg/mouse_deg', individual = F)
#makeDEGmaps(rat, res = DEG$ratEgg_v_ratZygote, main = 'Rat Egg vs Zygote', output_dir = './deg/rat_deg', individual = F)


#### As of Jan 16th, we use DEGs from MAST, input is CPM produced by HISAT2
#### Only protein-coding, transcribed_pseudogene and lncRNA are used
#### genes must be expressed in 10% of either zygote or oocyte samples
#### No other filtering is done through MAST
#### Genes are recorded in the vectors below


mouse_deg <- mast_res$mouseEgg_v_mouseZygote$DESig
rat_deg <- mast_rat$ratEgg_v_ratZygote$DESig

mouse_deg_up <- row.names(subset(mouse_deg, Log2FC > log2(2) & fdr < 0.05))
mouse_deg_down <- row.names(subset(mouse_deg, Log2FC < -log2(2) & fdr < 0.05))

rat_deg_up <- row.names(subset(rat_deg, Log2FC > log2(2) & fdr < 0.05))
rat_deg_down <- row.names(subset(rat_deg, Log2FC < -log2(2) & fdr < 0.05))




#CP_Plots <- clusterProfilerPlots(DEG_ORA_expr_filt, dir = './enrichment/')

#enrich_plot_grid <- function(plts, method = 'mast'){
  #cowplot::plot_grid(cowplot::plot_grid(plts$pcnet$isotonic_sorbitol_v_AA_starvation$WKP_ora[[method]]+theme(legend.position = "none", plot.title = element_blank()),
                #                        plts$pcnet$isotonic_sorbitol_v_Glucose_starvation$WKP_ora[[method]]+theme(legend.position = "none", plot.title = element_blank()), ncol = 2),
                #  cowplot::plot_grid(plts$pcnet$isotonic_sorbitol_v_AA_starvation$KEGG_ora[[method]]+theme(legend.position = "none", plot.title = element_blank()), 
                #                     plts$pcnet$isotonic_sorbitol_v_Glucose_starvation$KEGG_ora[[method]]+theme(legend.position = "none", plot.title = element_blank()), ncol = 2),
                #  cowplot::plot_grid(plts$pcnet$isotonic_sorbitol_v_AA_starvation$GO_BP_ora[[method]]+theme(legend.position = "none", plot.title = element_blank()), 
                #                   plts$pcnet$isotonic_sorbitol_v_Glucose_starvation$GO_BP_ora[[method]]+theme(legend.position = "none", plot.title = element_blank()), ncol = 2), nrow = 3, labels = c('A','B','C'), hjust = c(-0.2,-0.2, -0.2), label_fontface = 'bold', label_size = 18)
# ggsave(filename = paste('./results/Enrichment/yeast_',method,'_res.png', sep=''), width = 10, height = 15, units = 'in', dpi = 300)
#

perc_rat_genes_unspliced <- colSums(rat_intron$unspliced[rat_intron$genes,])/(colSums(rat_intron$spliced )+colSums(rat_intron$unspliced[rat_intron$genes,]))
num_rat_genes_unspliced <- colSums(rat_intron$unspliced[rat_intron$genes,] > 0)#/(colSums(mouse_intron$spliced >0 ))#+colSums(mouse_intron$unspliced[mouse_intron$genes,]))

perc_mouse_genes_unspliced <- colSums(mouse_intron$unspliced[mouse_intron$genes,])/(colSums(mouse_intron$spliced )+colSums(mouse_intron$unspliced[mouse_intron$genes,]))
num_mouse_genes_unspliced <- colSums(mouse_intron$unspliced[mouse_intron$genes,] > 0)#/(colSums(mouse_intron$spliced >0 ))#+colSums(mouse_intron$unspliced[mouse_intron$genes,]))

percentage_summary <- data.frame(percentages = c(perc_mouse_genes_unspliced, perc_rat_genes_unspliced),
           Species = c(rep('Mouse', 28), rep('Rat', 25)),
           cellType = c(mouse$meta$cellType, rat$meta$cellType)) %>% group_by(cellType) %>% summarise(Percentage = mean(percentages), sd_Perc = sd(percentages), Species = unique(Species))

num_summary <- data.frame(num_of_Genes = c(num_mouse_genes_unspliced, num_rat_genes_unspliced),
                                      Species = c(rep('Mouse', 28), rep('Rat', 25)),
                                      cellType = c(mouse$meta$cellType, rat$meta$cellType)) %>% group_by(cellType) %>% summarise(Num_of_Genes = mean(num_of_Genes), sd_num = sd(num_of_Genes), Species = unique(Species))


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







ggplot(data = unsplic_splic_log2FC, aes(x = x, y = y)) + 
  geom_point()+
  facet_wrap(~Species)+
  theme_classic()+
  ylab(expression('Mature Log'[2]*'FC'))+
  xlab(expression('Nascent Log'[2]*'FC'))+
  geom_text(data = data.frame(x = c(-2, -2), y = c(4,4), lab = paste('PCC:',c(round(cor(mouse_unsplic_mast$Log2FC, mast_res$mouseEgg_v_mouseZygote$DESig[row.names(mouse_unsplic_mast),'Log2FC'],use = 'na.or.complete', method = 'pearson'),3),
                                                              round(cor(rat_unsplic_mast$Log2FC, mast_rat$ratEgg_v_ratZygote$DESig[row.names(rat_unsplic_mast),'Log2FC'], use = 'na.or.complete',method = 'pearson'),3)), sep = ' '),
                                Species = c('Mouse', 'Rat')),mapping = aes(x = x, y = y, label = lab))
ggsave('nascent_mature_log2FC.png', height = 3, width =3.5)


  





dev.off()






sample_PCA(log(t(t(mouse_intron$unspliced[mouse_intron$genes,])/estimateSizeFactorsForMatrix(mouse_intron$spliced))+1), mouse$meta, umap =  T, dimension=2 ,main = "UMAP of Log Expression", labeling = F, point_size = 2, color_by = 'cellType', umap.config = list(n_neighbors = 15, min_dist = 0.5, metric='pearson'))
ggsave(paste(output_dir,"_PCA_label.tiff", sep=''), units="in", width=5, height=4, dpi=300, compression = 'lzw')
sample_PCA(log(t(t(rat_intron$unspliced[rat_intron$genes,])/estimateSizeFactorsForMatrix(rat_intron$spliced))+1), rat$meta, umap =  T, dimension=2 ,main = "UMAP of Log Expression", labeling = F, point_size = 2, color_by = 'cellType', umap.config = list(n_neighbors = 15, min_dist = 0.5, metric='pearson'))
ggsave(paste(output_dir,"_PCA_label.tiff", sep=''), units="in", width=5, height=4, dpi=300, compression = 'lzw')



### SUPPA2
{
## SUPPA2 results

## diffsplice parameters are abs(dPSI) > 0.2 and SUPPA2 pval < 0.05
### AS of Nov 28th, total of 5430 genes have diff splice events in rat data and 6275 genes have splice events in mouse
### total 358 sig diff splice in mouse, 82 genes sig diffsplice in rat
### events by type:
# mouse rat
#SE 220 36 
#MX 6 1
#AL 4 4
#AF 81  20
#RI 11  2
#A5 59  11
#A3 64  12

suppa_mouse_eve_res <- read.csv('./splicing/mouse/suppa/mouse_suppa_drimseq/suppa/diffsplice/per_local_event/MF-MU_local_diffsplice.dpsi', header = T, row.names = 1, sep = '\t')
suppa_mouse_eve_res$gene_short_name <- sapply(row.names(suppa_mouse_eve_res), FUN = function(x){strsplit(x,"\\;")[[1]][1]})
suppa_mouse_eve_res$event_type <- sapply(sapply(row.names(suppa_mouse_eve_res), FUN = function(x){strsplit(x,"\\;")[[1]][2]}), FUN = function(x){strsplit(x,"\\:")[[1]][1]})
suppa_mouse_eve_res <- suppa_mouse_eve_res[suppa_mouse_eve_res[,1] != 'NaN',]
suppa_mouse_eve_res$diff <- abs(suppa_mouse_eve_res[,1]) > 0.2 & suppa_mouse_eve_res[,2] < 0.05

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
names(suppa_mouse_eve_ora) <- unique(suppa_mouse_eve_res$event_type)
suppa_mouse_eve_ora$all <- enrich_CP(do.call(c, lapply(unique(subset(suppa_mouse_eve_res, diff)$gene_short_name), FUN =function(x){strsplit(x, '_and_')[[1]]})), 
                                     universe = do.call(c, lapply(unique(suppa_mouse_eve_res$gene_short_name), FUN =function(x){strsplit(x, '_and_')[[1]]})), 
                                     organisms = 'mouse')


suppa_rat_eve_res <- read.csv('./splicing/rat/suppa/rat_suppa_drimseq/suppa/diffsplice/per_local_event/RF-RU_local_diffsplice.dpsi', header = T, row.names = 1, sep = '\t')
suppa_rat_eve_res$gene_short_name <- sapply(row.names(suppa_rat_eve_res), FUN = function(x){strsplit(x,"\\;")[[1]][1]})
suppa_rat_eve_res$event_type <- sapply(sapply(row.names(suppa_rat_eve_res), FUN = function(x){strsplit(x,"\\;")[[1]][2]}), FUN = function(x){strsplit(x,"\\:")[[1]][1]})
suppa_rat_eve_res <- suppa_rat_eve_res[suppa_rat_eve_res[,1] != 'NaN',]
suppa_rat_eve_res$diff <- abs(suppa_rat_eve_res[,1]) > 0.2 & suppa_rat_eve_res[,2] < 0.05

suppa_rat_iso_res <- read.csv('./splicing/rat/suppa/rat_suppa_drimseq/suppa/diffsplice/per_isoform/RF-RU_transcript_diffsplice.dpsi', header = T, row.names = 1, sep = '\t')
suppa_rat_iso_res$gene_short_name <- sapply(row.names(suppa_rat_iso_res), FUN = function(x){strsplit(x,"\\;")[[1]][1]})
suppa_rat_iso_res <- suppa_rat_eve_res[suppa_rat_iso_res[,1] != 'NaN',]
suppa_rat_iso_res$diff <- abs(suppa_rat_iso_res[,1]) > 0.2 & suppa_rat_iso_res[,2] < 0.05

suppa_rat_iso_ora <- enrich_CP(subset(suppa_rat_iso_res, diff)$gene_short_name, universe = suppa_rat_iso_res$gene_short_name, organisms = 'rat')
suppa_rat_eve_ora <- lapply(unique(suppa_rat_eve_res$event_type), FUN = function(x){
  sub_df <- subset(suppa_rat_eve_res, event_type == x)
  diff_g = subset(sub_df, diff)$gene_short_name
  if(length(diff_g) <= 2){
    return(NULL)
  }else{
    enrich_CP(diff_g, universe = sub_df$gene_short_name, organisms = 'rat')
  }
})

names(suppa_rat_eve_ora) <- unique(suppa_rat_eve_res$event_type)
suppa_rat_eve_ora$all <- enrich_CP(do.call(c, lapply(unique(subset(suppa_rat_eve_res, diff)$gene_short_name), FUN =function(x){strsplit(x, '_and_')[[1]]}))
                                   , universe = do.call(c, lapply(unique(suppa_mouse_eve_res$gene_short_name), FUN =function(x){strsplit(x, '_and_')[[1]]})), 
                                   organisms = 'rat')



#Count genes
length(unique(suppa_mouse_eve_res$gene_short_name))
length(do.call(c, lapply(unique(subset(suppa_mouse_eve_res, local_MF.local_MU_p.val < 0.05 )$gene_short_name), FUN =function(x){strsplit(x, '_and_')[[1]]})))
tapply(row.names(suppa_mouse_eve_res), suppa_mouse_eve_res$event_type, FUN =function(x){length(unique(subset(suppa_mouse_eve_res[x,], diff)$gene_short_name))})

length(unique(suppa_rat_eve_res$gene_short_name))
length(unique(subset(suppa_rat_eve_res, diff)$gene_short_name))
tapply(row.names(suppa_rat_eve_res), suppa_rat_eve_res$event_type, FUN =function(x){length(unique(subset(suppa_rat_eve_res[x,], diff)$gene_short_name))})

}



### DRIMSEQ/STAGER
stager_rat <- read.csv('./dataset/salmon_splice/salmon/dexseq_dtu/results/stager/getAdjustedPValues.M-C.tsv', header = T,  sep = '\t')

#stager_rat_ora1 <- enrich_CP(subset(stager_rat, gene < 0.05)$geneID, universe = stager_rat$geneID, organisms = 'rat')

stager_mouse <- read.csv('./splicing/mouse/suppa/mouse_suppa_drimseq/dexseq_dtu/results/stager/getAdjustedPValues.MF-MU.tsv', header = T,  sep = '\t')

#stager_mouse_ora1 <- enrich_CP(subset(stager_mouse, gdene < 0.05)$geneID, universe = stager_mouse$geneID, organisms = 'mouse')




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
stager_mouse_ora1 <- enrich_CP(intersect(crc_prop_gene, subset(dexseq_crc_stageR, gene < 0.05)$geneID), universe = dexseq_crc_stageR$geneID, organisms = 'human', classic = T, Msig = c('H', 'C3-TFT:GTRD', 'C3-TFT:TFT_Legacy'), alpha = 1, combine = T, simple_combine = T, full_combine = T)



dotplot(stager_rat_ora$GO_BP_ora, showCategory=10, color='pvalue') +theme(axis.text.y = element_text(angle = 0, vjust = 0, hjust=0.5, size = 8), 
                                                                          axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, size = 8), 
                                                                          legend.text = element_text(size=8), legend.key.size = unit(1, 'cm'),  
                                                                          legend.key.height = unit(0.3, 'cm'),legend.key.width = unit(0.3, 'cm'),
                                                                          legend.title = element_text(size=8))+ggtitle('Rat DTU GO: BP')+ scale_fill_continuous(low="red", high="blue", name = 'pvalue', guide=guide_colorbar(reverse=TRUE), limits=c(0,0.002))+
                                                                              scale_size(range=c(3,8), limits = c(5,220))
ggsave('Rat_dex_dtu_GOBP.png', width = 4, height = 4.5)

dotplot(stager_mouse_ora$GO_BP_ora, showCategory=10, color='pvalue') +theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size = 8, colour = c( rep('black',4),rep('red', 6))), 
                                                                    axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, size = 8), 
                                                                    legend.text = element_text(size=8), legend.key.size = unit(1, 'cm'),  
                                                                    legend.key.height = unit(0.3, 'cm'),legend.key.width = unit(0.3, 'cm'),
                                                                    legend.title = element_text(size=8))+ggtitle('Mouse DTU GO: BP')+
                                                                    scale_fill_continuous(low="red", high="blue", name = 'pvalue', guide=guide_colorbar(reverse=TRUE), limits=c(0,0.002))+
                                                                    scale_size(range=c(3,8), limits = c(5,220))
ggsave('Mouse_dex_dtu_GOBP.png', width = 4, height = 4.5)


View(subset(dexseq_mouse_stageR, gene < 0.05))

table(tapply(subset(dexseq_mouse_stageR, gene < 0.05), subset(dexseq_mouse_stageR, gene < 0.05)$geneID, FUN = function(x){sum(x$transcript < 0.05)}))

table(tapply(subset(dexseq_rat_stageR, gene < 0.05), subset(dexseq_rat_stageR, gene < 0.05)$geneID, FUN = function(x){sum(x$transcript < 0.05)}))



dtu_deg_mouse <- list('Mouse Up DEG'= row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 1)),
                'Mouse Down DEG'= row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC < -1)),
                'Mouse DTU'= unique(subset(dexseq_mouse_stageR, gene < 0.05)$geneID))

UpSet(make_comb_mat(dtu_deg_mouse)[1:2], comb_col = c('black'))
UpSet(make_comb_mat(dtu_deg_rat)[1:2], comb_col = c('black'))

"#CD534CFF"
ggvenn(
  dtu_deg_mouse, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 3, text_size = 5, show_percentage = F
)

dtu_deg_mouse <- list('Mouse Up DEG'= row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 1)),
                      'Mouse Down DEG'= row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC < -1)),
                      'Mouse DTU'= unique(subset(dexseq_mouse_stageR, gene < 0.05)$geneID))
dtu_deg_rat <- list('Rat Up DEG'= row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC > 1)),
                      'Rat Down DEG'= row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC < -1)),
                      'Rat DTU'= unique(subset(dexseq_rat_stageR, gene < 0.05)$geneID))

"#CD534CFF"
ggvenn(
  dtu_deg_rat, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 3, text_size = 5, show_percentage = F
)

dtu_deg_mouse_ora_up <- enrich_CP(intersect(row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 1)), unique(subset(dexseq_mouse_stageR, gene < 0.05)$geneID)),
                                'mouse', universe = intersect(row.names(mast_res$mouseEgg_v_mouseZygote$DESig), unique(dexseq_mouse_stageR$geneID)))
dtu_deg_mouse_ora_down <- enrich_CP(intersect(row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC < -1)), unique(subset(dexseq_mouse_stageR, gene < 0.05)$geneID)),
                                  'mouse', universe = intersect(row.names(mast_res$mouseEgg_v_mouseZygote$DESig), unique(dexseq_mouse_stageR$geneID)))

dtu_deg_rat_ora_up <- enrich_CP(row.names(subset(mouse_intron_prop, qvalue < 0.05 & prop_diff > 0)), 'mouse', universe = row.names(mouse_intron_prop))
dtu_deg_rat_ora_down <- enrich_CP(row.names(subset(mouse_intron_prop, qvalue < 0.05 & prop_diff < 0)), 'mouse', universe = row.names(mouse_intron_prop))




plotProportions(mouse_drim$drim, gene_id = 'Anapc5', group_variable = "condition", plot_type = "ribbonplot")+theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, size = 8))
plotProportions(rat_drim$drim, gene_id = 'Msl3', group_variable = "condition", plot_type = "boxplot1")+theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, size = 8))


# plot dexseq proportions

plotDEXSeqDTU(rat_dtu$prop, 'Foxm1', rat_dtu$drimseq@samples, isProportion = T)
plotDEXSeqDTU(mouse_dtu$prop, 'Abi3bp', mouse_drim$samps, isProportion = T)







print(length(unique(subset(stager_mouse, gene < 0.05)$geneID)))
print(length(unique(stager_mouse$geneID)))
print(length(unique(subset(stager_rat, gene < 0.05)$geneID)))
print(length(unique(stager_rat$geneID)))

##Rmats results

subset_rmats_by_coverage <- function(rmats, min_cov = 20){
  is1 <- do.call(rbind, lapply(strsplit(rmats$IJC_SAMPLE_1, ','), as.numeric))
  is2 <- do.call(rbind, lapply(strsplit(rmats$IJC_SAMPLE_2, ','), as.numeric))
  ss1 <- do.call(rbind, lapply(strsplit(rmats$SJC_SAMPLE_1, ','), as.numeric))
  ss2 <- do.call(rbind, lapply(strsplit(rmats$SJC_SAMPLE_2, ','), as.numeric))
  new_rmats <- rmats[(rowMeans(is1) > 20 | rowMeans(ss1) > 20) & (rowMeans(is2) > 20 | rowMeans(ss2) > 20),]
  new_rmats$FDR <- p.adjust(new_rmats$PValue, 'BH')
  new_rmats
}



rmats_read_enrich <- function(res_dir,  org = 'mouse', gs_diff = 0.2, filter_low = T, filter_thresh = 10){
  rmats <- list()
  rmats[['a3ss']] <- read.csv(paste(res_dir, 'A3SS.MATS.JC.txt', sep = '/'), sep = '\t', header = T, row.names = 1)
  rmats[['a5ss']] <- read.csv(paste(res_dir, 'A5SS.MATS.JC.txt', sep = '/'), sep = '\t', header = T, row.names = 1)
  rmats[['se']] <- read.csv(paste(res_dir, 'SE.MATS.JC.txt', sep = '/'), sep = '\t', header = T, row.names = 1)
  rmats[['ri']] <- read.csv(paste(res_dir, 'RI.MATS.JC.txt', sep = '/'), sep = '\t', header = T, row.names = 1)
  rmats[['mxe']] <- read.csv(paste(res_dir, 'MXE.MATS.JC.txt', sep = '/'), sep = '\t', header = T, row.names = 1)
  rmats_gs <- list()
  rmats_f <- list()
  for(n in names(rmats)){
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


crc_rmats <- rmats_read_enrich('./dataset/gbam_splice/star_salmon/rmats/M-C/rmats_post/', 'human')
crc_rmats_f <- rmats_read_enrich('./dataset/gbam_splice/star_salmon/rmats/M-C/rmats_post/', 'human')
crc_rmats_5 <- rmats_read_enrich('./dataset/gbam_splice/star_salmon/rmats/M-C/rmats_post/', 'human', gs_diff = 0.5)
crc_rmats_ora <- enrich_CP(do.call(c, lapply(crc_rmats$rmats, FUN = function(x){subset(x, FDR < 0.05 & abs(IncLevelDifference) > 0.2)$GeneID})),
                           universe = do.call(c, lapply(crc_rmats$rmats, FUN = function(x){x$GeneID})) , organisms = 'human', classic = T, Msig = c('H','C3-TFT:GTRD', 'C3-TFT:TFT_Legacy', 'NCG'), 
                           simple_combine = F, combine = T, full_combine = T)
crc_rmats_ora_5 <- enrich_CP(do.call(c, lapply(crc_rmats$rmats, FUN = function(x){subset(x, FDR < 0.05 & abs(IncLevelDifference) > 0.5)$GeneID})),
                           universe = do.call(c, lapply(crc_rmats$rmats, FUN = function(x){x$GeneID})) , organisms = 'human', classic = T, Msig = c('H','C3-TFT:GTRD', 'C3-TFT:TFT_Legacy', 'NCG'), 
                           simple_combine = F, combine = T, full_combine = T)

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

mouse_dapars_sf_long_300_st25 <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_300_long_st25.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)
mouse_dapars_sf_long_300_st25$gene_res[mouse_dapars_sf_long_300_st25$gene_res$APA_dist <= 75,]$diff = F
utr_mouse_sf_long_300_st25_gs_up <- enrich_CP(unique(subset(mouse_dapars_sf_long_300_st25$gene_res, diff & mean.diff > 0.2)$gene_short_names), universe = unique(mouse_dapars_sf_long_300_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')
utr_mouse_sf_gs_long_300_st25_down <- enrich_CP(unique(subset(mouse_dapars_sf_long_300_st25$gene_res, diff & mean.diff < -0.2)$gene_short_names), universe = unique(mouse_dapars_sf_long_300_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')
utr_mouse_sf_gs_long_300_st25_gs <- enrich_CP(unique(subset(mouse_dapars_sf_long_300_st25$gene_res, diff & abs(mean.diff) >= 0.2)$gene_short_names), universe = unique(mouse_dapars_sf_long_300_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')


{"
mouse_dapars_sf_nodiff <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_short_long_nodiff.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)
mouse_dapars_sf_min_long <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_min_long.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)
mouse_dapars_sf_utr_100 <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_utr_100.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)
mouse_dapars_sf_utr_100_med <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_utr_100_median.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)
mouse_dapars_sf_utr_100_min_long_median <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_utr_100_min_long_median.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)
mouse_dapars_sf_min_long_median <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_min_long_median.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)

mouse_dapars_sf_long_100 <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_100_long.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)
mouse_dapars_sf_long_200 <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_200_long.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)
mouse_dapars_sf_long_300 <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_300_long.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)

mouse_dapars_sf_long_300$gene_res[mouse_dapars_sf_long_300$gene_res$APA_dist == 100,]$diff = F
mouse_dapars_sf_long_200$gene_res[mouse_dapars_sf_long_200$gene_res$APA_dist == 100,]$diff = F
mouse_dapars_sf_long_100$gene_res[mouse_dapars_sf_long_100$gene_res$APA_dist == 100,]$diff = F

mouse_dapars_sf$gene_res[mouse_dapars_sf$gene_res$APA_dist == 100,]$diff = F
utr_mouse_sf_gs_long_100_down <- enrich_CP(unique(subset(mouse_dapars_sf_long_100$gene_res, diff & mean.diff < -0.2)$gene_short_names), universe = unique(mouse_dapars_sf_long_100$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')
utr_mouse_sf_gs_long_200_down <- enrich_CP(unique(subset(mouse_dapars_sf_long_200$gene_res, diff & mean.diff < -0.2)$gene_short_names), universe = unique(mouse_dapars_sf_long_200$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')
utr_mouse_sf_gs_long_300_down <- enrich_CP(unique(subset(mouse_dapars_sf_long_300$gene_res, diff & mean.diff < -0.2)$gene_short_names), universe = unique(mouse_dapars_sf_long_300$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')



mouse_dapars_sf <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)
mouse_dapars_sf$gene_res[mouse_dapars_sf$gene_res$APA_dist == 100,]$diff = F
utr_mouse_sf_gs_up <- enrich_CP(unique(subset(mouse_dapars_sf$gene_res, diff & mean.diff > 0.2)$gene_short_names), universe = unique(mouse_dapars_sf$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')
utr_mouse_sf_gs_down <- enrich_CP(unique(subset(mouse_dapars_sf$gene_res, diff & mean.diff < -0.2)$gene_short_names), universe = unique(mouse_dapars_sf$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')


mouse_dapars_sf_long_100_st25 <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_100_long_st25.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)
mouse_dapars_sf_long_100_st25$gene_res[mouse_dapars_sf_long_100_st25$gene_res$APA_dist <= 75,]$diff = F



utr_mouse_sf_utr_100_gs_up <- enrich_CP(unique(subset(mouse_dapars_sf_utr_100$gene_res, diff & mean.diff > 0.2)$gene_short_names), universe = unique(mouse_dapars_sf_utr_100$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')
utr_mouse_sf_utr_100_gs_down <- enrich_CP(unique(subset(mouse_dapars_sf_utr_100$gene_res, diff & mean.diff < -0.2)$gene_short_names), universe = unique(mouse_dapars_sf_utr_100$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')

utr_mouse_sf_utr_100_min_long_median_gs_up <- enrich_CP(unique(subset(mouse_dapars_sf_utr_100_min_long_median$gene_res, diff & mean.diff > 0.2)$gene_short_names), universe = unique(mouse_dapars_sf_utr_100_min_long_median$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')
utr_mouse_sf_utr_100_min_long_median_gs_down <- enrich_CP(unique(subset(mouse_dapars_sf_utr_100_min_long_median$gene_res, diff & mean.diff < -0.2)$gene_short_names), universe = unique(mouse_dapars_utr_100_min_long_median$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')
utr_mouse_sf_min_long_gs_up <- enrich_CP(unique(subset(mouse_dapars_sf_min_long$gene_res, diff & mean.diff > 0.2)$gene_short_names), universe = unique(mouse_dapars_sf_min_long$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')
utr_mouse_sf_min_long_gs_down <- enrich_CP(unique(subset(mouse_dapars_sf_min_long$gene_res, diff & mean.diff < -0.2)$gene_short_names), universe = unique(mouse_dapars_sf_min_long$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')

  
"}



intersect(subset(mouse_dapars_sf_st25$gene_res, diff)$gene_short_names, subset(rat_dapars_sf_st25$gene_res, diff)$gene_short_names)
intersect(subset(mouse_dapars_sf_long_300_st25$gene_res, diff)$gene_short_names, subset(rat_dapars_sf_long_300_st25$gene_res, diff)$gene_short_names)

crc_dapars_sf_st25 <- enrich_CP(intersect(subset(mouse_dapars_sf_st25$gene_res, fdr < 0.05 & mean.diff < -0.1)$gene_short_names, subset(rat_dapars_sf_st25$gene_res, fdr < 0.05 & mean.diff < -0.1)$gene_short_names),
          universe = intersect(mouse_dapars_sf_st25$gene_res$gene_short_names,rat_dapars_sf_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')

both_utr_gs_up <- enrich_CP(intersect(subset(mouse_dapars_sf_st25$gene_res, fdr < 0.05 & mean.diff > 0.2)$gene_short_names, subset(rat_dapars_sf_st25$gene_res, fdr < 0.05 & mean.diff > 0.2)$gene_short_names),
                              universe = intersect(mouse_dapars_sf_st25$gene_res$gene_short_names,rat_dapars_sf_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')

both_utr_gs <- enrich_CP(intersect(subset(mouse_dapars_sf_st25$gene_res, fdr < 0.05 & abs(mean.diff) > 0.1)$gene_short_names, subset(rat_dapars_sf_st25$gene_res, fdr < 0.05 & abs(mean.diff) > 0.1)$gene_short_names),
                              universe = intersect(mouse_dapars_sf_st25$gene_res$gene_short_names,rat_dapars_sf_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')



utr_deg_mouse <- list('Mouse Up DEG'= row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > log2(2))),
                      'Mouse Down DEG'= row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC < -log2(2))),
                        'Mouse Long UTR'= unique(subset(mouse_PAS_dist$PAS_motif, mean.diff > 0.2 & fdr < 0.05 & num_motif > 0)$gene_short_name),
                      'Mouse Short UTR' =unique(subset(mouse_PAS_dist$PAS_motif, mean.diff < -0.2 & fdr < 0.05 & num_motif > 0)$gene_short_name))

utr_deg_rat <- list('Rat Up DEG'= row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC > 1)),
                      'Rat Down DEG'= row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC < -1)),
                      'Rat Long UTR'= unique(subset(rat_PAS_dist$PAS_motif, mean.diff > 0.2 & fdr < 0.05 & num_motif > 0)$gene_short_name),
                      'Rat Short UTR' =unique(subset(rat_PAS_dist$PAS_motif, mean.diff < -0.2 & fdr < 0.05 & num_motif > 0)$gene_short_name))

png('utr_deg_mouse.png', width = 3, height = 2.6, res = 300, units = 'in')
UpSet(make_comb_mat(utr_deg_mouse)[1:4], comb_col = c('black', 'red', 'black','black'))
dev.off()

png('utr_deg_rat.png', width = 2.6, height = 2.6, res = 300, units = 'in')
UpSet(make_comb_mat(utr_deg_rat)[2:4], comb_col = c('black', 'black','black'))
dev.off()

setdiff(row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 1)), subset(mouse_dapars_sf$gene_res, mean.diff < -0.1)$gene_short_names)


setdiff(row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 1)), subset(mouse_dapars_sf$gene_res, mean.diff < -0.1)$gene_short_names)



mast_go_noployA_up <- enrich_CP(setdiff(row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 1)), 
                                subset(mouse_dapars_sf$gene_res, mean.diff < -0.2 & diff)$gene_short_names),
                                universe = row.names(mast_res$mouseEgg_v_mouseZygote$DESig),
                                logFC = sapply(setdiff(row.names(mast_res$mouseEgg_v_mouseZygote$DESig), subset(mouse_dapars_sf$gene_res, mean.diff < -0.1)$gene_short_names), 
                                FUN =function(x){mast_res$mouseEgg_v_mouseZygote$DESig[x,'Log2FC']}), 
                                organisms = 'mouse',
                                GSE = F)

mast_go_short_utr_up <- enrich_CP(intersect(row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 0.05)), 
                                        subset(mouse_dapars_sf$gene_res, mean.diff < -0.05 & diff)$gene_short_names),
                                universe = intersect(row.names(mast_res$mouseEgg_v_mouseZygote$DESig), 
                                                    mouse_dapars_sf$gene_res$gene_short_names),
                                logFC = sapply(setdiff(row.names(mast_res$mouseEgg_v_mouseZygote$DESig), subset(mouse_dapars_sf$gene_res, mean.diff < -0.1)$gene_short_names), 
                                               FUN =function(x){mast_res$mouseEgg_v_mouseZygote$DESig[x,'Log2FC']}), 
                                organisms = 'mouse',
                                GSE = F)


mast_gse_noployA <- gse_CP(setdiff(row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 1)), 
                                         subset(mouse_dapars_sf$gene_res, mean.diff < -0.05)$gene_short_names),
                                 logFC = sapply(setdiff(row.names(mast_res$mouseEgg_v_mouseZygote$DESig), subset(mouse_dapars_sf$gene_res, mean.diff < -0.05)$gene_short_names), 
                                                FUN =function(x){mast_res$mouseEgg_v_mouseZygote$DESig[x,'Log2FC']}), 
                                 organisms = 'mouse',
                                 GSE = T)

mast_gse_rat <- gse_CP(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & abs(Log2FC) > 1)),
                       logFC = sapply(row.names(mast_rat$ratEgg_v_ratZygote$DESig), FUN =function(x){mast_rat$ratEgg_v_ratZygote$DESig[x,'Log2FC']}), organisms = 'rat', GSE = T)


rat_dapars_sf$pdui <- rat_dapars_sf$df[,grepl('PDUI', colnames(rat_dapars_sf$df))]
colnames(rat_dapars_sf$pdui) <- colnames(rat_dapars_sf$long)
rat_dapars_sf$pdui <- rat_dapars_sf$pdui[row.names(rat_dapars_sf$deg),]
rat_dapars_sf$pdui <- t(apply(rat_dapars_sf$pdui, 1, FUN = function(x){x[is.na(x)] = mean(x, na.rm =T); x}))


mouse_dapars_sf$pdui <- mouse_dapars_sf$df[,grepl('PDUI', colnames(mouse_dapars_sf$df))]
colnames(mouse_dapars_sf$pdui) <- colnames(mouse_dapars_sf$long)
mouse_dapars_sf$pdui <- mouse_dapars_sf$pdui[row.names(mouse_dapars_sf$deg),]
mouse_dapars_sf_st25$pdui <- t(apply(mouse_dapars_sf_st25$pdui, 1, FUN = function(x){x[is.na(x)] = mean(x, na.rm =T); x}))



sample_PCA(log(rat_dapars_sf_st25$pdui_imp+1), rat$meta, umap =  T, dimension=2 ,
           main = "PDUI Umap", labeling = F, point_size = 2, color_by = 'cellType',
           umap.config = list(n_neighbors = 25, min_dist = 0.2, metric='cosine'))+
            theme(axis.text.y = element_text( size = 8), axis.title.x = element_text(size = 8),
          axis.text.x = element_text( size = 8) , axis.title.y = element_text(size = 8))+guides(color=FALSE)
ggsave('rat_umap_DAP.png', units = 'in', dpi = 300, height = 2.5, width = 2.5)
sample_PCA(log(mouse_dapars_sf_st25$pdui_imp+1), mouse$meta, umap =  T, dimension=2 ,
           main = "PDUI Umap", labeling = F, point_size = 2, color_by = 'cellType', 
           umap.config = list(n_neighbors = 28, min_dist = 0.4, metric='cosine'))+
          theme(axis.text.y = element_text( size = 8), axis.title.x = element_text(size = 8),
        axis.text.x = element_text( size = 8) , axis.title.y = element_text(size = 8))+guides(color=FALSE)
ggsave('mouse_umap_DAP.png', units = 'in', dpi = 300, height = 2.5, width = 2.5)


mouse_dapars_sf_st25$gene_res = mouse_dapars_sf_st25$gene_res[unique(subset(mouse_PAS_dist$PAS_motif, mean.diff < -0.2 & fdr < 0.05 & num_motif > 0)$gene_short_name),]

test_cols =  rep('black', length(subset(mouse_dapars_sf_st25$gene_res, mean.diff < -0.2 & fdr < 0.05)$mean.diff))
names(test_cols) <- row.names(subset(mouse_dapars_sf_st25$gene_res, mean.diff < -0.2 & fdr < 0.05))
test_alpha = rep(0.3, length(subset(mouse_dapars_sf_st25$gene_res, mean.diff < -0.2 & fdr < 0.05)$mean.diff))
names(test_alpha) <- row.names(subset(mouse_dapars_sf_st25$gene_res, mean.diff < -0.2 & fdr < 0.05))

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

#### Orthology
utr_ortho <- list('Rat Long UTR'= row.names(subset(rat_dapars_sf_st25$gene_res, mean.diff > 0.2 & fdr < 0.05)),
                      'Rat Short UTR' =row.names(subset(rat_dapars_sf_st25$gene_res, mean.diff < -0.2 & fdr < 0.05)),
                      'Mouse Long UTR'= row.names(subset(mouse_dapars_sf_st25$gene_res, mean.diff > 0.2 & fdr < 0.05)),
                      'Mouse Short UTR' =row.names(subset(m$gene_res,mean.diff < -0.2 & fdr < 0.05)))

utr_ortho <- list('Rat Long UTR'= unique(subset(rat_PAS_dist$PAS_motif, mean.diff > 0.2 & fdr < 0.05 & num_motif > 0)$gene_short_name),
                  'Rat Short UTR' =unique(subset(rat_PAS_dist$PAS_motif, mean.diff < -0.2 & fdr < 0.05 & num_motif > 0)$gene_short_name),
                  'Mouse Long UTR'= unique(subset(mouse_PAS_dist$PAS_motif, mean.diff > 0.2 & fdr < 0.05 & num_motif > 0)$gene_short_name),
                  'Mouse Short UTR' =unique(subset(mouse_PAS_dist$PAS_motif, mean.diff < -0.2 & fdr < 0.05 & num_motif > 0)$gene_short_name))

unique(subset(mouse_PAS_dist$PAS_motif, mean.diff < -0.2 & fdr < 0.05 & num_motif > 0)$gene_short_name)

png('utr_ortho.png', width = 3.6, height = 2.6, res = 300, units = 'in')
UpSet(make_comb_mat(utr_ortho)[2:5], comb_col = c('black'))
dev.off()
length(intersect(row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, abs(Log2FC) > log2(1.25) & fdr < 0.05)), row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, abs(Log2FC) > log2(1.25) & fdr < 0.05))))

deg_ortho <- list('Rat Up DEG'= row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, Log2FC > log2(1.25) & fdr < 0.05)),
                  'Rat Down DEG' =row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, Log2FC < -log2(1.25)  & fdr < 0.05)),
                  'Mouse Up DEG'= row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, Log2FC > log2(1.25) & fdr < 0.05)),
                  'Mouse Down DEG' =row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, Log2FC < -log2(1.25)  & fdr < 0.05)))

mu_ru <- intersect(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, Log2FC > log2(1.25) & fdr < 0.05)), row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, Log2FC > log2(1.25) & fdr < 0.05)))
md_ru <- intersect(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, Log2FC > log2(1.25) & fdr < 0.05)), row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, Log2FC < -log2(1.25) & fdr < 0.05)))
mu_rd <- intersect(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, Log2FC < -log2(1.25) & fdr < 0.05)), row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, Log2FC > log2(1.25) & fdr < 0.05)))
md_rd <- intersect(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, Log2FC < -log2(1.25) & fdr < 0.05)), row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, Log2FC < -log2(1.25) & fdr < 0.05)))
m_r <- intersect(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, abs(Log2FC) > log2(1.25) & fdr < 0.05)), row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, abs(Log2FC) > log2(1.25) & fdr < 0.05)))

m_r_state <- c(rep('Mouse Up / Rat Up', length(mu_ru)),rep('Mouse down / Rat Up', length(md_ru)),rep('Mouse up / Rat down', length(mu_rd)),rep('Mouse down / Rat down', length(md_rd)))

names(m_r_state) <- c(mu_ru, md_ru,mu_rd,md_rd)

mouse_rat_all_deg <- intersect(row.names(mast_rat$ratEgg_v_ratZygote$DESig), row.names(mast_res$mouseEgg_v_mouseZygote$DESig))

deg_ortho_ora <- list('mu_ru' <- enrich_CP(mu_ru,universe = mouse_rat_all_deg, logFC = NULL, organisms = 'mouse'),
                      'mu_rd'<- enrich_CP(mu_rd,universe = mouse_rat_all_deg, logFC = NULL, organisms = 'mouse'),
                      'md_ru' <- enrich_CP(md_ru,universe = mouse_rat_all_deg, logFC = NULL, organisms = 'mouse'),
                      'md_rd'<- enrich_CP(md_rd,universe = mouse_rat_all_deg, logFC = NULL, organisms = 'mouse'),
                      'm_r'<- enrich_CP(md_rd,universe = mouse_rat_all_deg, logFC = NULL, organisms = 'mouse')
)
m_r_ora <- enrich_CP(m_r,universe = mouse_rat_all_deg, logFC = NULL, organisms = 'mouse')
m_r_ora$combined <- m_r_ora$GO_BP_ora
m_r_ora$combined@result <- rbind(m_r_ora$GO_BP_ora@result, m_r_ora$WKP_ora@result, m_r_ora$KEGG_ora@result, m_r_ora$REACT_ora@result)
m_r_ora$combined@geneSets <- c(m_r_ora$GO_BP_ora@geneSets, m_r_ora$REACT_ora@geneSets,m_r_ora$KEGG_ora@geneSets, m_r_ora$WKP_ora@geneSets)
m_r_ora_cnetplots <- custom_cnet_plot(m_r_ora$combined, category = c('R-MMU-1428517', 'R-MMU-927802', 'R-MMU-72172', 'R-MMU-72202', 'R-MMU-69620', 'R-MMU-72312', 'GO:0008380', 'GO:0002181', 'GO:0016574','R-MMU-4551638'),
                         gene_color = sapply(mast_rat$ratEgg_v_ratZygote$DESig$features, FUN = function(x){mast_rat$ratEgg_v_ratZygote$DESig[x,'Log2FC']}), 
                         gene_color2 = sapply(mast_res$mouseEgg_v_mouseZygote$DESig$features, FUN = function(x){mast_res$mouseEgg_v_mouseZygote$DESig[x,'Log2FC']}), 
                         layout = 'fr', color_cat_pval = F)

ggsave(plot = m_r_ora_cnetplots$plot1, 'm_r_ora_mouse.png', width = 5, height = 2)
ggsave(plot = m_r_ora_cnetplots$plot2, 'm_r_ora_rat.png', width = 5, height = 2)



png('deg_ortho.png', width = 3.6, height = 2.6, res = 300, units = 'in')
UpSet(make_comb_mat(deg_ortho)[1:4], comb_col = c('black'))
dev.off()


dnp_ortho <- list('Mouse Increase DNP'= row.names(subset(test, prop_diff > 0 & qvalue < 0.05)),
                    'Mouse Decrease DNP' =row.names(subset(test, prop_diff < 0 & qvalue < 0.05)),
                    'Rat Increase DNP'= row.names(subset(test1, prop_diff > 0 & qvalue < 0.05)),
                    'Rat Decrease DNP' =row.names(subset(test1, prop_diff < 0 & qvalue < 0.05)))
png('dnp_ortho.png', width = 3.6, height = 2.6, res = 300, units = 'in')
UpSet(make_comb_mat(dnp_ortho)[1:4], comb_col = c('black'))
dev.off()

dtu_ortho <- list( 'Rat DTU'= unique(subset(dexseq_rat_stageR, gene < 0.05)$geneID),
                      'Mouse DTU'= unique(subset(dexseq_mouse_stageR, gene < 0.05)$geneID))

dtu_ortho_ora <- enrich_CP(intersect(unique(subset(dexseq_rat_stageR, gene < 0.05)$geneID), unique(subset(dexseq_mouse_stageR, gene < 0.05)$geneID)), 
                                            universe = intersect(dexseq_rat_stageR$geneID, dexseq_mouse_stageR$geneID), organisms = 'mouse')



png('dtu_ortho.png', width = 2.6, height = 2.6, res = 300, units = 'in')

ggvenn(
  dtu_ortho, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 0, text_size = 5, show_percentage = F, auto_scale = T
)

dev.off()


png('utr_deg_rat.png', width = 2.6, height = 2.6, res = 300, units = 'in')
UpSet(make_comb_mat(utr_deg_rat)[1:3], comb_col = c('black', 'black','black'))
dev.off()

ggvenn(
  utr_ortho, 
  fill_color = c("#0073C2FF", "#EFC000FF","#CD534CFF",  "#868686FF"),
  stroke_size = 0.5, set_name_size = 5, text_size = 5, show_percentage = F
)


### Plotting Genes Regions From R

test_vp <- plot_range_coverage('NC_000082.7:87,251,316-87,281,194', strand = '+', mouse_gtf, c('./dataset/igv/MF.sf.average.sorted.bw', './dataset/igv/MU.sf.average.sorted.bw'), log = F, cols = DOT_COLOR[c('mouseZygote', 'mouseEgg')], cell_names = c('mouseZygote', 'mouseEgg'), file_name_suffix = 'mouse', gene_name = 'Usp16')
test_vp <- plot_range_coverage('NC_051349.1:35,194,767-35,197,558', strand = '-', rat_gtf, c('./dataset/igv/RF.sf.average.sorted.bw', './dataset/igv/RU.sf.average.sorted.bw'), log = T, cols = DOT_COLOR[c('ratZygote', 'ratEgg')], cell_names = c('ratZygote', 'ratEgg'), file_name_suffix = 'rat', 'Oog1')
test_vp <- plot_range_coverage('NC_000073.7:15,130,624-15,132,520', strand = '+', mouse_gtf, c('./dataset/igv/MF.sf.average.sorted.bw', './dataset/igv/MU.sf.average.sorted.bw'), log = T, cols = DOT_COLOR[c('mouseZygote', 'mouseEgg')], cell_names = c('mouseZygote', 'mouseEgg'), file_name_suffix = 'mouse')





test_vp <- plot_utr_coverage('Cdc25a', utr_res = rat_dapars_sf_st25$gene_res, rat_gtf, c('./dataset/igv/RF.sf.average.sorted.bw', './dataset/igv/RU.sf.average.sorted.bw'), cols= c(DOT_COLOR['ratZygote'], DOT_COLOR['ratEgg']),cell_names = c('Rat Zygote', 'Rat Egg'), file_name_suffix = 'rat')
test_vp <- plot_utr_coverage('Cnot6l', utr_res = rat_dapars_sf_st25$gene_res, rat_gtf, c('./dataset/igv/RF.sf.average.sorted.bw', './dataset/igv/RU.sf.average.sorted.bw'), cols= c(DOT_COLOR['ratZygote'], DOT_COLOR['ratEgg']),cell_names = c('Rat Zygote', 'Rat Egg'), file_name_suffix = 'rat', loci = 'NC_051349.1:13495326-13502511')
test_vp <- plot_utr_coverage('Nek7', utr_res = mouse_dapars_sf_st25$gene_res, mouse_gtf, c('./dataset/igv/MF.sf.average.sorted.bw', './dataset/igv/MU.sf.average.sorted.bw'), cols= c(DOT_COLOR['mouseZygote'], DOT_COLOR['mouseEgg']),cell_names = c('Mouse Zygote', 'Mouse Egg'), file_name_suffix = 'mouse')
test_vp <- plot_utr_coverage('Btg4', utr_res = mouse_dapars_sf_st25$gene_res, mouse_gtf, c('./dataset/igv/MF.sf.average.sorted.bw', './dataset/igv/MU.sf.average.sorted.bw'), cols= c(DOT_COLOR['mouseZygote'], DOT_COLOR['mouseEgg']),cell_names = c('Mouse Zygote', 'Mouse Egg'), file_name_suffix = 'mouse')




# Comparisons with previous studies on ZGA in mice
## Park, 2013, 2015 DBMTEE
park_data <- do.call(rbind, lapply(1:25,FUN = function(x){
  df <- read.csv(paste('~/Downloads/tables/hclust_result/CLS_norm_',x,'.dat', sep = ''), sep = '\t', skip = 1);
  df$clust <- x
  if(x %in% c(1,2,4,12,13,17,18,19,20)){
    df$pattern <- 'maternal'
  }else if(x %in% c(5,6,7,8,9,10,11,14,22,23)){
    df$pattern <- 'mZGA'
  }else{
    df$pattern <- 'Other'
  }
  df}))

common_genes_park <- intersect(row.names(mast_res$mouseEgg_v_mouseZygote$DESig), row.names(park_data))
park_data_sub <- park_data[common_genes_park,]
mast_data_sub <- mast_res$mouseEgg_v_mouseZygote$DESig[common_genes_park,]

dnp_ortho <- list('Our data mZGA'= row.names(subset(mast_data_sub, Log2FC > 1 )),
                  'Our data maternal' =row.names(subset(mast_data_sub, Log2FC < -1)),
                  'DMTEE Maternal'= row.names(subset(park_data_sub, pattern == 'maternal')),
                  'DMTEE mZGA ' =row.names(subset(park_data_sub, pattern == 'mZGA')))
UpSet(make_comb_mat(dnp_ortho )[1:4], comb_col = c('black'))



######Intergenic and close-to-gene reads from QoRTs
mouse_intergenic <- lapply(list.files('./dataset/mouse_data/QCData/hisat/'), function(x){
  df <- read.csv(paste('./dataset/mouse_data/QCData/hisat/', x, 'QC.geneCounts.detailed.txt.gz', sep ='/'), header = T, sep = '\t', row.names = 1)
  df$far_AMBIG.. <- as.numeric(str_remove_all(df$far_AMBIG..,  "[()]"))
  df
})
names(mouse_intergenic) <- list.files('./dataset/mouse_data/QCData/hisat/')


mouse_qorts_gene <- do.call(cbind, lapply(mouse_intergenic, function(x){x[,'COUNT']}))
mouse_close_ig <- do.call(cbind, lapply(mouse_intergenic, function(x){x[,'nearby']}))
row.names(mouse_close_ig) <- row.names(mouse_intergenic$MF1)
mouse_far_ig <- do.call(cbind, lapply(mouse_intergenic,function(x){x[,'far']}))
row.names(mouse_far_ig) <- row.names(mouse_intergenic$MF1)
mouse_far_ig_amb <- do.call(cbind, lapply(mouse_intergenic,function(x){x[,'far_AMBIG..']}))
mouse_qorts_intron <- do.call(cbind, lapply(mouse_intergenic, function(x){x[,'intronic']}))
row.names(mouse_qorts_intron) <- row.names(mouse_intergenic$MF1)

intersect(row.names(mouse_qorts_intron), row.names(mouse_intron$unspliced))



rat_intergenic <- lapply(list.files('./dataset/rat_data/QCData/hisat/'), function(x){
  df <- read.csv(paste('./dataset/rat_data/QCData/hisat/', x, 'QC.geneCounts.detailed.txt.gz', sep ='/'), header = T, sep = '\t', row.names = 1)
  df$far_AMBIG.. <- as.numeric(str_remove_all(df$far_AMBIG..,  "[()]"))
  df
})
names(rat_intergenic) <- list.files('./dataset/rat_data/QCData/hisat/')

rat_close_ig <- do.call(cbind, lapply(rat_intergenic, function(x){x[,'nearby']}))
rat_far_ig <- do.call(cbind, lapply(rat_intergenic,function(x){x[,'far']}))
rat_far_ig_amb <- do.call(cbind, lapply(rat_intergenic,function(x){x[,'far_AMBIG..']}))
rat_qorts_intron <- do.call(cbind, lapply(rat_intergenic, function(x){x[,'intronic']}))

## Intergenic from mostly 1000bp segments in intergenic regions


rat_hisat_intergenic <- read.csv('./dataset/rat_intergenic.hisat.ct.txt', row.names = 1, header =T, sep = '\t')
colnames(rat_hisat_intergenic) <- sapply(strsplit(colnames(rat_hisat_intergenic), '_'), function(x){x[1]})
mouse_hisat_intergenic <- read.csv('./dataset/mouse_intergenic.hisat.ct.txt', row.names = 1, header =T, sep = '\t')
colnames(mouse_hisat_intergenic) <- sapply(strsplit(colnames(mouse_hisat_intergenic), '_'), function(x){x[1]})


rat_hisat_intergenic1 <- read.csv('./dataset/mouse_data/rat_intergenic_1000.hisat.ct.txt', row.names = 1, header =T, sep = '\t')
colnames(rat_hisat_intergenic1) <- sapply(strsplit(colnames(rat_hisat_intergenic1), '_'), function(x){x[1]})
mouse_intergenic1 <- read_htseq_intergenic(samples = colnames(mouse$ct$bio), cov_thresh = 0)
rat_intergenic1 <- read_htseq_intergenic('./dataset/mouse_data/rat_intergenic_1000.hisat.ct.txt', samples = colnames(rat$ct$bio), 0)


rat_intergenic1$df %>% melt('cellType') %>% ggbarplot( x = "variable", y = "value",color = "cellType",fill='cellType', add.params = list(color = 'black'), palette = 'rickandmorty', add = "mean_se",position = position_dodge(0.8))+ 
  stat_compare_means(aes(group = cellType), method = 't.test', label = "p.signif", label.y = c(5000,4000,8000,8000,3000,20000))+xlab('Location')+ylab('Number of Regions')


rat_intergenic1$df %>% melt('cellType') %>% ggbarplot( x = "variable", y = "value",color = "cellType",fill='cellType', add.params = list(color = 'black'), palette = 'rickandmorty', add = "mean_se",position = position_dodge(0.8))+ 
  stat_compare_means(aes(group = cellType), method = 't.test', label = "p.signif", label.y = c(5000,4000,8000,8000,3000,20000))+xlab('Location')+ylab('Number of Regions')





read_htseq_intergenic <- function(fn='./dataset/mouse_data/mouse_intergenic_1000.hisat.ct.txt', ct, meta, cov_thresh = 5){
  samples = colnames(ct)
  
  intergenic <- read.csv(fn, row.names = 1, header =T, sep = '\t')
  intergenic1 <- intergenic[-c(1:5),]
  intergenic1_u <- intergenic1[, grepl('.u.', colnames(intergenic1))]
  intergenic1_n <- intergenic1[, grepl('int_n', colnames(intergenic1))]
  intergenic1_s <- intergenic1[, grepl('int_s', colnames(intergenic1))]
  colnames(intergenic1_u) <- sapply(strsplit(colnames(intergenic1_n), '_'), function(x){x[1]})
  colnames(intergenic1_n) <- sapply(strsplit(colnames(intergenic1_n), '_'), function(x){x[1]})
  colnames(intergenic1_s) <- sapply(strsplit(colnames(intergenic1_n), '_'), function(x){x[1]})
  intergenic1_u <- intergenic1_u[,samples]
  intergenic1_s <- intergenic1_s[,samples]
  intergenic1_n <- intergenic1_n[,samples]
  u_1k = grepl('_u_0', row.names(intergenic1)) & sapply(strsplit(row.names(intergenic1), '_'), length) == 3
  
  d_1k = grepl('_d_0', row.names(intergenic1)) & sapply(strsplit(row.names(intergenic1), '_'), length) == 3
  
  s_1k = grepl('_0', row.names(intergenic1)) & sapply(strsplit(row.names(intergenic1), '_'), length) > 3
  
  u_2to10k = grepl('_u', row.names(intergenic1)) & !grepl('_0', row.names(intergenic1)) & !grepl('far', row.names(intergenic1)) & sapply(strsplit(row.names(intergenic1), '_'), length) == 3
  d_2to10k = grepl('_d', row.names(intergenic1)) & !grepl('_0', row.names(intergenic1)) & !grepl('_far', row.names(intergenic1)) & sapply(strsplit(row.names(intergenic1), '_'), length) == 3
  
  s_2to10k = !grepl('_0', row.names(intergenic1)) & !grepl('far', row.names(intergenic1)) & sapply(strsplit(row.names(intergenic1), '_'), length) > 3
  
  far = grepl('_far', row.names(intergenic1)) 
  
  num_genes = colSums(ct > 0)
  print(num_genes)
  num_reads = colSums(intergenic1_u)
  print(num_reads)
  
  norm_df =data.frame(row.names = colnames(intergenic1_u), 
                      up1k = colSums(intergenic1_u[u_1k,] > cov_thresh)/num_genes,
                      down1k = colSums(intergenic1_u[d_1k,] > cov_thresh)/num_genes,
                      u_2to10k = colSums(intergenic1_u[u_2to10k,] > cov_thresh)*100/num_reads,
                      d_2to10 = colSums(intergenic1_u[d_2to10k,] > cov_thresh)*100/num_reads,
                      s_1k = colSums(intergenic1_u[s_1k,] > cov_thresh)/num_genes,
                      far = colSums(intergenic1_u[far,] > cov_thresh)*1000/num_reads,
                      extend_u = colSums((intergenic1_n[u_1k,] - intergenic1_s[u_1k,]) > cov_thresh)/num_genes,
                      extend_d = colSums((intergenic1_n[d_1k,] - intergenic1_s[d_1k,]) > cov_thresh)/num_genes,
                      cellType = meta$cellType
  )
  
  df = data.frame(row.names = colnames(intergenic1_u), 
                  up1k = colSums(intergenic1_u[u_1k,] > cov_thresh),
                  down1k = colSums(intergenic1_u[d_1k,] > cov_thresh),
                  u_2to10k = colSums(intergenic1_u[u_2to10k,] > cov_thresh),
                  d_2to10 = colSums(intergenic1_u[d_2to10k,] > cov_thresh),
                  s_1k = colSums(intergenic1_u[s_1k,] > cov_thresh),
                  far = colSums(intergenic1_u[far,] > cov_thresh),
                  extend_u = colSums((intergenic1_n[u_1k,] - intergenic1_s[u_1k,]) > cov_thresh),
                  extend_d = colSums((intergenic1_n[d_1k,] - intergenic1_s[d_1k,]) > cov_thresh),
                  all = colSums(intergenic1_u > cov_thresh),
                  cellType = meta$cellType
                  )

  
  small_df = data.frame(row.names = colnames(intergenic1_u), 
                        far = colSums(intergenic1_u[far,] > cov_thresh),
                        all = colSums(intergenic1_u > cov_thresh),
                        cellType = meta$cellType
  )
  return(list(df = df, norm_df = norm_df, small_df = small_df, raw = intergenic, u=intergenic1_u, s=intergenic1_s, n=intergenic1_n, u_1k = u_1k, d_1k=d_1k, s_1k=s_1k,u_2to10k =u_2to10k,d_2to10k=d_2to10k,s_2to10k=s_2to10k,far=far ))

}


mouse_intergenic1 <- read_htseq_intergenic('./dataset/mouse_data/mouse_intergenic_500.hisat.ct.txt', ct = mouse$ct$bio, meta = mouse$meta, cov_thresh = 1)
rat_intergenic1 <- read_htseq_intergenic('./dataset/mouse_data/rat_intergenic_500.hisat.ct.txt', ct = rat$ct$bio, meta = rat$meta, 1)
mouse_rat_int_df <- rbind(rat_intergenic1$small_df, mouse_intergenic1$small_df)
mouse_rat_int_df$organism <- c(rep('Rat', 25), rep('Mouse', 28))

mouse_rat_int_df %>% melt(c('cellType', 'organism')) %>% ggbarplot( x = "variable", y = "value", color = "cellType",fill='cellType', facet.by = 'organism', add.params = list(color = '#ffa500'), 
                                                                    palette = DOT_COLOR, add = "mean_se",position = position_dodge(0.8))+
  stat_compare_means(aes(group = cellType), method = 't.test', label = "p.signif")+xlab('Location')+ylab('Number of Regions')+theme_classic2()
ggsave('intergenic_regions.png', width = 5,height = 4)


only_up_1000_regions = grepl('_u_0', row.names(mouse_hisat_intergenic1)) & sapply(strsplit(row.names(mouse_hisat_intergenic1), '_'), length) == 3

only_down_1000_regions = grepl('_d_0', row.names(mouse_hisat_intergenic1)) & sapply(strsplit(row.names(mouse_hisat_intergenic1), '_'), length) == 3

shared_1000_regions = grepl('_0', row.names(mouse_hisat_intergenic1)) & sapply(strsplit(row.names(mouse_hisat_intergenic1), '_'), length) > 3

close_up_regions = grepl('_u', row.names(mouse_hisat_intergenic1)) & !grepl('_0', row.names(mouse_hisat_intergenic1)) & !grepl('far', row.names(mouse_hisat_intergenic1)) & sapply(strsplit(row.names(mouse_hisat_intergenic1), '_'), length) == 3
close_down_regions = grepl('_d', row.names(mouse_hisat_intergenic1)) & !grepl('_0', row.names(mouse_hisat_intergenic1)) & !grepl('_far', row.names(mouse_hisat_intergenic1)) & sapply(strsplit(row.names(mouse_hisat_intergenic1), '_'), length) == 3

close_shared_regions = !grepl('_0', row.names(mouse_hisat_intergenic1)) & !grepl('far', row.names(mouse_hisat_intergenic1)) & sapply(strsplit(row.names(mouse_hisat_intergenic1), '_'), length) > 3

far_regions = grepl('_far', row.names(mouse_hisat_intergenic1)) 

sum(only_up_1000_regions)

grepl('int_s', colnames(mouse_hisat_intergenic1)), grepl('int_n', colnames(mouse_hisat_intergenic1)), grepl('.u.', colnames(mouse_hisat_intergenic1))


num_mouse_genes_intergenic <- colSums(mouse_hisat_intergenic1[grepl('_0', row.names(mouse_hisat_intergenic1)),] > 0)[row.names(mouse$meta)]
num_rat_genes_intergenic <- colSums(rat_hisat_intergenic1[grepl('_0', row.names(rat_hisat_intergenic1)),] > 0)[row.names(rat$meta)]






intergenic_num_summary <- data.frame(percentages = c(num_mouse_genes_intergenic, num_rat_genes_intergenic),
                                 Species = c(rep('Mouse', 28), rep('Rat', 25)),
                                 cellType = c(mouse$meta$cellType, rat$meta$cellType)) %>% group_by(cellType) %>% summarise(Num_of_genes = mean(percentages), sd_Perc = sd(percentages), Species = unique(Species))



ggplot(intergenic_num_summary, aes(Species, Num_of_genes, fill = cellType)) + geom_bar(stat="identity", color="black",position=position_dodge()) + 
  geom_errorbar(aes(ymin = Num_of_genes - sd_Perc, ymax = Num_of_genes + sd_Perc), color = '#ffa500',width = 0.2, position=position_dodge(.9))+
  theme_classic()+scale_fill_manual(values=DOT_COLOR)+
  annotate("text", x = 1, y = 0.023, label = "", size = 5) +annotate("text", x = 2, y = 0.025, label = "**", size = 5)+ggtitle('Upstream Intergenic Read Number of Genes')+theme(plot.title = element_text(hjust = 0.5, size = 12))
ggsave('up_intergenic.png', width = 3,height = 2.5)









