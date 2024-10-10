#Making Plots and Analysis of figures  in Chen's Dissertation Paper
umap_config= list(n_neighbors = 15, min_dist = 0.5, metric='pearson')
sample_PCA <- function(cds,
           meta,
           color_by = 'cellType',
           umap = F,
           labeling = FALSE,
           point_size = 4,
           dimension = 2,
           reduce_noise = FALSE,
           return_matrix = FALSE,
           highlight = NULL,
           main = 'plot',
           umap.config = umap_config)
  {
    if (!umap) {
      test <- rsvd::rpca(t(cds), center = T, scale = T)
      percentVar <-
        c(
          100 * round(test$sdev[1] ^ 2 / sum(test$sdev ^ 2), 3),
          100 * round(test$sdev[2] ^ 2 / sum(test$sdev ^ 2), 3),
          100 * round(test$sdev[3] ^ 2 / sum(test$sdev ^ 2), 3)
        )
      axis_x = paste0('PC1: ', percentVar[1], "% variance")
      axis_y = paste0('PC2: ', percentVar[2], "% variance")
      axis_z = paste0('PC3: ', percentVar[3], "% variance")
      test <- test$x
    }
    else{
      umap_test <-
        umap(
          t(cds),
          n_components = dimension,
          n_neighbors = umap.config[['n_neighbors']],
          min_dist = umap.config[['min_dist']],
          metric = umap.config[['metric']], 
          seed = umap.config[['seed']]
        )
      axis_x = paste0('Component 1')
      axis_y = paste0('Component 2')
      axis_z = paste0('Component 3')
      test <- umap_test$layout
    }
    alpha = 1
    if (dimension == 2)
    {
      theme0 <-
        theme_bw() + theme(
          plot.title = element_text(hjust = 0.5, size = 14),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          aspect.ratio = 1,
          axis.text.y = element_text(size = 16),
          axis.text.x = element_text(size = 16),
        )
      scale_color_manual()
      if (color_by == 'cellType') {
        dot_color <- DOT_COLOR
        colors_used <- dot_color[as.character(unique(meta[[color_by]]))]
        meta$condition <- CONDITIONS_NAMES[as.character(meta[[color_by]])]
        alpha = 1
      }
      test <- data.frame(test[, c(1, 2)], condition = as.factor(meta[['cellType']]), org = as.factor(meta[['organism']]))
      colnames(test) <- c('PC1', 'PC2', 'condition', 'org')
      if (labeling == TRUE)
      {
        plot0 <-
          ggplot(test,
                 label = T,
                 aes(PC1, PC2, color = condition, label = row.names(test), shape = NULL)) + geom_text(size = 2) + geom_point(size = 0) + xlab(axis_x) +
          ylab(axis_y) + coord_fixed() + theme0 + ggtitle(main)
      }
      else
      {
        plot0 <-
          ggplot(test, label = T, aes(PC1, PC2, color = condition, shape = NULL)) + geom_point(size =point_size, alpha = alpha) + xlab(axis_x) +
          ylab(axis_y) + coord_fixed() + theme0 + ggtitle(main)
      }
      if (color_by == 'cellType') {
        plot0 <- plot0 + scale_color_manual(values = colors_used)
      } else{
        plot0 <- plot0
      }
      
    }
    else
    {
      test <-
        data.frame(test[, c(1, 2, 3)], condition = as.factor(meta[[color_by]]))
      colnames(test) <- c('PC1', 'PC2', 'PC3', 'condition')
      if (labeling == TRUE)
      {
        plot0 <- plot_ly(
          test,
          x = ~ PC1,
          y = ~ PC2,
          z = ~ PC3,
          color = ~ condition,
          text = rownames(test)
        ) %>% add_text() %>% layout(scene = list(
          xaxis = list(title = axis_x),
          yaxis = list(title = axis_y),
          zaxis = list(title = axis_z)
        ))
      }
      else
      {
        plot0 <- plot_ly(
          test,
          x = ~ PC1,
          y = ~ PC2,
          z = ~ PC3,
          color = ~ condition
        ) % >% layout(scene = list(
          xaxis = list(title = axis_x),
          yaxis = list(title = axis_y),
          zaxis = list(title = axis_z)
        ))
      }
    }
    return(plot0)
  }

mouse_gtf <- makeTxDbFromGFF('./dataset/mouse_data/mouse_spike.gtf', 'gtf')
rat_gtf <- makeTxDbFromGFF('./dataset/rat_data/rat_spike.gtf', 'gtf')

options(ucscChromosomeNames=FALSE)

plot_utr_coverage <- function(gene, utr_res, txdb, bw_file_list, cols, cell_names, file_name_suffix = '', loci = NULL){
  if(is.null(loci)){
    loci = utr_res[gene,'loci']
  }
  chrom_range <- strsplit(loci, split = ':')[[1]]
  chrom <- chrom_range[1]
  rg <- strsplit(chrom_range[2], '-')[[1]]
  strt <- as.numeric(rg[1])
  nd <- as.numeric(rg[2])
  strd <- utr_res[gene,'strand']
  
  if(strd == '+'){
    start <- strt-200
    end <- nd+10
  }else{
    start <- strt-10
    end <- nd+200
  }
  
  gr <- GRanges(chrom, IRanges(start, end), strand=strd)
  #grW <- parse2GRanges(utr_res[gene,'loci'])
  ids <- getGeneIDsFromTxDb(gr, txdb)
  symbols <- ids
  genes <- geneTrack(ids, txdb, 
                     symbols, asList=FALSE)
  #temp_score <- importScore(bw_file_list[1],
                            #bw_file_list[2],
                         #format="BigWig", ranges = gr)
  temp_score <- sapply(bw_file_list, FUN = function(x){importScore(x, format="BigWig", ranges = gr)})
  #setTrackStyleParam(temp_score, "color", c(DOT_COLOR['mouseZygote'], DOT_COLOR['mouseEgg']))
  #strand(trackList[['Mouse']]@dat) <- '+'
  #strand(trackList[['Mouse']]@dat2) <- '-'
  
  #temp <- geneModelFromTxdb(mouse_gtf, org.Mm.eg.db, gr = gr)
  trackList <- trackList(c(genes, temp_score))
  #names(trackList) <- c(gene, 'Mouse')
  
  optSty <- optimizeStyle(trackList, theme="bw")
  #viewerStyle <- trackViewerStyle()
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  names(trackList) <- c(paste(gene, ' (',strd,')',sep=''), cell_names[1], cell_names[2])
  setTrackStyleParam(trackList[[2]], "color", cols[1])
  setTrackStyleParam(trackList[[3]], "color", cols[2])
  setTrackViewerStyleParam(viewerStyle, "xaxis", T)
  #setTrackViewerStyleParam(viewerStyle, "xgp", list(cex = 1.3, col = 'black'))
  setTrackViewerStyleParam(viewerStyle, "margin", c(0.14, 0.15, -0.25, 0.05))
  setTrackXscaleParam(trackList[[1]], "draw", TRUE)
  setTrackStyleParam(trackList[[1]], "ylabpos", "bottomright")
  setTrackStyleParam(trackList[[2]], "ylabpos", "bottomright")
  setTrackStyleParam(trackList[[3]], "ylabpos", "bottomright")
  setTrackStyleParam(trackList[[3]], "ylabgp", list(cex=1.5, col="black"))
  setTrackStyleParam(trackList[[2]], "ylabgp", list(cex=1.5, col="black"))
  setTrackStyleParam(trackList[[1]], "ylabgp", list(cex=0, col="black"))
  setTrackStyleParam(trackList[[1]], "height", 0.1)
  trackList[[1]]@style@yaxis@main <- T
  trackList[[2]]@style@yaxis@main <- T
  trackList[[3]]@style@yaxis@main <- T
  trackList[[2]]@style@yaxis@gp$cex <- 0.8
  trackList[[3]]@style@yaxis@gp$cex <- 0.8
  #trackList
  png(paste(paste(gene,file_name_suffix,sep = '_'), '.png', sep =''), width = 4.2, height = 2.5, units = 'in', res= 300)
  vp <-  viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  addGuideLine(c(strt, nd), vp=vp, col = 'black')
  grid.text(paste(gene, ' (',strd,')',sep=''), 0.05, 0.65, rot = 90, gp = gpar(fontsize = 18, fontface = "bold"))
  dev.off()
  vp
}


plot_range_coverage <- function(range, strand, txdb, bw_file_list, log = T, cols, cell_names,file_name_suffix = '', gene_name = NULL){
  #chrom_range <- strsplit(range, split = ':')[[1]]
  #chrom <- chrom_range[1]
  #rg <- strsplit(chrom_range[2], '-')[[1]]
  #strt <- as.numeric(rg[1])
  #nd <- as.numeric(rg[2])
  strd <- strand
  
  #if(strd == '+'){
    #start <- strt-200
    #end <- nd+10
  #}else{
   # start <- strt-10
  #  end <- nd+200
  #}
  
  #gr <- GRanges(chrom, IRanges(start, end), strand=strd)
  grW <- parse2GRanges(range)
  ids <- getGeneIDsFromTxDb(grW, txdb)
  gene <- ifelse(is.null(gene_name), ids[1], gene_name)
  symbols <- gene
  #gene <- ifelse(is.null(gene_name), ids[1], gene_name)
  genes <- geneTrack(gene, txdb, 
                     gene, asList=FALSE)
  #temp_score <- importScore(bw_file_list[1],
  #bw_file_list[2],
  #format="BigWig", ranges = gr)
  temp_score <- sapply(bw_file_list, FUN = function(x){importScore(x, format="BigWig", ranges = grW)})
  #setTrackStyleParam(temp_score, "color", c(DOT_COLOR['mouseZygote'], DOT_COLOR['mouseEgg']))
  #strand(trackList[['Mouse']]@dat) <- '+'
  #strand(trackList[['Mouse']]@dat2) <- '-'
  
  #temp <- geneModelFromTxdb(mouse_gtf, org.Mm.eg.db, gr = gr)
  trackList <- trackList(c(genes, temp_score))
  if(log){
  trackList[[2]]@dat$score <- log10(trackList[[2]]@dat$score+1)
  trackList[[3]]@dat$score <- log10(trackList[[3]]@dat$score+1)
  }
  #names(trackList) <- c(gene, 'Mouse')
  
  optSty <- optimizeStyle(trackList, theme="bw")
  #viewerStyle <- trackViewerStyle()
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  names(trackList) <- c(paste(gene, ' (',strd,')',sep=''), cell_names[1], cell_names[2])
  setTrackStyleParam(trackList[[2]], "color", cols[1])
  setTrackStyleParam(trackList[[3]], "color", cols[2])
  setTrackViewerStyleParam(viewerStyle, "xaxis", T)
  setTrackViewerStyleParam(viewerStyle, "xgp", list(cex = 0.8, col = 'black'))
  setTrackViewerStyleParam(viewerStyle, "margin", c(0.14, 0.15, -0.25, 0.05))
  setTrackXscaleParam(trackList[[1]], "draw", TRUE)
  setTrackStyleParam(trackList[[1]], "ylabpos", "bottomright")
  setTrackStyleParam(trackList[[2]], "ylabpos", "bottomright")
  setTrackStyleParam(trackList[[3]], "ylabpos", "bottomright")
  setTrackStyleParam(trackList[[3]], "ylabgp", list(cex=1.5, col="black"))
  setTrackStyleParam(trackList[[2]], "ylabgp", list(cex=1.5, col="black"))
  setTrackStyleParam(trackList[[1]], "ylabgp", list(cex=0, col="black"))
  setTrackStyleParam(trackList[[1]], "height", 0.1)
  trackList[[1]]@style@yaxis@main <- T
  trackList[[2]]@style@yaxis@main <- T
  trackList[[3]]@style@yaxis@main <- T
  trackList[[2]]@style@yaxis@gp$cex <- 0.8
  trackList[[3]]@style@yaxis@gp$cex <- 0.8
  #trackList[[1]]@style@xaxis@gp$cex <- 0.8
  #return(trackList)
  png(paste(paste(gene,file_name_suffix,sep = '_'), '.intron.png', sep =''), width = 4.2, height = 2.5, units = 'in', res= 300)
  vp <- viewTracks(trackList, gr=grW, viewerStyle=viewerStyle)
  grid.text(paste(gene, ' (',strd,')',sep=''), 0.05, 0.65, rot = 90, gp = gpar(fontsize = 18, fontface = "bold"))
  dev.off()
  vp
}



poolingCurve <- function(ct,
                         meta,
                         color_by = 'condition',
                         main = 'Pooling Curves',
                         file = 'pooling_curves.pdf') {
  n_cond <- table(meta$condition)
  max_pool <- max(n_cond)
  pdf(file = file, width = 5, height = 5)
  par(
    mfrow = c(1,1),
    mgp=c(3,1,0),
    cex.lab = 2,
    font.main = 2,
    font.lab = 2,
    cex.main = 2,
    cex.axis = 1.7,
    oma = c(4.8,3.8,2,0),
    mar = c(0,0,0.8,0.8) 
  )
  require(zoo)
  plot(NULL,
    xlim = c(1 , max_pool),
    ylim = c(0, 6000),
    xlab = "",
    ylab = "",
    xaxt="n", 
    yaxt="n",
    
  )
  title(
    line = 0,
    main = main, outer = T
  )
  title(
    xlab = expression('Pooled Cells'),
    line = 2,outer = T
  )
  axis(1, tcl=.5, labels = T, padj=-.5)
  axis(2, tcl=.5, labels  = F, padj=1.5)
  
  #axis( 1, 0:1e9, c(  c(  "0", "1", "2", "3", "4" ,"5")) )
  #axis( 2, 0:7000, c( "-2","-1","0","1","2"), las=2 )
  for (cond in unique(meta[[color_by]])) {
    dot_color <- DOT_COLOR[cond]
    samples <- which(meta[[color_by]] == cond)
    sampled_reads <- poolingSampler(ct[, samples])
    #points(sampled_reads[,ncol(sampled_reads)], rollmean(sampled_reads[,i], 5, fill=NA), col= alpha(dot_color, 0.2), lwd=1)
    lines(
      as.integer(names(sampled_reads)),
      sampled_reads,
      type = 'o',
      col = alpha(dot_color, 1),
      lwd = 3,
      lty = "solid",
      pch = 17
    )
    
  }
  dev.off()
}

densityPlots <- function(meta,
           plot = T,
           color_by = 'condition',
           dtype = 'num_read',
           bounds = 'lower',
           main = 'title',
           save = NULL,
           bound_line = T) {
    aspect <- meta[[dtype]]
    
    upper_bound <- 10 ^ (mean(log10(aspect)) +
                           2 * sd(log10(aspect)))
    lower_bound <- 10 ^ (mean(log10(aspect)) -
                           2 * sd(log10(aspect)))
    if (!plot) {
      meta <- meta[aspect > lower_bound &
                     aspect < upper_bound, ]
      return(meta)
    }
    theme0 <-
      theme_bw() + theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12),
        
        axis.text.x = element_text(size = 16))
    data = meta
    data[[color_by]] <- as.factor(data[[color_by]])
    colnames(data)[colnames(data) == color_by] = 'condition0'
    if (color_by == 'condition') {
      dot_color <-
        c("#8C4900",
          "#008331",
          "#005FAD" ,
          "#ffa31a",
          "#B20000",
          "#6304EC",
          "#8218f2")
      names(dot_color) <-
        CONDITIONS_NAMES[c(
          'hypotonic',
          'isotonic',
          '-aa',
          '-aa_population',
          '-glucose',
          'hypotonic_C1',
          'isotonic_hypotonic'
        )]
      colors_used <-
        dot_color[CONDITIONS_NAMES[as.character(unique(data[['condition0']]))]]
      data$condition0 <-
        CONDITIONS_NAMES[as.character(data[['condition0']])]
    }
    plot0 <-
      ggplot(data, aes_string(x = dtype , fill = 'condition0')) + geom_density(alpha =0.5) +
      guides(fill = guide_legend(title = color_by)) + theme0 + ggtitle(main)
    if (bound_line) {
      plot0 <- plot0 + geom_vline(xintercept = lower_bound, linetype = "dashed") + geom_vline(xintercept = upper_bound, linetype = "dashed")
    }
    if (color_by == 'condition') {
      plot0 <- plot0 + scale_fill_manual(values = colors_used)
    }
    if (!is.null(save)) {
      ggsave(
        paste('./', save, '.png', sep = ''),
        dpi = 300,
        width = 5,
        height = 4
      )
    }
    plot0 <- plot0+theme(plot.margin = unit(c(-0.3, 0.2, 0.5, 0.5), "cm"),
           plot.title = element_text(hjust = 0.5, vjust = -3, face = 'bold'),
           axis.title.x = element_blank(),
           axis.title.y = element_blank(), 
           axis.ticks.length=unit(-0.1, "cm"), 
           legend.position = "none", 
           axis.text.y = element_blank(),
           axis.text.x = element_text(size = 16, vjust = -1), 
           aspect.ratio = NULL
       )
    return(plot0)
  }




sample_heatmap <- function(cds,
           color = "Spectral",
           distmethod = "spearman",
           main = "sample distance heatmap",
           fontsize = 5, legend = T,
           annotations)
  {
    if (distmethod %in% c("pearson", "spearman", "kandall", "cosine", "mcd", "ogk"))
    {
      require(MKmisc)
      cdist <- corDist(t(cds), method = distmethod)
      print("calculating correlation distance")
    }
    else
    {
      cdist <- dist(t(cds), method = distmethod)
      print("calculating non-correlation distance")
    }
    matrix <- as.matrix(cdist)
    colors <- colorRampPalette(rev(brewer.pal(9,  "Spectral")))(255)
    anno_color <-
      list(
        cellType = c(
          
          'Mouse Egg' = "#293acc", 
          'Mouse Zygote' = "#B20000", 
          'Rat Egg' = "#008331", 
          'Rat Zygote' = "#000000"
            
          )[unique(as.character(annotations$cellType))]
      )
    pheatmap::pheatmap(
      matrix,
      clustering_distance_rows = cdist,
      clustering_distance_cols = cdist,
      treeheight_row = 0,
      annotation_col = annotations,
      annotation_colors = anno_color,
      treeheight_col = 10,
      color = colors,
      main = main,
      fontsize = fontsize,
      fontsize_row = 6,
      annotation_names_col = F,
      show_colnames  = F,
      legend = legend, 
      clustering_method = 'complete', annotation_legend = F
    )
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
                        annotations, clust_n = 1, mmethod = 'ward.D2', label_genes = NULL) {
  
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
  ha <- NULL
  if(!is.null(label_genes)){
    if(type(label_genes) == 'character'){
      label_genes <- names(ENSEMBL2ALIAS[match(label_genes, ENSEMBL2ALIAS)])
      label_genes <- label_genes[!is.na(label_genes)]
      label_ind = match(label_genes, row.names(exprs_mat))
      label_genes <- label_genes[!is.na(label_ind)]
      label_ind = label_ind[!is.na(label_ind)]
      ha <- rowAnnotation(foo = anno_mark(at = label_ind, 
                                    labels = ENSEMBL2ALIAS[label_genes]))
    }
  else{
    label_genes <- sample(row.names(exprs_mat), 20)
    label_ind = match(label_genes, row.names(exprs_mat))
    label_genes <- label_genes[!is.na(label_ind)]
    label_ind = label_ind[!is.na(label_ind)]
    ha <- rowAnnotation(foo = anno_mark(at = label_ind, 
                                  labels = ENSEMBL2ALIAS[label_genes]))
  }
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
                      row_dend_width = unit(5, "mm"),
                      show_row_dend = T,
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
                      column_names_gp = gpar(fontsize = 20),
                      column_names_rot = 0,
                      column_names_centered = FALSE,
                      
                      top_annotation = HeatmapAnnotation(condition = annotations$condition,
                                                         col = list(condition = c(control='#0f8330', treated ='#b21500')), show_legend = F),
                     
                      bottom_annotation = HeatmapAnnotation( condition = anno_text(colnames(cds), rot = 0, 
                                                                                 gp = gpar(fontsize = 12), location = 0.5, 
                                                                                 just = "center"), which = c('column')),
                      left_annotation = NULL,
                      right_annotation = ha,
                      
                      km = 1,
                      split = NULL,
                      row_km = 1,
                      row_km_repeats = 1,
                      row_split = NULL,
                      column_km = 1,
                      column_km_repeats = 1,
                      column_split = factor(annotations$condition),
                      gap = unit(1, "mm"),
                      row_gap = unit(3, "mm"),
                      column_gap = unit(1, "mm"),
                      
                      heatmap_width = unit(0.7, "npc"),
                      width = NULL,
                      heatmap_height = unit(1, "npc"),
                      height = NULL,
                      
                      show_heatmap_legend = TRUE,
                      heatmap_legend_param = list(title = 'Z score')
                      
  )
  return(heat_obj)
}

addTermsHeatmap <- function(ct, genes, meta, 
                            exprs_mat=NULL, 
                            dist_row = 'pearson', 
                            dist_col = 'euclidean', 
                            hmap_obj=NULL, 
                            org = 'mouse', 
                            k=NULL, 
                            universe = NULL, 
                            plot = F, 
                            tpm = F){
  genes <- intersect(row.names(ct), genes)
  if(!is.null(exprs_mat)){
    genes <- intersect(row.names(exprs_mat), genes)
  }else{
    if(!tpm){
      exprs_mat <- t(t(ct)/estimateSizeFactorsForMatrix(round(ct)))
    }else{
      exprs_mat <- ct
    }
  }
  mat <- ct[genes,]
  exprs_mat <- exprs_mat[genes,]
  anno <- meta[,'cellType',drop =F]
  if(is.null(hmap_obj)){
    hmap_obj <- geneHeatmap(mat, exprs_mat = t(scale(t(log2(exprs_mat+1)))), annotations = anno, tpm = tpm, dist_row  = dist_row,  dist_col = dist_col, clust_n = 2)
    return(hmap_obj)
  }else{
    hmap_obj <- geneHeatmap(mat,exprs_mat  = t(scale(t(log2(exprs_mat+1)))),  annotations = anno, tpm = tpm, dist_col = dist_col, dist_row  = dist_row, clust_n = k)
  }

  hm_order = row_order(draw(hmap_obj))
  clust_l <- rep(0, length(genes))
  names(clust_l) <- 1:length(genes)
  for(i in 1:length(hm_order)){
   clust_l[hm_order[[i]]] <- i
  }
  if(is.null(universe)){
    universe = genes
  }
  
  test_gene_terms <- lapply(1:length(hm_order), FUN = function(x){
    sig_terms <- enrich_CP(genes[hm_order[[x]]], universe = universe, organisms = org, logFC = NULL, GSE = F, GO_BP_only = T)
    if('GO_BP_ora' %in% names(sig_terms))
    {
      sig_terms <- subset(sig_terms$GO_BP_ora@result, qvalue < 0.05 | pvalue < 0.0001)$Description
    }else{
      return(NULL)
    }
    if(length(sig_terms) > 0){
      return(sig_terms[1:min(3, length(sig_terms))])
    }else{
      return(NULL)}
  })
  clust_l <- factor(clust_l)
  names(test_gene_terms) <-  levels(clust_l)
  test_gene_terms[sapply(test_gene_terms, is.null)] <- NULL
  if(length(intersect(names(clust_l), names(test_gene_terms))) == 0)
  {return(hmap_obj)}
  heatmap_obj <- hmap_obj+rowAnnotation(textbox = anno_textbox(clust_l, test_gene_terms,  gp = gpar(fontsize = 10), word_wrap = T, max_width = unit(60, "mm")))
  if(plot != F){
    png(file=plot, res = 300, width = 7, height = 7, units = 'in')
    draw(heatmap_obj)
    dev.off()
  }
  heatmap_obj
}



remove_heatmap_element <- function(hmap, element = 'legend'){
  num <- which(hmap$gtable$layout$name == element)
  if(length(num) == 0) return(hmap)
  hmap$gtable$widths <- embed_vis$aa$hmap$gtable$widths[-num]
  hmap$gtable$heights <- embed_vis$aa$hmap$gtable$heights[-num]
  hmap$gtable$grobs[[num]] <- NULL
  hmap$gtable$layout <- embed_vis$aa$hmap$gtable$layout[-num,]
  return(hmap)
}

diff_violin <- function(df, gene, label) {
  dot_color <- c(
    "#8C4900",
    "#008331", 
    "#005FAD" , 
    "#f48cf9", 
    "#B20000", 
    "#6304EC", 
    "#8218f2",
    "#8cccf9", 
    "#8cb4f9", 
    "#8c9cf9", 
    "#a08cf9", 
    "#cf8cf9",
    "#f48cf9"
  )
  names(dot_color) <- CONDITIONS_NAMES[c(
    'hypotonic',
    'isotonic',
    '-aa',
    '-aa_population',
    '-glucose',
    'hypotonic_C1',
    'isotonic_hypotonic',
    '5p',
    '10p',
    '20p',
    '100p',
    '1n',
    '10n'
  )]
  theme0 <-
    theme_bw() + theme(
      plot.title = element_text(hjust = 0.5, size = 17),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      aspect.ratio = 1,
      axis.text.y = element_text(size = 15),
      axis.text.x = element_text(size = 9)
    )
  gene_exprs <- df[gene, ]
  test <-
    data.frame(expression = log2(gene_exprs + 1),
               condition = as.factor(CONDITIONS_NAMES[as.character(label)]))
  View(test)
  p <- ggplot(test, aes(x = condition, y = expression, fill = condition)) +
    geom_violin(trim = FALSE) + ggtitle(gene)
  p + scale_fill_manual(values = dot_color[unique(test$condition)]) + theme0
}


# for sampling of detection rate
detectionSampler <- function(ct, num, non_ct){
  det_rate <- do.call(rbind, lapply(num, function(N_reads) {
    reads <- sapply(1:ncol(ct), FUN = function(n){
      x <- ct[,n]
      x <- c(x, non_ct[n]-sum(x))
      p = x/sum(x)
      N_reads = N_reads
      sample_ct <- rmultinom(n=1, size =N_reads, prob=p)
      sum(sample_ct[1:(length(sample_ct)-1)] != 0)
    })
    
    reads
  }))
  det_rate <- data.frame(det_rate)
  det_rate['sample_rate'] <- num
  det_rate
}

# for pooling of detection rate
poolingSampler <- function(ct){
  sums <- colSums(ct)
  #ct <- ct[,order(sums)]
  pool <- c()
  det = ct[,1]
  for( i in 1:ncol(ct)){
    temp <- c()
    for(j in 1:100){
      samples <- sample(ncol(ct), i, replace =F)
      det = sum(rowSums(ct[,samples, drop = F]) != 0)
      temp <- c(temp, det)
    }
    pool <- c(pool, mean(temp))
  }
  pool_rate <- pool
  names(pool_rate) <- 1:ncol(ct)
  pool_rate
}

saturationCurve <- function(ct,
                            meta,
                            color_by = 'condition',
                            sep_plots = F,
                            main = 'Saturation Curves',
                            file = 'saturation_curves.pdf',
                            y_label = '',
                            x_label = '',
                            xticks = T) {
  rate <- round(1e2 + 2 ** c(1:22))
  max_rate <- max(rate)
  rate <- c(rate[2:length(rate)], 8e6, 1.2e7)
  pdf(file = file)
  par(
    mfrow = c(1, 3),
    mgp=c(3,1,0),
    cex.lab = 2,
    font.main = 2,
    font.lab = 2,
    cex.main = 2,
    cex.axis = 1,
    oma = c(4.8,3.8,2,0),
    mar = c(0,0,0.8,0.8) 
  )
  if (color_by == 'condition') {
    sep_plots = T
  }
  require(zoo)
  if (!sep_plots) {
    par(
      mfrow = c(1, 3),
      pty = "s",
      mgp=c(3,1,0),
      cex.lab = 2,
      font.main = 2,
      font.lab = 2,
      cex.main = 2,
      cex.axis = 1,
      oma = c(4.8,3.8,2,0),
      mar = c(0,0,0.8,0.8) 
    )
    plot(
      NULL,
      xlim = c(min(rate) , max(max_rate)),
      ylim = c(0, 6000),
      xlab = x_label,
      ylab = "",
      xaxt="n", 
      yaxt="n",
    )
    axis(1, tcl=.5, labels = xticks, padj=-1.5)
    axis(2, tcl=.5, padj=1.5)
    title(ylab = y_label, line = 2)
  }else{
    par(
      pty = "s",
      mfrow = c(1, 3),
      mgp=c(3,1,0),
      cex.lab = 2,
      font.main = 2,
      font.lab = 2,
      cex.main = 2,
      cex.axis = 1,
      oma = c(3,3.8,2,0),
      mar = c(0,0,0.8,0.8) 
    )
  }
  #axis( 1, 0:1e9, c(  c(  "0", "1", "2", "3", "4" ,"5")) )
  #axis( 2, 0:7000, c( "-2","-1","0","1","2"), las=2 )
  j <- 1
  for (cond in unique(meta[[color_by]])) {
    xtckl = xticks; ytckl = T;
    if(j == 3 | j == 2 | j == 4){
      ytckl = F
    }
    print(j)
    if (sep_plots) {
      plot(
        NULL,
        xlim = c(min(rate) , max(max_rate)),
        ylim = c(0, 6000),
        xlab = "",
        xaxt="n", 
        yaxt="n",
        ylab = "", 
      )
      axis(1, tcl=.5, labels = xtckl, padj=-1.5)
      axis(2, tcl=.5, labels  = ytckl, padj=1.5)
      #title(ylab = expression("Detected Genes"),xlab = expression('Reads'),line = 2.4,main = CONDITIONS_NAMES[cond])
    }
    dot_color <- DOT_COLOR[cond]
    samples <- which(meta[[color_by]] == cond)
    non_ct <- meta[samples, 'num_read']
    sampled_reads <- detectionSampler(ct[, samples], rate, non_ct)
    for (i in 1:length(samples)) {
      #points(sampled_reads[,ncol(sampled_reads)], rollmean(sampled_reads[,i], 5, fill=NA), col= alpha(dot_color, 0.2), lwd=1)
      lines(
        sampled_reads[, ncol(sampled_reads)],
        rollmean(sampled_reads[, i], 5, fill = NA),
        col = alpha(dot_color, 0.2),
        lwd = 1.5,
        lty = "solid"
      )
    }
    mean_rate <- rowMeans(sampled_reads[, 1:length(samples)])
    
    lines(
      sampled_reads[, ncol(sampled_reads)],
      rollmean(mean_rate, 5, fill = NA),
      col = alpha(dot_color, 1),
      lwd = 3,
      lty = "solid"
    )
    j <- j+1
  }
  title(ylab =y_label, xlab =x_label, outer = TRUE, line = 2)
  title(main = main, line = 0, outer = T)
  dev.off()
}


