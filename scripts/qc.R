library(ggplot2)
library(ggfortify)
library(QoRTs)
library(plotly)
library(MASS)
library(rrcov)
library(robustbase)
#library(SMNN)
library(reticulate)
use_python("/Users/syc042813/anaconda3/bin/python")
librarOUTLIERS <- c("Y206","Y204","Y207","Y702","Y213", "Y211", "Y217", "Y214", "Y220", "Y208", "Y205", "Y219", "Y203", "Y221", "Y216", "Y213", "Y714", "Y713","Y720", "Y212", "Y209", "Y201", "Y218", "Y310", "Y306","Y328", "Y323", "H6P11", "H6P14", "H6P15" ,"H6P16")

NormalFilter <- function(meta, factor, attribs,  pval=0.023, density="normal", method='union')
{
  
  factors <- unique(meta[, factor])
  if (method == 'union')
  {
    unwanted_samples = c()
  }
  for(i in factors)
  {
    if(method == 'intersect'){
      unwanted_samples = row.names(meta)
    }
    for (a in names(attribs))
    {
      this_group <- log10(meta[meta[,factor] == i,][,a]+1)
      names(this_group) <- row.names(meta[meta[,factor] == i,])
      distr <- MASS::fitdistr(this_group, densfun = density)
      probs <- pnorm(this_group, distr$estimate[1], sd = distr$estimate[2])
      if (attribs[[a]] == 'u')
      {
        filtered_samples <- names(this_group[which(probs >= 1-pval)])
      }else if (attribs[[a]] == 'l'){
        filtered_samples <- names(this_group[which(probs <= pval)])
      }
      unwanted_samples <- eval(parse(text=method))(unwanted_samples, filtered_samples)
    }
    meta <- meta[!row.names(meta) %in% unwanted_samples,]
  }
  return(meta)
}

ManualFilter <- function(meta, list_threshold=NULL,list_samples=NULL)
{
  if (!is.null(list_samples)){
    meta <- meta[!rownames(meta) %in% list_samples,]
    return(meta)
  }
  else{
    if (is.null(list_threshold)){
      list_threshold = list(complexity=c(0.1, 'l'), continuity =c( 0.6, 'u'), evenness = c(2, 'u'), non_spikein_uniq_rate = c(0, 'l'), correlation = c(0.4, 'l'))
    }
    for (i in 1:length(list_threshold)){
      feature = names(list_threshold)[i]
      if (list_threshold[[i]][2] == 'u')
      {
        meta <- meta[which(meta[,feature] < as.numeric(list_threshold[[i]][1])), ]
      }
      else
      {
        meta <- meta[which(meta[,feature] >  as.numeric(list_threshold[[i]][1])), ]
      }
    }
  }
  return(meta)
}

MCDFilter <- function(meta, alpha){
  filter <- robustbase::covMcd(meta[,c(3,5,6,7,8)], alpha = 0.9)
  meta <- meta[row.names(meta)[filter$best],]
  return(meta)
}




mouse.qorts <- subset(read.csv('./dataset/mouse_data/meta.tsv', row.names = 1, header = T, sep = '\t'), organism == 'mouse')
rat.qorts <- subset(read.csv('./dataset/rat_data/meta.tsv', row.names = 1, header = T, sep = '\t'), organism == 'rat')
mouse.qorts$input.read.pair.count <- get_qorts_summary('./dataset/mouse_data/hisat/QCData/')['PREALIGNMENT_READ_CT', row.names(mouse.qorts)]
rat.qorts$input.read.pair.count <- get_qorts_summary('./dataset/rat_data/hisat/QCData/')['PREALIGNMENT_READ_CT', row.names(rat.qorts)]

write.table(mouse.qorts,'./dataset/mouse_data/meta.tsv', sep = '\t', quote = F, row.names = T, col.names = T )
write.table(rat.qorts,'./dataset/rat_data/meta.tsv', sep = '\t', quote = F, row.names = T, col.names = T )



mouse.res.pair <- read.qc.results.data("./dataset/mouse_data/hisat/", decoder=mouse.qorts, autodetectMissingSamples = TRUE,  debugMode = F, calc.DESeq2 = T, calc.edgeR = T)
rat.res.pair <- read.qc.results.data("./dataset/rat_data/hisat/", decoder=rat.qorts, autodetectMissingSamples = TRUE, debugMode = F,  calc.DESeq2 = T, calc.edgeR = T)



plot_qorts <- function(res=NULL,
                       dirc=NULL, 
                       qorts=NULL, 
                       plot_types = c('biotype.rates',
                                      'chrom.type.rates', 
                                      'clipping', 
                                      'dropped.rates', 
                                      'gene.assignment.rates', 
                                      'genebody.coverage',
                                      'genebody.coverage.UMQuartile',
                                      'mapping.rates',
                                      'insert.size'),
                       sep = 'condition',
                       plot_path = './qort_plot.pdf'){
  if(!is.null(res)){
    res = res
    qorts = res@decoder
  }else{
    res <- read.qc.results.data(infile.dir = dirc, decoder=qorts, autodetectMissingSamples = TRUE, debugMode = F, calc.DESeq2 = T, calc.edgeR = T)
  }
  #print(qorts$condition)
  colorby <- qorts[res@decoder$unique.ID,sep]
  names(colorby) <- res@decoder$unique.ID
  print(colorby)
  all_plot <- build.plotter.advanced(res, colorBy = as.character(colorby), color.title = 'cellType', plotter.params = list(contrasting.colors = DOT_COLOR[sort(as.character(unique(colorby)))]))
  group_list = list()
  for(c in unique(colorby)){
    print(c)
    group_list[[c]] = build.plotter.advanced(res, plotter.params = list(std.color = DOT_COLOR[c]), highlightBy = colorby, highlight = c , highlightTitle.singular = 'Condition', outgroup.title = 'Others')
  }
  pdf(plot_path)
  plot.new()
  makePlot.legend.box(all_plot)
  for(p in plot_types){
    if(p == 'biotype.rates'){
      params = '(plotter, count.type = "unambigOnly", showTypes = c("protein_coding", "ncRNA", "rRNA", "pseudogene", "UNK"))'
    }else{
      params = '(plotter)'
    }
    par(mfrow=c(1,1))
    plotter = all_plot
    eval(parse(text=paste("makePlot.", p,params, sep="")))
    
    par(mfrow=c(2,2))
    for(c in names(group_list)){
      plotter = group_list[[c]]
      eval(parse(text=paste("makePlot.", p, params, sep="")))
    }
  }
  
  #makePlot.biotype.rates(all_plot)
  #makePlot.chrom.type.rates(all_plot)
  #makePlot.clipping(all_plot)
  #makePlot.dropped.rates(all_plot)
  #makePlot.gene.assignment.rates(all_plot)
  #makePlot.genebody.coverage(all_plot)
  #makePlot.genebody.coverage.UMQuartile(all_plot)
  #makePlot.mapping.rates(all_plot)
  dev.off()
}

qortsplot <- function(res, by="basic", highlight=FALSE, plotterparam=list(), plottype="insert.size", plotparam="", legendpos="topright")
{
  biotype <- res@qc.data[["biotype.counts"]][[1]][which(res@qc.data[["biotype.counts"]][[1]]$COUNT > 5000 ),]$BIOTYPE
  mkplotter <- function(res, type=by, highlight=FALSE)
  {
    if (highlight != FALSE)
    {
      if(type=="Lane")
      {
        return(build.plotter.highlightSample.colorByLane(curr.sample = highlight, res=res, plotter.params = plotterparam))
      }
      else 
      {
        return(build.plotter.highlightSample(curr.sample=highlight, res, plotter.params = plotterparam))
      }
    }
    else if(type != "")
    {
      return( eval(parse(text=paste("build",".plotter.colorBy", by, '(res, plotter.params=plotterparam)', sep=""))))
    }
    else
    {
      return(build.plotter.basic(res, plotter.params=plotterparam))
    }
  }
  
  plotter=mkplotter(res, type=by, highlight=highlight)
  
  if (plotparam == "")
  {
    eval(parse(text=paste("makePlot.", plottype,'(plotter)', sep="")))
    makePlot.legend.over(legendpos,plotter);
  }
  else
  {
    eval(parse(text=paste("makePlot.",plottype,'(plotter,', plotparam,')', sep="")))
    makePlot.legend.over(legendpos,plotter);
  }
}
qortsplot(worms.res.pair, by="Group", plottype = "genebody", legendpos="topleft")
qortsplot(worms.res, by="Lane", plottype = "genebody.coverage.lowExpress", legendpos="bottomright")

filterSummaryQorts <- function(res){
  res <- res[rowVars(res) != 0,]
  res <- res[rowMin(res) >= 0,]
  res <- res[!grepl('length', row.names(res), ignore.case = T) & 
               !grepl('BENCH', row.names(res)) & 
               !grepl('NumberOfChromosomes', row.names(res)) &
               !grepl('SINGLE', row.names(res)), ]
#test[grepl('read', row.names(test), ignore.case = T) & grepl('pair', row.names(test), ignore.case = T),] <- t(t(test[grepl('read', row.names(test), ignore.case = T) & grepl('pair', row.names(test), ignore.case = T),])/test['PREALIGNMENT_READ_CT',])
#test[rowMax(test) > 10000,] <- t(t(test[rowMax(test) > 10000,])/test['PREALIGNMENT_READ_CT',])
return(res)
}


PCA <- function(df, meta, features, labeling=FALSE, point_size=3, dimension=2,  return_matrix=FALSE, main="PCA", log_t=T)
{
  theme0 <- theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.border = element_blank(), panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio=1)
  scale_color_manual()
  dot_color <- c("#1B9E77","#000000", "#8ce6f9" , "#ffa31a", "#ff0000", "#6304EC" )
  names(dot_color) <- c('hypertonic', 'isotonic_sorbitol', 'AA_starvation', 'AA_starvation_population', 'Glucose_starvation', 'isotonic_sorbitol_C1')
  colors_used <- dot_color[as.character(unique(meta$condition))]
  df <- meta[,features]
  if(log_t){
  test <-prcomp(log(df+1))
  }else{
    test <-prcomp(df)
  }
  if (dimension==2)
  {
    percentVar <- c(100*round(test$sdev[1]^2/sum(test$sdev^2),3), 100*round(test$sdev[2]^2/sum(test$sdev^2),3),100*round(test$sdev[3]^2/sum(test$sdev^2),3))
    test <- data.frame(test$x, condition=meta[[condition]])
    if (labeling == TRUE)
    {
      ggplot(test, label=T, aes(PC1, PC2, color=condition, label=row.names(test)))+geom_text(size=2)+scale_color_manual(values=colors_used)+geom_point(size=0)+xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() +ggtitle(main)+coord_fixed() +theme0
    }
    else
    {
      ggplot(test, label=T, aes(PC1, PC2, color=condition))+scale_color_manual(values=colors_used)+geom_point(size=point_size)+xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() +ggtitle(main)+coord_fixed() +theme0
      }
  }
  else
  {
    
    percentVar <- c(100*round(test$sdev[1]^2/sum(test$sdev^2),3), 100*round(test$sdev[2]^2/sum(test$sdev^2),3),100*round(test$sdev[3]^2/sum(test$sdev^2),3))
    test <- data.frame(test$x, condition=meta$condition)
    if (labeling == TRUE)
    {
      plot_ly(test, x = ~PC1, y = ~PC2, z = ~PC3, color = ~condition, text = rownames(test))%>%add_text() %>% layout(title=main,scene=list(xaxis = list(title = paste0("PC1: ",percentVar[1],"% variance")), yaxis = list(title = paste0("PC2: ",percentVar[2],"% variance")), zaxis=list(title = paste0("PC3: ",percentVar[3],"% variance"))))
    }
    else
    {plot_ly(test, x = ~PC1, y = ~PC2, z = ~PC3, color = ~condition)%>% layout(title=main,scene=list(xaxis = list(title = paste0("PC1: ",percentVar[1],"% variance")), yaxis = list(title = paste0("PC2: ",percentVar[2],"% variance")), zaxis=list(title = paste0("PC3: ",percentVar[3],"% variance"))))
    }
  }
}

dot_color <- c("#1B9E77","#000000", "#8ce6f9" , "#ffa31a", "#ff0000", "#6304EC" )
names(dot_color) <- c('hypertonic', 'isotonic_sorbitol', 'AA_starvation', 'AA_starvation_population', 'Glucose_starvation', 'isotonic_sorbitol_C1')

plotbars <- function(meta, main_cex=2, cex_names = 1, width=1000, height=3000, filename='qc_barplots', x=0,y=0, colorby='condition')
{

pos <- order(meta[,colorby])
if(colorby == 'condition'){
  cols <-  DOT_COLOR[as.character(meta[pos,colorby])]
}
else{
  cols=meta[pos,colorby]
}
cols <-  DOT_COLOR[as.character(meta[pos,colorby])]
png(filename, width=width, height=height, res=300)
par(mfrow=c(7,1), xpd=NA, oma = c(4, 1, 1, 1) )
barplot(meta[pos,]$complexity, cex.names = cex_names, las=2, main = "Complexity", cex.main = main_cex, col = cols)
barplot(meta[pos,]$evenness,  cex.names = cex_names, las=2, main="Evenness of Coverage", cex.main = main_cex, col = cols)
barplot(meta[pos,]$gap, cex.names = cex_names, las=2, main="Number of gaps", cex.main = main_cex, col = cols)
barplot(meta[pos,]$sensitivity, cex.names = cex_names, las=2, main = "Sensitivity", cex.main = main_cex, col = cols)
barplot(meta[pos,]$eff.reads.ratio, cex.names = cex_names, las=2, main = "Effective reads ratio", cex.main = main_cex, col = cols, names.arg = row.names(meta[pos,]) )
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend('bottom', cex=2, ncol=1, bty = 'n', inset=c(x,y),
       x.intersp = 1,y.intersp = 1,box.lwd = 0,
       legend = unique(meta[pos,colorby]),
       fill = unique(cols),xpd=TRUE)
dev.off()
par(mfrow=c(1,1))
}
plotbars(worms.meta, main_cex = 2, filename = 'batch_0',20,0.05)
plotbars(worms.meta[worms.meta$batch_date == 1,], main_cex = 2, filename = 'batch_1',20,0.05)
plotbars(worms.meta[worms.meta$batch_date == 2,], main_cex = 2, filename = 'batch_2',20,0.05)
plotbars(worms.meta[worms.meta$batch_date == 3,], main_cex = 2, filename = 'batch_3',20,0.05)
plotbars(worms.meta[worms.meta$batch_date == 4,], main_cex = 2, filename = 'batch_4',20,0.05)
plotbars(worms.meta[worms.meta$batch_date == 5,], main_cex = 2, filename = 'batch_5',20,0.05)

sample_PCA(worms.meta, TRUE)
sapply(colnames(worms.salmon.tpm), 2, FUN = function(x){cor(worms.salmon.tpm[,x], worms.rsem.tpm[row.names(worms.salmon.tpm), x])})
barplot(sapply(row.names(worms_raw$meta), FUN = function(x){cor(worms.salmon.tpm[,x], worms.rsem.tpm[row.names(worms.salmon.tpm), x])})
        , cex.names = 1, las=2, main = "RSEM vs Salmon Cell Correlation", cex.main = 2, col = worms_raw$meta$color)


#################################################################################################
# Comparison with Chen's
chen_correlation <- read.table("chen_correlation.tsv", row.names=1)
chen_sample <- c("H6P11","H6P12","H6P14","H6P15","H6P16","D1","D19","D20","D22","D23","D24","D26X","D27","D28","D29","D3","D30",
                 "D7","D9","D91X","D93X","D94X","D95X","D96X","D97X","C80X","C81X","F20X","F21X","F22X","F23X","F24X","F30X",
                 "F3X","F40X","F55X","F9X","E10X","E11X","G36","G50","G53","G55","G56","G57","A23","A25","A26","A27","A28","A29","A34","A35","A44","A47","A7")
worms_meta_chen <- worms.meta[chen_sample,]
worms_meta_chen$correlation <- chen_correlation[row.names(worms_meta_chen),]

# manual removal based on certain cutoffs and inspection on PCA results of quality attributes

meta_MANUAL <- manual_filter(worms.meta)
# gaussian filtering on each attribute
# at 0.05 this method leaves 162 samples, however includes Y200s which should all be removed

meta_Normal <- NormalFiltering(worms.meta, 'condition', c("non_spikein_uniq_rate", "evenness","correlation", "continuity",  "complexity"), 0.01)



# minimum covariance determinant remove outliers using robust distance 
filter <- covMcd(worms.meta_all[,c(3,5,6,7,8)], alpha = 0.75)
#covPlot(test$X)
meta_MCD_all <- worms.meta_all[row.names(worms.meta_all)[filter$best],]

# Testing NormalFilter Threshold and factor
num_samples=c()
threshold_val=seq(0.5,0,-0.005)
for(i in 1:length(alpha_val)){
  new <- NormalFiltering(worms.meta, factor='group.ID', c("non_spikein_uniq_rate", "evenness","correlation", "continuity",  "complexity"), threshold_val[i])
  num_samples[i] = nrow(new[new$condition != 'AA_starvation_population',])
}
plot(threshold_val, num_samples, pch=20)

#Testing ManualFilter Thresholds
num_samples=c()

threshold_val=seq(1,0,-0.002)
for(i in 1:length(threshold_val)){
  #list_threshold = list(complexity=c(threshold_val[i], 'l'), continuity =c( max(worms.meta$continuity), 'u'), evenness = c(max(worms.meta$evenness), 'u'), non_spikein_uniq_rate = c(min(worms.meta$non_spikein_uniq_rate), 'l'), correlation = c(min(worms.meta$correlation), 'l'))
  new <- NormalFilter(worms.meta, factor = 'group.ID', attribs = list(non_spikein_uniq_rate='u', evenness='u', continuity='u',  complexity='l', sensitivity='l'), method='union', pval = threshold_val[i])
  num_samples[i] = nrow(new)#[new$condition != 'AA_starvation_population',])
}
names(num_samples) <- as.character(round(threshold_val, 3))
plot(threshold_val, num_samples, pch=20)

meta_Normal <- ManualFilter(worms.meta, list_threshold = list(sensitivity=c(500, 'l')))
meta_Normal <- NormalFilter(meta_Normal, 'group.ID', list(non_spikein_uniq_rate='l', evenness='u', continuity='u',  complexity='l'), method='union',0.05)
