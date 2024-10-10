read_expression <- function(dir, mode = 'salmon', tx2gene=NULL, tpmType = 'no', dropInfReps = T){
  if(mode %in% c('star', 'hisat')){
    files <- list.files(dir, pattern = paste(mode, 'htseq.ct', sep = '.'), full.names = T)
    sample_names <- as.character(data.frame(strsplit(list.files(dir, pattern = paste(mode, 'htseq.ct', sep = '.')), split = '_'))[1,])
    cts <- do.call(cbind, lapply(files, FUN = function(x){read.table(x, row.names = 1, header = F)}))
    colnames(cts) <- sample_names
    #cts <- cts[1:(nrow(cts)-5),]
    #cts <- cts[rowSums(cts != 0) > 0,]
    return(cts)
  }else if(mode  == 'salmon'){
    tx2gene <- read.csv(tx2gene,sep='\t' ,header = T, col.names = c("TXNAME", "GENEID"))
    salmon_files <- list.files(dir, recursive = T, pattern='*quant.sf', full.names = T)
    names(salmon_files) <- as.character(as.data.frame(strsplit(salmon_files, '/'))[length(strsplit(salmon_files, '/')[[1]])-1,])
    gene.salmon <- tximport(salmon_files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = tpmType, dropInfReps = dropInfReps, importer=read.delim)
    #gene.salmon$counts <- gene.salmon$counts[rowSums(gene.salmon$counts != 0) > 0,]
    #gene.salmon$abundance <- gene.salmon$abundance[rowSums(gene.salmon$abundance != 0) > 0,]
    transcript.salmon <- tximport(salmon_files, type = "salmon", txOut = TRUE, countsFromAbundance = tpmType, dropInfReps = dropInfReps, importer=read.delim)
    #transcript.salmon$counts <- transcript.salmon$counts[rowSums(transcript.salmon$counts != 0) > 0,]
    #transcript.salmon$abundance <- transcript.salmon$abundance[rowSums(transcript.salmon$abundance != 0) > 0,]
    return(list(gene = gene.salmon, transcript = transcript.salmon))
  }else if(mode == 'rsem'){
    file_names <- list.files(dir, recursive = T, pattern='*genes.results', full.names = F)
    gene_files <- list.files(dir, recursive = T, pattern='*genes.results', full.names = T)
    names(gene_files) <- as.character(as.data.frame(strsplit(file_names, '\\.'))[1,])
    gene.rsem <- tximport(gene_files, type = "rsem", txIn = FALSE, txOut = FALSE, importer=read.delim)
    #gene.rsem$counts <- gene.rsem$counts[rowSums(gene.rsem$counts != 0) > 0,]
    #gene.rsem$abundance <- gene.rsem$abundance[rowSums(gene.rsem$abundance != 0) > 0,]
    file_names <- list.files(dir, recursive = T, pattern='*isoforms.results', full.names = F)
    transcript_files <- list.files(dir, recursive = T, pattern='*isoforms.results', full.names = T)
    names(transcript_files) <- as.character(as.data.frame(strsplit(file_names, '\\.'))[1,])
    transcript.rsem <- tximport(transcript_files, type = "rsem", txIn = TRUE, txOut = TRUE, importer=read.delim)
    #transcript.rsem$counts <- transcript.rsem$counts[rowSums(transcript.rsem$counts != 0) > 0,]
    #transcript.rsem$abundance <- transcript.rsem$abundance[rowSums(transcript.rsem$abundance != 0) > 0,]
    return(list(gene = gene.rsem, transcript = transcript.rsem))
  }
  else{
    cat('only support salmon, rsem, hisat (htseq) and star (htseq)\n')
  }
}





get_qorts_summary <- function(directory, file_name = 'QC.summary.txt'){
  rnames <- row.names(read.csv(paste(list.dirs(directory, recursive = F)[[1]], file_name,sep = '/'), row.names = 1, sep = '\t', header = F))
  summary0 <- do.call(cbind, lapply(list.dirs(directory, recursive = F), FUN = function(x){
    if(length(list.files(x)) > 2){
      new <- read.csv(paste(x,file_name,sep = '/'), row.names = 1, sep='\t', header = F)
      rnames <- intersect(row.names(new), rnames)
      new <- new[rnames,1]
      return(new)
    }else{
      print(x)
      return()
    }
  }
  ))
  row.names(summary0) <- rnames

  colnames(summary0) <- do.call(c, lapply(list.dirs(directory, recursive = F), FUN = function(x){strsplit(x, split = '/')[[1]][length(strsplit(x, split = '/')[[1]])]}))
  
  return(summary0)
}


MeanInGroupCorrelation <- function(meta, mat, groupby = 'cellType', log=F){
  ct <- t(t(mat)/edgeR::calcNormFactors(mat)) #using TMM normalization
  ct <- ct[rowSums(ct > 0) > ncol(ct)*0.1,]
  
  #ct <- mat[rowMeans(mat) >= quantile(rowMeans(mat)[rowMeans(mat) > 0], 0.025) & rowMeans(mat) <= quantile(rowMeans(mat)[rowMeans(mat) > 0], 0.975), ]
  if(log){ct <- log(ct+1)}
  group <- meta[[groupby]]
  all_cor <- cor(ct)
  all_rank <- apply(all_cor, 1, rank)
  meta[['intraCor']] <- Reduce(c, tapply(X = 1:ncol(ct), INDEX = as.character(group), FUN = function(x){rowMeans(cor(ct[,x]))}))[row.names(meta)]
  meta[['interCor']] <- Reduce(c, tapply(X = 1:ncol(ct), INDEX = as.character(group), FUN = function(x){rowMeans(all_cor[x,-x])}))[row.names(meta)]
  meta[['intraRank']] <- Reduce(c, tapply(X = 1:ncol(ct), INDEX = as.character(group), FUN = function(x){rowMeans(all_rank[x,x])}))[row.names(meta)]
  meta[['interRank']] <- Reduce(c, tapply(X = 1:ncol(ct), INDEX = as.character(group), FUN = function(x){rowMeans(all_rank[x,-x])}))[row.names(meta)]
  
  meta[['intraDiffinter']] <- meta[['intraCor']] - meta[['interCor']]
  meta[['intraDiffinterRank']] <- meta[['intraRank']] - meta[['interRank']]
  return(meta)
}

prepareMeta <- function(cds, dirc, filter_conditions=vector(), filter_samples=c(), QC=F)
{
  dot_color <- DOT_COLOR
  stat <- read.table(paste(FILES, dirc, "meta.tsv",sep='/'), header=TRUE, row.names = 1)
  cells <- row.names(stat)
  ct <- get_qorts_summary(paste(FILES,dirc,'/QCData/hisat/',sep='/'), file_name = 'QC.geneCounts.formatted.for.DESeq.txt.gz')
  ct <- ct[1:(nrow(ct)-5),]
  spike.reads <- ct[which(grepl("ERCC", row.names(ct))),]
  bio.reads <- ct[which(!grepl("ERCC", row.names(ct))),]
  spike.rate <- round(colSums(spike.reads)/colSums(ct), 5)
  names(spike.rate) <- colnames(ct)
  stat$spike.rate <- spike.rate[row.names(stat)]
  eff.reads.ratio <- (colSums(ct) - colSums(spike.reads))/(stat[colnames(ct),]$input.read.pair.count*(1-stat[colnames(ct),]$spike.rate))
  names(eff.reads.ratio) <- colnames(ct)
  stat$eff.reads.ratio <- eff.reads.ratio[row.names(stat)]
  stat$num_bio_reads <- colSums(bio.reads)[row.names(stat)]
  stat$sensitivity <- colSums(bio.reads > 0)[row.names(stat)]
  meta <- data.frame(stat[cells,])
  remove_samples <- union(filter_samples, row.names(meta[meta$condition %in% filter_conditions,]))
  meta <- meta[!row.names(meta) %in% remove_samples,]
  meta$color <- dot_color[meta$cellType]
  cds$bio.reads <- bio.reads
  cds$meta <- meta
  if(QC){
  cds$meta <- subset(cds$meta, sensitivity > 500 & complexity > 0.01 & gap < 1.5 & eff.reads.ratio > 0.05 )
  cds$meta <- subset(cds$meta, intraDiffinterRank > 1)
  }
  return(cds)
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


monocleBatchCorrect <- function(cds){
  
  ct <- prepareMonocle(cds, tpm2abs = F, F)
  pData(ct)$batch_date <- as.character(pData(ct)$batch_date)
  cds$ct$bio <- round(customDRBatchCorrection2(ct)[, colnames(cds$ct$bio)])
  cds$ct$bio <- apply(cds$ct$bio,2,function(x) {storage.mode(x) <- 'integer'; x})
  tpm <- prepareMonocle(cds, tpm2abs =F, tpm =T, ct2tpm = F)
  pData(tpm)$batch_date <- as.character(pData(tpm)$batch_date)
  tpm@expressionFamily@vfamily  <- 'Tobit'
  cds$tpm$bio <- customDRBatchCorrection2(tpm, tpm = T)[, colnames(cds$tpm$bio)]
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
  # retrieve data from biomaRt
  #ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = org_type, host="www.ensembl.org")
  #terms <- getBM(attributes = c("entrezgene_accession", 
   #                             "ensembl_gene_id",
  #                              "external_gene_name", 
   #                             "go_id", "name_1006",
  #                              "namespace_1003",
   #                             "definition_1006", 
  #                              "gene_biotype"
  #), filters = "entrezgene_accession", values = rownames(cds$tpm$bio), mart = ensembl)
  gtf = rtracklayer::import(gtf_file)
  gtf = data.frame(gtf[gtf$type == 'gene',])
  gtf <- gtf[,colSums(is.na(gtf)) != nrow(gtf)]
  gtf <- data.frame(row.names = gtf$gene_id, gtf)
  gtf$gene_short_name <- gtf$gene
  gtf <- gtf[genes,]
  cds$features <- gtf
  cds$features_tpm <- cds$features[row.names(cds$tpm$bio),]
  cds$features_ct <- cds$features[row.names(cds$ct$bio),]
  
  # cell cycle ids
  #cc.id <- terms[grepl("cell cycle", terms$definition_1006) | grepl("cell cycle", terms$name_1006), ]
  #cc.id <- cc.id[!grepl("meio", cc.id$name_1006) & !grepl("meio", cc.id$definition_1006), ]
  #cc.id <- cc.id[!grepl("cycle-independent", cc.id$definition_1006), ]
  #cc.id <- cc.id[cc.id$namespace_1003 == "biological_process",]
  #cyclebase <- read.table(paste(FILES,"cerevisiae_periodic.tsv",sep=''), header=TRUE, row.names=2)
  #cyclebase <- cyclebase[cyclebase$rank <= 600, ]
  #cycle.genes <- cc.id[cc.id$go_id == "GO:0007049",]$entrezgene_accession
  #cds$cycle.genes <- cycle.genes
  # different enrichment environments
  # go environment #
  #go <- getBM(attributes = c("entrezgene_accession", "namespace_1003","go_id", "name_1006"), filters = "entrezgene_accession", values = genes, mart = ensembl)
  #go <- go[go$namespace_1003 == "biological_process",]
  #go$term <- go$name_1006
  #s = split(go$ensembl_gene_id, paste(go$go_id, go$term))
  #go.env <- list2env(s)
  #go.env <- list2env(tapply(paste(go$ensembl_gene_id, go$term, go$go_id, sep="_"), 
                            #go$go_id, 
                            #function(x) list(name = paste(strsplit(x[1],"_")[[1]][3], strsplit(x[1],"_")[[1]][2]), 
                                             #genes = unname(sapply(x, function(a) strsplit(a,"_")[[1]][1])))
  #))
  
  # biocyc envrionment #
  #biocyc <- read.table(paste(FILES,"biocyc_gene.txt",sep=''), header=TRUE)
  #biocyc_id <- read.delim(paste(FILES,"biocyc_id.txt",sep=''), header=TRUE, sep="\t", stringsAsFactors = F)
  #ids <- as.character(biocyc_id[,2])
  #names(ids) <- as.character(biocyc_id[,1])
  #cds$"biocyc_gene2path" <- biocyc
  #biocyc.env <- list2env(tapply(paste(cds$biocyc_gene2path$gene, cds$biocyc_gene2path$pathway), 
                                #cds$biocyc_gene2path$pathway, 
                                #function(x) list(name = ids[strsplit(x[1]," ")[[1]][2]], genes = unname(sapply(x, function(a) strsplit(a," ")[[1]][1])))))
  
  # kegg envrionment #
  mapPathway2Name2Gene <- function(organism) 
  {
    KEGG_PATHWAY_GENE <- "http://rest.kegg.jp/link/pathway/"
    KEGG_PATHWAY_LIST_BASE <- "http://rest.kegg.jp/list/pathway/"
    pathway_list_REST_url <- paste(KEGG_PATHWAY_LIST_BASE, organism, sep="")
    pathway_list_REST_url_gene <- paste(KEGG_PATHWAY_GENE, organism, sep="")
    pathway_id_gene <- list()
    for (line in readLines(pathway_list_REST_url)) 
    {
      tmp <- strsplit(line, "\t")[[1]]
      pathway_id <- strsplit(tmp[1], organism)[[1]][2]
      pathway_name <- tmp[2]
      pathway_name <- strsplit(pathway_name, "\\s+-\\s+")[[1]][1]
      pathway_id_gene[[pathway_id]] =list(name=paste(pathway_id,pathway_name, sep=":"), genes=vector())
    }
    #colnames(pathway_id_name)[1] <- "pathway_name"
    for (line in readLines(pathway_list_REST_url_gene))
    {
      tmp <- strsplit(line,"\t")[[1]]
      gene_id <- strsplit(tmp[1], ":")[[1]][2]
      pathway_id <- strsplit(strsplit(tmp[2], ":")[[1]][2], "sce")[[1]][2]
      pathway_id_gene[[pathway_id]][["genes"]][length(pathway_id_gene[[pathway_id]][["genes"]])+1] <- gene_id
    }
    return(pathway_id_gene)
  }
  
  #kegg <- mapPathway2Name2Gene(kegg_org)
  #kegg.env <- list2env(kegg)
  #cds$"go.env" <- go.env
  #cds$"biocyc.env" <- biocyc.env
  #cds$"kegg.env" <- kegg.env
  
  # Transcription Factors
  #TFgenes <- read.csv(paste(FILES,dirc, "TF.txt",sep='/'), header = T, sep ='\t')$Ensembl
  #TFensembl.ids <- unique(cds$features[cds$features$external_gene_name %in% TFgenes | 
                                         #cds$features$ensembl_gene_id %in% TFgenes ,][,c("entrezgene_accession", "external_gene_name", "ensembl_gene_id")])
  #TFensembl.ids[which(TFensembl.ids$external_gene_name == ""),]$external_gene_name <- TFensembl.ids[which(TFensembl.ids$external_gene_name == ""),]$ensembl_gene_id
  #cds$TFgenes <- TFensembl.ids
  return(cds)
}



# filtered yeast data with ribosomal genes in analysis
mouse2rat <- read.csv('./dataset/mouse2rat.txt', sep ='\t', header = T)
rat2mouse <- read.csv('./dataset/rat2mouse.txt', sep ='\t', header = T)
mouse2rat <- subset(mouse2rat, !is.na(X.id..target.Rat.gene.identical.to.query.gene) & mouse2rat$Rat.orthology.confidence..0.low..1.high. == 1)
rat2mouse <- subset(rat2mouse, !is.na(X.id..target.Mouse.gene.identical.to.query.gene) & rat2mouse$Mouse.orthology.confidence..0.low..1.high. == 1)
one2one <- subset(mouse2rat, Rat.homology.type == 'ortholog_one2one')


# yeast data with all cells except population data


