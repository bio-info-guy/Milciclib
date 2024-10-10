library('reticulate')
library('universalmotif')
use_python("/Users/ysu13/opt/anaconda3/bin/python")
source_python('./scripts/getUpstreamSeq.py')
yeast_upstream = get_gene_upstream('./dataset/YEAST/S288C_reference_genome_R64-2-1_20150113/S288C_reference_sequence_R64-2-1_20150113.fsa', 
                                   './dataset/YEAST/S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113.gff',
                                   'sc.db')
get_upstream_seqs(row.names(iso_v_aa$bio.reads), './all_upstreams.fa', yeast_upstream)
get_upstream_seqs(row.names(DEG$isotonic_sorbitol_v_AA_starvation$mast), './mast_aa_upstreams.fa', yeast_upstream)
get_upstream_seqs(row.names(DEG$isotonic_sorbitol_v_Glucose_starvation$mast), './mast_glu_upstreams.fa', yeast_upstream)


system('./scripts/prosampler -i ./mast_aa_upstreams.fa -o mast_aa -t 4.9 -w 4.6 -b 3' )
memes::check_meme_install()


degUpstreamMotiffEnrich <- function(genes, seqs, method = c('streme', 'prosampler'), motif_db, prosampler_path = './scripts/prosampler', temp_fasta = './temp.fa', meme_name = 'meme', control = NULL ){
  get_upstream_seqs(genes,temp_fasta, seqs)
  
    if(method == 'streme'){
      if(is.null(control)){
      ctrl = 'shuffle'
      }else{
        ctrl = control
      }
    meme = memes::runStreme('./temp.fa', control = ctrl)
  }else{
    if(is.null(control)){
      flat_seq <- paste(do.call(c, lapply(genes, FUN = function(x){seqs[[x]]})), collapse = '', sep='')
      kmer_ct <- length(count_kmer(flat_seq))
      sigz = qnorm(0.01/kmer_ct, lower.tail = F); sub_sigz = qnorm(0.05/kmer_ct, lower.tail = F)
      ctrl = 3
    }else{
        ctrl = control
        sigz = 8; sub_sigz = 4.5
    }
    print(sigz)
    system(paste(prosampler_path, '-i', temp_fasta, '-o', meme_name, '-t', sigz, '-w', sub_sigz, '-b', ctrl))
    meme = paste(meme_name, 'meme', sep = '.')
  }
  tomtom_res <- memes::runTomTom(meme, database = motif_db, thresh = 0.05, evalue = F,  min_overlap = 8)
  tomtom_res$tomtom <- lapply(tomtom_res$tomtom, FUN = function(x){x[match(unique(x$match_altname), x$match_altname),]})
  subset(tomtom_res, !is.na(best_match_strand))
}

tomtom2Terms <- function(tt_res){
  tfs <- unique(do.call(c, lapply(tt_res$tomtom, FUN = function(res){
     x_match = subset(res, match_qval < 0.05)$match_altname
     sapply(x_match, FUN = function(x){str_to_upper(substr(x, 1,nchar(x)-1))})
     })))
  nsites <- rep(0, length(tfs))
  names(nsites) <- tfs
  for(i in 1:nrow(tt_res)){
    x_match <- tt_res$tomtom[[i]]$match_altname
    cap_names <- sapply(x_match, FUN = function(x){str_to_upper(substr(x, 1,nchar(x)-1))})
    for(n in cap_names){
      nsites[n] <- max(c(nsites[n], tt_res[i,'nsites']))
    }
  }
  define <- do.call(rbind, lapply(tfs, FUN = function(x){
    tryCatch({
    orf <- unique(as.character(org.Sc.sgdCOMMON2ORF[x]))[[1]]
    as.character(org.Sc.sgdDESCRIPTION[orf])
    }, error = function(cond) { return(NA)})
  }))
  row.names(define) <- tfs
  colnames(define) <- 'definition'
  define <- data.frame(define, nsites = nsites)
  define
}

genes2Terms <- function(deg, top_n= 10, up = TRUE, sig.genes = T){
  if(sig.genes){
  deg <- subset(deg, !is.na(Log2FC) & deg$fdr < 0.05)
  }
  if(up){
    deg <- subset(deg[rev(order(deg$Log2FC))[1:top_n],], Log2FC > 0)
  }else{
    deg <- subset(deg[order(deg$Log2FC)[1:top_n],], Log2FC < 0)
  }
  
  orfs <- row.names(deg)
  print(orfs)
  define <- do.call(rbind, lapply(orfs, FUN = function(x){
    tryCatch({
      gene <- unique(as.character(org.Sc.sgdGENENAME[x]))[[1]]
      description <- as.character(org.Sc.sgdDESCRIPTION[x])
      return(c(gene, description))
    }, error = function(cond) {
      gene <- x
      description <- as.character(org.Sc.sgdDESCRIPTION[x])
      return(c(gene, description))
      })
  }))
  define <- data.frame(define, Log2FC = as.numeric(round(deg$Log2FC, digits = 2)))
  row.names(define) <- orfs
  colnames(define) <- c('gene_name', 'definition', 'Log2FC')
  define
}

plot_motifs <- function(res, top_n = 5, file_name = 'plot', title = ''){
  for(i in 1:length(res$motif)){
    
    motif_name = res$motif[[i]]@name
    print(motif_name)
  high <- min(1+length(res$tomtom[[i]]$match_motif), 1+top_n)
  plt <- view_motifs(c(list(res$motif[[i]]), res$tomtom[[i]]$match_motif[1:top_n]))+ggtitle(title)+theme(
  plot.title = element_text(hjust = 0.5, size = 20))
  ggsave(paste(file_name, motif_name, 'png', sep = '.'),plot = plt, dpi = 200, height = high+0.3, width = 6)
  }
}

add_res <- function(res, df = NULL, name = 'new'){
  if(!is.null(df)){
    new_df <- tomtom2Terms(res)
    df[,name] <- rep('-', nrow(df))
    for(n in row.names(new_df)){
      if(n %in% row.names(df)){
        df[n,name] <- new_df[n,'nsites']
        next}else{
        df[n,] <- rep('-', ncol(df))
        df[n,c(1)] <- new_df[n,c(1)]
        df[n,name] <- new_df[n,'nsites']
        
      }
      
    }
  }else{
    df <- tomtom2Terms(res)
    df[[name]] <- df$nsites
    df$nsites <- NULL
  }
  return(df)
}

