
post_dapars_pas_filter <- function(dapars_res, fa_seq_file, up_range = 400, down_range = 400, offset = 0){
  require(Rsamtools)
  require(BSgenome)
  PAS_SIG = c('AATAAA','ATTAAA','TATAAA','AGTAAA',
              'AATACA','CATAAA','AATATA','GATAAA',
              'AATGAA','AATAAT','AAGAAA','ACTAAA',
              'AATAGA','ATTACA','AACAAA','ATTATA',
              'AACAAG','AATAAG')
  A_Rich_SIG = c('AAAAAA')
  T_Rich_SIG = c('TTTTTT')
  
  
  UTR = do.call(rbind, lapply(1:nrow(dapars_res), FUN = function(i){
    strand = dapars_res$strand[i]
    #pas = ifelse(strand == '+', dapars_res$predicted_p_APA[i]+offset, dapars_res$predicted_p_APA[i]-offset)
    utr = as.numeric(strsplit(strsplit(dapars_res$loci[i], ':')[[1]][2], '-')[[1]])
    #random sample a PAS
    #pas = sample(as.numeric(utr[1]+25):as.numeric(utr[2]-25), 1)
    st = as.numeric(utr[1])
    en = as.numeric(utr[2])
    
    
    c(st, en)
  }))
  
  
  start_end = do.call(rbind, lapply(1:nrow(dapars_res), FUN = function(i){
    strand = dapars_res$strand[i]
    pas = ifelse(strand == '+', dapars_res$predicted_p_APA[i]+offset, dapars_res$predicted_p_APA[i]-offset)
    utr = as.numeric(strsplit(strsplit(dapars_res$loci[i], ':')[[1]][2], '-')[[1]])
    #random sample a PAS
    #pas = sample(as.numeric(utr[1]+25):as.numeric(utr[2]-25), 1)
    utr_len = abs(as.numeric(utr[1])-as.numeric(utr[2]))

    if(strand == '-'){
      new_st = max(as.numeric(utr[1])+60, pas-down_range)
      new_en = min(as.numeric(utr[2]), pas+up_range)
    }
    else{
      new_st = max(as.numeric(utr[1]), pas-up_range)
      new_en = min(as.numeric(utr[2])-60, pas+down_range)
    }
    pas_loc = ifelse(strand == '+', pas-new_st, new_en-pas)
    if(new_en <= new_st){
      return(c(dapars_res$predicted_p_APA[i]-2, dapars_res$predicted_p_APA[i]+2, 3))
    }
    c(new_st, new_en, pas_loc)
  }))
  utr_grange <- makeGRangesFromDataFrame(data.frame(
    chrom=sapply(strsplit(dapars_res$loci, ':'), FUN =function(x){x[1]}), 
    start = UTR[,1], 
    end = UTR[,2], 
    strand = dapars_res$strand))

  test_grange <- makeGRangesFromDataFrame(data.frame(
                      chrom=sapply(strsplit(dapars_res$loci, ':'), FUN =function(x){x[1]}), 
                      start = start_end[,1], 
                      end = start_end[,2], 
                      strand = dapars_res$strand))

  test_fa = FaFile(fa_seq_file)
  test_seqs = getSeq(test_fa, test_grange)
  
  UTR_seq = as.character(getSeq(test_fa, utr_grange))
  #A_rich_range <- test_grange
  #A_rich_range@ranges@start <- A_rich_range@ranges@start+100
  seqs = as.character(test_seqs)
  motif_match_dist = as.data.frame(do.call(rbind, lapply(1:length(seqs), FUN = function(i){
    s = seqs[i]
    locs = do.call(c, lapply(str_locate_all(s, PAS_SIG), rowMeans))
    if(length(locs) == 0){return(c(num_motif = 0, min_dist=NA))}
    dist_to_sig = start_end[i,3] - locs
    min_d = min(abs(dist_to_sig))
    c(num_motif = length(locs), min_dist = dist_to_sig[abs(dist_to_sig) == min_d][1])
  }
  
  
  
  )))
  
  
  
  priming_match_dist = as.data.frame(do.call(rbind, lapply(1:length(UTR_seq), FUN = function(i){
    s = UTR_seq[i]
    estimatePAS=start_end[i,3]
    
    A_locs = do.call(c, lapply(str_locate_all(s, A_Rich_SIG), rowMeans))
    
    
    T_locs = do.call(c, lapply(str_locate_all(s, T_Rich_SIG), rowMeans))
    if(length(A_locs) == 0)
      { dist_to_Asig=NA
        
        }else{
          dist_to_Asig = min(abs(start_end[i,3] - A_locs))
        }
    if(length(T_locs) == 0){
      dist_to_Tsig=NA
    }else{
      dist_to_Tsig = min(abs(start_end[i,3] - T_locs))
        }
    num_A = length(A_locs)
    num_T= length(T_locs)
    
    c(num_A = num_A, num_T = num_T, dist_to_Asig = dist_to_Asig, dist_to_Tsig = dist_to_Tsig)
  }
  )))
  
  
  row.names(motif_match_dist) <- row.names(dapars_res)
  row.names(priming_match_dist) <- row.names(dapars_res)
  final_df = cbind(dapars_res, motif_match_dist, priming_match_dist )
  final_df$diff <- final_df$num_motif > 0 & final_df$fdr < 0.05 & abs(final_df$mean.diff) > 0.2
  list(ranges = test_grange, PAS_motif = final_df)
}

mouse_PAS_dist <- post_dapars_pas_filter(mouse_dapars_sf_st25$deg, './dataset/igv/mouse/mouse_spike.fa', up_range = 80, down_range = 120, offset = 0)
rat_PAS_dist <- post_dapars_pas_filter(rat_dapars_sf_st25$deg, './dataset/igv/rat/rat_spike.fa', up_range = 80, down_range = 120, offset = 0)

mouse_PAS_dist_1000 <- post_dapars_pas_filter(mouse_dapars_sf_st25$deg, './dataset/igv/mouse/mouse_spike.fa', up_range = 500, down_range = 500, offset = 0)
rat_PAS_dist_1000 <- post_dapars_pas_filter(rat_dapars_sf_st25$deg, './dataset/igv/rat/rat_spike.fa', up_range = 500, down_range = 500, offset = 0)


mouse_dapars_sf_st25_deg$PAS_dist= c('+'=1,'-'=-1)[mouse_dapars_sf_st25_deg$strand]*mouse_dapars_sf_st25_deg$PAS_dist

test_motif_match_u30_d70 = sapply(as.character(test_seqs), FUN = function(s){sum(sapply(PAS_SIG, FUN = function(x){grepl(x, s)}))})
test_motif_match_u10_d90 = sapply(as.character(test_seqs), FUN = function(s){sum(sapply(PAS_SIG, FUN = function(x){grepl(x, s)}))})
test_motif_match_800= sapply(as.character(test_seqs), FUN = function(s){sum(sapply(PAS_SIG, FUN = function(x){grepl(x, s)}))})
test_motif_match = sapply(as.character(test_seqs), FUN = function(s){sum(sapply(PAS_SIG, FUN = function(x){grepl(x, s)}))})





mouse_dapars_sf_st25_deg <- mouse_dapars_sf_st25$deg
mouse_dapars_sf_st25_deg$PAS_signal <- test_motif_match
mouse_dapars_sf_st25_deg$PAS_signal_u80_d20 <- test_motif_match_u80_d20
mouse_dapars_sf_st25_deg$PAS_signal_u30_d70 <- test_motif_match_u30_d70
mouse_dapars_sf_st25_deg$PAS_signal_u_20_d120 <- test_motif_match_u_20_d120
mouse_dapars_sf_st25_deg$PAS_signal_u10_d90 <- test_motif_match_u10_d90
mouse_dapars_sf_st25_deg$PAS_signal_800 <- test_motif_match_800
utr_mouse_sf_st25_gs_stricter <- enrich_CP(unique(subset(mouse_PAS_dist$PAS_motif, num_motif > 0 & diff & abs(mean.diff) >= 0.2 )$gene_short_names), universe = unique(subset(mouse_PAS_dist$PAS_motif, num_motif > 0)$gene_short_names), logFC = NULL, organisms = 'mouse')

utr_rat_sf_st25_gs_strict_long <- enrich_CP(unique(subset(rat_PAS_dist$PAS_motif, num_motif > 0 & diff & mean.diff >= 0.2 )$gene_short_names), universe = unique(subset(rat_PAS_dist$PAS_motif, num_motif > 0)$gene_short_names), logFC = NULL, organisms = 'rat')
utr_rat_sf_st25_gs_strict_short <- enrich_CP(unique(subset(rat_PAS_dist$PAS_motif, num_motif > 0 & diff & mean.diff <= -0.2 )$gene_short_names), universe = unique(subset(rat_PAS_dist$PAS_motif, num_motif > 0)$gene_short_names), logFC = NULL, organisms = 'rat')

utr_mouse_sf_st25_gs_strict_long <- enrich_CP(unique(subset(mouse_PAS_dist$PAS_motif, num_motif > 0 & diff & mean.diff >= 0.2 )$gene_short_names), universe = unique(subset(mouse_PAS_dist$PAS_motif, num_motif > 0)$gene_short_names), logFC = NULL, organisms = 'mouse')
utr_mouse_sf_st25_gs_strict_short <- enrich_CP(unique(subset(mouse_PAS_dist$PAS_motif, num_motif > 0 & diff & mean.diff <= -0.2 )$gene_short_names), universe = unique(subset(mouse_PAS_dist$PAS_motif, num_motif > 0)$gene_short_names), logFC = NULL, organisms = 'mouse')

utr_rat_sf_st25_gs_strict <- enrich_CP(unique(subset(rat_PAS_dist$PAS_motif, num_motif > 0 & diff & abs(mean.diff) >= 0.2 )$gene_short_names), universe = unique(subset(rat_PAS_dist$PAS_motif, num_motif > 0)$gene_short_names), logFC = NULL, organisms = 'rat')
utr_mouse_sf_st25_gs_strict <- enrich_CP(unique(subset(mouse_PAS_dist$PAS_motif, num_motif > 0 & diff & abs(mean.diff) >= 0.2 )$gene_short_names), universe = unique(subset(mouse_PAS_dist$PAS_motif, num_motif > 0)$gene_short_names), logFC = NULL, organisms = 'mouse')

