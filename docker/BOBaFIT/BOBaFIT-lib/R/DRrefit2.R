#!/usr/bin/Rscript

ComputeWeightedCN <- function(seg_df){
  
  return(as.data.frame(seg_df %>% 
                         group_by(chrarm) %>%
                         summarise(weighted_mean_CN = weighted.mean(CN, width), .groups = "drop")))
}

FindRefCluster <- function(CN_weighted_df, normChr) {
  
  ref_test <- NULL
  
  test <- try({
    
    cluster_res <- NbClust(CN_weighted_df$weighted_mean_CN, distance = "euclidean", method = "ward.D2", index = "all", 
                           min.nc = 2, max.nc = 6)
    invisible(capture.output(cluster_res, type = c("output", "message")))
    
    cluster_df <- data.frame(chr = CN_weighted_df$chrarm, 
                             cluster = cluster_res$Best.partition,
                             stringsAsFactors = FALSE)
    
    cluster_df_sorted <- cluster_df %>%
      filter(chr %in% normChr) %>%
      group_by(cluster) %>%
      summarise(num = n()) %>%
      arrange(desc(num))
    
    ref_cluster <- cluster_df_sorted$cluster[1]
    ref_df <- cluster_df %>% filter(cluster %in% ref_cluster)
    cluster_normChr <- ref_df$chr
    
  }, silent = TRUE)
  
  if (is(test, "try-error")){
    ref_test <- list(out = "FAIL", chrs = normChr, n_clust = NA)
  } else {
    ref_test <- list(out = "SUCCESS", chrs = cluster_normChr, n_clust = max(cluster_res$Best.partition))
  }
  
  return(ref_test)
}

WriteSampleRep <- function(sample_name, ref_test) {
  
  return(data.frame(sample = sample_name,
                    clustering = ref_test$out,
                    ref_chrs = paste0(ref_test$chrs, collapse = ","), 
                    n_cluster = ref_test$n_clust))
  
}

DRrefit <- function(seg_df, ref_test, sample_rep, maxCN = 6) {

  BOB_out <- list()
  
  seg_df$CN[seg_df$CN == Inf] <- maxCN
  ref_seg_df <- seg_df %>% filter(chrarm %in% ref_test$chrs)

  CN_chr_ref <- as.data.frame(ref_seg_df %>%
                                group_by(chrarm) %>%
                                summarise(weighted_mean_CN = weighted.mean(CN, width), .groups = "drop"))
  
  real_diploid_reg <- median(CN_chr_ref$weighted_mean_CN)
  CF <- 2 - real_diploid_reg
  
  sample_rep$correction_factor <- CF
  sample_rep$correction_class <- ifelse(abs(CF) > 0.5, "REFITTED", ifelse(abs(CF) <= 0.1, "NO CHANGES", "RECALIBRATED"))
  
  seg_df$CN_corrected <- seg_df$CN + CF
  seg_df$CN_corrected <- ifelse(seg_df$CN_corrected < 0, 0.001, seg_df$CN_corrected)
  
  BOB_out[["corrected_segments"]] <- seg_df
  BOB_out[["report"]] <- sample_rep
  
  return(BOB_out)
  
}
