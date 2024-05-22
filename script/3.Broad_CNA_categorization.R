library(dplyr)

# ------------------------ SNParray ---------------------------------------
corr_seg <- read.csv("data/ASCATseg_file.csv", row.names = 1)

corr_seg2 <- corr_seg %>% group_by(ID, chrarm) %>% summarise(weighted.CN=weighted.mean(CN,w=width))
corr_seg2$chrarm <- as.factor(corr_seg2$chrarm )

wide_arms <-data.table::dcast(corr_seg2, ID ~ chrarm)


arms <- colnames(wide_arms[-1])

CNA_db <- data.frame(id=wide_arms$ID)
a=1

for (a in seq_along(arms)) {
  chr <- arms[a]
  
  del <- ifelse(wide_arms[,arms[a]] <= 1.50, 1,0)
  amp <- ifelse(wide_arms[,arms[a]] >= 2.50, 1,0)
  
  alt <- data.frame(amp,del)
  colnames(alt) <- c(paste0("amp_",chr),paste0("del_",chr))
  CNA_db <- cbind(CNA_db,alt)
  
}

write.csv(CNA_db,"data/Broad_CNA_SNParray.csv")

#-------------------------------- ULP-WGS -----------------------------------
ulp_seg <- read.csv("data/ULPseg_file.csv", row.names = 1)

ulp_seg2 <- ulp_seg %>% group_by(sample_sheet_ID, chrarm) %>% summarise(weighted.CN=weighted.mean(CN,w=width))
ulp_seg2$chrarm <- as.factor(ulp_seg2$chrarm )

ulp_wide_arms <-data.table::dcast(ulp_seg2, sample_sheet_ID ~ chrarm)

arms <- colnames(ulp_wide_arms[-1])

ulp_CNA_db <- data.frame(id=ulp_wide_arms$sample_sheet_ID)
a=1

for (a in seq_along(arms)) {
  chr <- arms[a]
  
  del <- ifelse(ulp_wide_arms[,arms[a]] <= 1.50, 1,0)
  amp <- ifelse(ulp_wide_arms[,arms[a]] >= 2.50, 1,0)
  
  alt <- data.frame(amp,del)
  colnames(alt) <- c(paste0("amp_",chr),paste0("del_",chr))
  ulp_CNA_db <- cbind(ulp_CNA_db,alt)
  
}


# ------------- mpc ----------------------

mpc <- read.csv("data/name_codes.csv", sep = ";")

ulp_CNA_db2 <- left_join(ulp_CNA_db,mpc, by=c("id" = "exp_name"))
ulp_CNA_db2$Method <- 1
ulp_CNA_db3 <- ulp_CNA_db2[,-1]

CNA_db2 <- left_join(CNA_db,mpc, by=c("id" = "exp_name"))
CNA_db2$Method <- 2
CNA_db3 <- CNA_db2[,-1]

alt <- colnames(CNA_db2[,-c(1,80:81)]) %>% stringr::str_sort(numeric = T)

ord <- c("MPC", "Method", alt)

ulp_CNA_db3 <- ulp_CNA_db3 %>%  select(ord)

genomic_db <- rbind(ulp_CNA_db3, CNA_db3)

write.csv(genomic_db, "script/final_genomic_dataset.csv")
