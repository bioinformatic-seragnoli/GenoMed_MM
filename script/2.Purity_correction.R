library(data.table)
library(dplyr)

ascat_pur <- read.csv("data/ascat_pur.csv", row.names = 1)
bob_seg <- read.csv("data/BOBseg_file.csv", row.names = 1)

samples <- unique(bob_seg$ID)
j=samples[1]

ascat_segs <- data.frame()

for (j in samples) {
  
  sample_seg <- bob_seg %>% filter(ID %in% j)
  
  purity <- ascat_pur$purity[ascat_pur$ID == paste0("sample_",j)]
  
  if (rlang::is_empty(purity) == FALSE) {
    sample_seg <- sample_seg %>% mutate(CN_purity_correction = (((CN_corrected-2) / purity) + 2),
                                        ascat_purity = purity)
    
  } else {
    sample_seg <- sample_seg %>% mutate(CN_purity_correction = CN_corrected) %>% mutate(ascat_purity = NA)
    
  }
  
  ascat_segs <- rbind(ascat_segs,sample_seg) 
  
}

write.csv(ascat_segs, "data/ASCATseg_file.csv")  
