library(dplyr)
library(BOBaFIT)

seg <- read.csv("data/seg_file.csv", row.names = 1)

seg <-
  seg %>% filter(Chromosome != "chrX" &
                   Chromosome != "chrY" & Chromosome != "chrM")

seg <- seg %>% mutate(CN = 2 ^ (logR / 0.55 + 1))

seg$CN <- ifelse(seg$CN > 6, 6, seg$CN)

seg_form <-seg %>% 
  mutate(chr = stringr::str_remove(pattern = "chr", Chromosome),
         width = End - Start) %>%
  select(ID, chr, start = Start, end = End, width, CN)


seg_form2 <- Popeye(seg_form)

chr_list <- computeNormalChromosomes(seg_form2)

DR_results <- DRrefit(seg_form2, chr_list, verbose = F)

write.csv(as.data.frame(DR_results$corrected_segments), "data/BOBseg_file.csv")
