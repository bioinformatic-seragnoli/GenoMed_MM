#!/usr/bin/Rscript

library(optparse)

option_list <- list(
  make_option(c("--input", "-i"), type="character", help="Path to input file. Required."), 
  make_option(c("--normChr"), type="character", default="c(\"1p\",\"2p\",\"2q\",\"4p\",\"8p\",\"8q\",\"10p\",\"10q\",\"12p\",\"12q\",\"16p\",\"17p\",\"17q\",\"18p\",\"18q\")", help = "Vector of normal chromosome. Default: [%default]"), 
  make_option(c("--refGenome"), type="character", default="hg19", help="Human reference genome. Default: [%default]"),
  make_option(c("--libdir"), type="character", help="Path to BOBaFIT scripts"),
  make_option(c("--output", "-o"), default="./", type="character", help="Path to output directory. Default: [%default]")
)
 
parseobj <- OptionParser(option_list = option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen = 0, stringsAsFactors = FALSE)

## defining command line args
input <- opt$input
normChr <- as.character(eval(parse(text = opt$normChr)))
refGen <- opt$refGenome
libdir <- opt$libdir
outDir <- opt$output 

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(GenomicRanges)))
suppressWarnings(suppressMessages(library(magrittr)))

## importing BOBaFIT R scripts
source(paste0(libdir, "/R/Popeye2.R"))
source(paste0(libdir, "/R/DRrefit2.R"))

## creating output directory
dir.create(outDir, recursive = TRUE)

## BOBaFIT: diploid region correction
message("Diploid region correction with BOBaFIT on ULP segments: initiated!\n")

seg_tmp_df <- read.delim(input)
seg_tmp_df <- seg_tmp_df %>% filter(chrom != "X" & chrom != "Y" & chrom != "M")
seg_tmp_df$logR_Copy_Number <- ifelse(seg_tmp_df$logR_Copy_Number > 6, 6, seg_tmp_df$logR_Copy_Number)

seg_df <- seg_tmp_df %>% 
  mutate(width = end - start) %>%
  select(ID, chr=chrom, start, end, width, CN=logR_Copy_Number)

seg_df <- Popeye2(libdir, seg_df, refGenome = refGen, removeXY = TRUE)
sample_name <- unique(seg_df$ID)

CN_weighted_df <- ComputeWeightedCN(seg_df)
ref_test <- FindRefCluster(CN_weighted_df, normChr)
sample_rep <- WriteSampleRep(sample_name, ref_test)
DR_res <- DRrefit(seg_df, ref_test, sample_rep)

write.csv(as.data.frame(DR_res$corrected_segments), paste0(outDir, "/", sample_name, "_BOB_seg.tsv"), row.names = FALSE)
write.table(DR_res$report, paste0(outDir, "/", sample_name, "_BOB_report.tsv"), row.names = FALSE)

message("Diploid region correction with BOBaFIT on ULP segments: done!\n")

## Categorization at chromosome arm level
message("Chromosome arm categorization: initiated!\n")

ulp_class_df <- as.data.frame(DR_res$corrected_segments) %>%
  group_by(ID, chrarm) %>%
  summarise(weighted_CN = weighted.mean(CN_corrected, width))

ulp_class_df$chrarm <- as.factor(ulp_class_df$chrarm)
ulp_class_df$amp <- as.factor(ifelse(ulp_class_df$weighted_CN >= 2.50, 1, 0))
ulp_class_df$del <- as.factor(ifelse(ulp_class_df$weighted_CN <= 1.50, 1, 0))

write.csv(ulp_class_df, paste0(outDir,"/", sample_name, "_chrarmClass.tsv"), row.names = FALSE)

message("Chromosome arm categorization: done!\n")
