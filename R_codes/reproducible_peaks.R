library(GenomicAlignments)
library(data.table)
library(dplyr)
library(GenomicRanges)

setwd("/Users/ziqi_fu/Documents/R_scripts/Ecoli/gapr_76_chip/")

peak_file <- c('SRR13308441_peaks.narrowPeak','SRR13308442_peaks.narrowPeak')

peak_data <- lapply(peak_file, function(x){
  tmp <- as.data.frame(fread(x))
  GRanges(seqnames=tmp[,1],ranges=IRanges(start=tmp[,2],end=tmp[,3]), FDR = tmp[,9])
})


##minoverlap specifies the minimum overlap between two peaks
overlap <- findOverlaps(peak_data[[1]], peak_data[[2]], minoverlap=10L)

peak_data_filter <- list()
peak_data_filter[[1]] <- peak_data[[1]][queryHits(overlap)]
peak_data_filter[[2]] <- peak_data[[2]][subjectHits(overlap)]

peak_com <- Reduce(c,peak_data_filter) %>% reduce()

write.csv(as.data.frame(peak_com), file="reproducible_peak_gapr76.csv", row.names=FALSE)

my_peak_data = peak_data #for my manipulation













