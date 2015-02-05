#!/usr/bin/Rscript

# From positions bed file generate pdf file with plots of pairwise distances, coverage of postions on every chromosome. Used for control libraries.
# Argument - bed file with postions 

args <- commandArgs(trailingOnly = TRUE)
# input data
pos_df <- read.table(args[1])
chr_sizes <- read.table(args[2])
work_chr <- chr_sizes$V1
# create lists chromosome by chromosome
cov_list <- list()
pd_list <- list()
for (chr in work_chr) {
  chr_pos_df <- subset(pos_df, V1 == chr) # subset positions in chromosome
  # coverage 
  chr_cov_list <- list(chr_pos_df$V4) # generate chromosome coverage list
  cov_list <- c(cov_list, chr_cov_list) # append to coverage list
  # pairwise distances
  r_start <- chr_pos_df$V2[2:nrow(chr_pos_df)]
  l_end <- chr_pos_df$V3[1:nrow(chr_pos_df)-1]
  chr_pd_list <- list(r_start - l_end)
  pd_list <- c(pd_list, chr_pd_list)
}
basename <- as.character(args[1])
basename <- substr(basename,1,nchar(basename)-8)
pdf(paste (basename, ".chrom.pdf", sep = ""))
options(scipen=10000) # writes 1000000 insetead of 1e+06
# plot coverage
boxplot(cov_list, log='y', varwidth = TRUE, pch=16, cex = 0.5, xlab='Chromosome', ylab='Position Coverage', cex.axis =0.4, cex.lab=1.1, yaxt="n")
axis(2,cex.axis=0.7)
title(paste (as.character(args[1]),"coverage per chromosome"))
# plot pairwise distances
boxplot(pd_list, log='y', varwidth = TRUE, pch=16, cex = 0.5, xlab='Chromosome', ylab='Position Pairwise Distance', cex.axis =0.4, cex.lab=1.1, yaxt="n")
axis(2,cex.axis=1.1)
title(paste (as.character(args[1]),"pairwise distances per chromosome"))
# plot reverse pairwise distances
boxplot(lapply(pd_list, function(x) 1000000/x), log='y', varwidth = TRUE, pch=16, cex = 0.5, xlab='Chromosome', ylab='Positions per Mb', cex.axis =0.4, cex.lab=1.1, yaxt="n")
axis(2,cex.axis=1.1)
title(paste (as.character(args[1]),"reverse pd per chromosome"))
dev.off()
