#!/usr/bin/Rscript


#This script creates 3 pdf graphs describing control libraries,
#as first argument it uses name of bam which was used to create bed files per chromosome,
#as second argument it uses number of chromosomes in reference genome

args <- commandArgs(trailingOnly = TRUE)

y <- NULL
x <- NULL
z <- NULL
for (n in 1:as.numeric(args[2]) )
{
  x = read.table(paste (as.character(args[1]), '.pos.chr', as.character(n),'.bed', sep = ""))
  y[as.character(n)] <- data.frame(x$V4)
}
x = read.table (paste (as.character(args[1]), '.pos.chr', 'X.bed', sep = ""))
y['X'] <- data.frame(x$V4)
for (n in 1:as.numeric(args[2]) )
{
  x = read.table(paste (as.character(args[1]), '.pos.chr', as.character(n),'.bed', sep = ""))
  z[as.character(n)] <- data.frame(x$V2[-1]) - head(data.frame(x$V3),-1)
}
x = read.table (paste (as.character(args[1]), '.pos.chr', 'X.bed', sep = ""))
z['X'] <- data.frame(x$V2[-1]) - head(data.frame(x$V3),-1)
f=function(x) 1000000/x
c=lapply(z,f)
options(scipen=10000)

pdf(paste (as.character(args[1]), "Position_coverage.pdf", sep = "."))
boxplot (y, log='y', varwidth = TRUE, pch=16, cex = 0.5, xlab='Chromosome', ylab='Position Coverage', cex.axis =0.4, cex.lab=1.1, yaxt="n")
axis(2,cex.axis=0.7)
dev.off()

pdf(paste (as.character(args[1]), "Position_pairwise_dist.pdf", sep = "."))
boxplot (z, log='y', varwidth = TRUE, pch=16, cex = 0.5, xlab='Chromosome', ylab='Position Pairwise Distance', cex.axis =0.4, cex.lab=1.1, yaxt="n")
axis(2,cex.axis=1.1)
dev.off()

pdf(paste (as.character(args[1]), "Positions_per_Mb.pdf", sep = "."))
boxplot (c, log='y', varwidth = TRUE, pch=16, cex = 0.5, xlab='Chromosome', ylab='Positions per Mb', cex.axis =0.4, cex.lab=1.1, yaxt="n")
axis(2,cex.axis=1.1)
dev.off()
