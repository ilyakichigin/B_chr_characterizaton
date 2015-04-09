#!/usr/bin/Rscript

# From positions bed file generate bed-like file with regions and their statistics and pdf plot for all chromosomes listed in .sizes file. 
# Steps:
# Calculate pairwise distances between positions
# Call regions based on distances to the position on the left
# Plot chromosomes to '.reg.pdf' file
# Correct regions - shift one position to the left if seg.mean.left>seg.mean.right
# Calculate statistics of distances and coverage for each region
# Write resulting table to '.reg.tsv' file

# Argument1 - '.pos.bed' file with postions
# Argument2 - '.sizes' file with sizes of chromosomes for which regions are called

# Package installation:
#source("http://bioconductor.org/biocLite.R")
#biocLite("DNAcopy")
library(DNAcopy)

# Parameter section
args <- commandArgs(trailingOnly = TRUE)
pos_file <- args[1] # file with positions
size_file <- args[2] # file with chromosome sizes
plot.type <- 's' # 'w' for whole genome in one picture
left_dist <- TRUE # Calculate distance to left position. Correction for right dist not implemented.
draw_plot <- TRUE # FALSE to skip plotting
log_dist <- TRUE # calculate regions based on log(pairvise distances)
# pdf size 
pdf_width <- 7
pdf_height <- 8

# Inputs
pos_df <- read.table(pos_file)
size_df <- read.table(size_file)
# Output filename
id <- unlist(strsplit(pos_file, '/')) # save pdf to current folder
id <- id[length(id)]
id <- substr(id,1,nchar(id)-8) # remove '.pos.bed' from filename

# Calculate distances between positions
# Only chromosomes listed in .sizes file are left 
flt_pos_list <- by(size_df, 1:nrow(size_df), function(chr_data) { # iterate over rows in size_df
  chr_pos_df <- subset(pos_df, V1 == as.character(chr_data$V1)) # chr_data$V1 subset positions in chromosome
  if (left_dist) {
    # dist to position on the left
    r_start <- chr_pos_df$V2 # start coordinates of all positions
    l_end <- c(0, chr_pos_df$V3[1:nrow(chr_pos_df)-1]) # end coordinates of prev position, 0 for first 
  } else {
    # dist to position on the right. Not supported in downstream analysis.
    r_start <- c(chr_pos_df$V2[2:nrow(chr_pos_df)],chr_data$V2) # start coordinates of next position, end of chromosome for last
    l_end <- chr_pos_df$V3 # end coordinates of all positions  
  }
  chr_pos_df$V5 <- r_start-l_end # add distance to df
  chr_pos_df
})
flt_pos_df <- do.call('rbind', flt_pos_list) # convert list of df's in df
flt_pos_df$V1 <- factor(flt_pos_df$V1) # remove unused factors

# Region calling
# Plotting version - log scale y
if (draw_plot) {
  CNA.object <- CNA(log(flt_pos_df$V5),flt_pos_df$V1,flt_pos_df$V2,
                    data.type = 'logratio',sampleid = id)
  smoothed.CNA.object <- smooth.CNA(CNA.object)
  segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)
  pdf(paste(id,'.reg.pdf',sep=''), width = pdf_width, height = pdf_height)
  plot(segment.smoothed.CNA.object, plot.type=plot.type, xmaploc=T,  
       ylim=c(0,max(smoothed.CNA.object[,3]))) # 'w' type for all chromosomes
  dev.off()
}

# File output version - normal scale y
if (log_dist) {
  CNA.object <- CNA(
    log(flt_pos_df$V5),
    flt_pos_df$V1,flt_pos_df$V2,
    data.type = 'logratio',sampleid = id)
} else {
  CNA.object <- CNA(
      flt_pos_df$V5,
      flt_pos_df$V1,flt_pos_df$V2,
      data.type = 'logratio',sampleid = id)  
}
smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)
outdata <- segment.smoothed.CNA.object$output

# loc.end - change from beginnings to ends of positions  
outdata$loc.end <- apply(outdata[,c('chrom','loc.end')], 1, function(y) {
  flt_pos_df[flt_pos_df$V1==y['chrom'] & flt_pos_df$V2 == as.numeric(y['loc.end']),]$V3
})

# Remove col1 (ID), initialize new output variable 
outdata1 <- outdata[,-c(1)]

# Correct for l distance: shift one position to the left, if seg.mean.left>seg.mean.right
corr_start_end <- c()
for (i in 1:nrow(outdata1)){
  #start of first position in segment
  start_index <- which(flt_pos_df$V1 == as.character(outdata[i,'chrom']) 
                       & flt_pos_df$V2 == as.numeric(outdata[i,'loc.start']))
  #end of last position in segment
  end_index <- which(flt_pos_df$V1 == as.character(outdata[i,'chrom']) 
                     & flt_pos_df$V3 == as.numeric(outdata[i,'loc.end']))
  #correction
  corr_start <- ifelse(i == 1 || outdata[i-1,'chrom']!=outdata[i,'chrom'], # start of chr
                       flt_pos_df[start_index,'V2'],
                       ifelse(outdata[i-1,'seg.mean']>outdata[i,'seg.mean'], # from lower to higher density?
                              flt_pos_df[start_index-1,'V2'], # only case to correct
                              flt_pos_df[start_index,'V2']))
  corr_end <- ifelse(i == nrow(outdata1) || outdata[i,'chrom']!=outdata[i+1,'chrom'], # end of chr
                     flt_pos_df[end_index,'V3'],
                     ifelse(outdata[i,'seg.mean']>outdata[i+1,'seg.mean'], # from lower to higher density?
                            flt_pos_df[end_index-1,'V3'], # only case to correct
                            flt_pos_df[end_index,'V3']))
  
  corr_start_end <- c(corr_start_end, corr_start, corr_end)               
}
corr_se_matrix <- matrix(corr_start_end,ncol=2, byrow=T)
outdata1$loc.start <- corr_se_matrix[,1]
outdata1$loc.end <- corr_se_matrix[,2]

# Calculate number of markers, read coverage mean and sd
add_stats <- apply(outdata1[,c('chrom','loc.start','loc.end')], 1, function(y) {
  # Subset of positions for each region
  flt_pos_sub <- subset(flt_pos_df, V1==y['chrom']
                        & V2>=as.numeric(y['loc.start']) 
                        & V3<=as.numeric(y['loc.end']) )
  # distances between positions within region - all but first
  in_dist <- flt_pos_sub$V5[2:nrow(flt_pos_sub)]
  # chrom size from .sizes file 
  chr_size <- as.numeric(subset(size_df,V1==y['chrom']))[2] 
  return(c(nrow(flt_pos_sub),mean(in_dist),sd(in_dist),
           mean(flt_pos_sub$V4),sd(flt_pos_sub$V4),chr_size,
           mean(flt_pos_sub$V5),sd(flt_pos_sub$V5),flt_pos_sub$V5[1]))
})
# add to output
outdata1$num.mark <- add_stats[1,]
outdata1$in.dist.mean <- round(add_stats[2,], digits=0)
outdata1$in.dist.sd <- round(add_stats[3,], digits=0)
outdata1$reads.mean <- round(add_stats[4,], digits=2)
outdata1$reads.sd <- round(add_stats[5,], digits=2)
outdata1$chrom.size <- add_stats[6,]
outdata1$mean.test <- add_stats[7,]
outdata1$sd.test <- add_stats[8,]
outdata1$test <- add_stats[9,]
outdata1$dist.by.size <- round(outdata1$in.dist.mean/outdata1$chrom.size, digits=6) # mean l_dist divided by chromosome size

# Write tsv file
write.table(outdata1,file=paste(id,'.reg.tsv',sep=''),quote=F,sep='\t',
            row.names=F,col.names=T)
