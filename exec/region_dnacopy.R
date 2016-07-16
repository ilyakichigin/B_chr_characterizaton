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
# Argument3 - plot width (optional)
# Argument4 - plot height (optional)

# Package installation:
#source("http://bioconductor.org/biocLite.R")
#biocLite("DNAcopy")
library(DNAcopy)

# Parameter section
args <- commandArgs(trailingOnly = TRUE)
pos_file <- args[1] # bed file with positions, additional column with coverage
size_file <- args[2] # file with chromosome sizes
if (length(args) == 4) { # pdf size
    pdf_width <- as.numeric(args[3])
    pdf_height <- as.numeric(args[4])} else {
    pdf_width <- 20
    pdf_height <- 20}
plot.type <- 's' # 'w' for whole genome in one picture
left_dist <- TRUE # Calculate distance to left position. Correction for right dist not implemented.
draw_plot <- TRUE # FALSE to skip plotting
max_chr_plot <- 42 # maximum number of chromosomes to do plotting - 7*6 plot 20*20 inch
log_dist <- TRUE # calculate regions based on log10(pairvise distances)
smooth <- TRUE # remove outlier values

# Inputs
pos_df <- read.table(pos_file)
size_df <- read.table(size_file)

# Output filename
basename <- unlist(strsplit(pos_file, '/')) # save pdf to current folder
basename <- basename[length(basename)]
basename <- substr(basename,1,nchar(basename)-8) # remove '.pos.bed' from filename

# Calculate distances between positions
# Leaves only positions on chromosomes listed in .sizes file
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
flt_pos_df <- do.call('rbind', flt_pos_list) # convert list of df's to df
flt_pos_df$V1 <- factor(flt_pos_df$V1) # remove unused chromosomes from factor levels
flt_pos_df$V6 <- log10(flt_pos_df$V5)
flt_pos_df$V7 <- mapply(gsub, 'chr', '', flt_pos_df$V1) 
colnames(flt_pos_df) <- c('chrom','start','end','cov','l.dist','log10.l.dist','num.chrom')

# Plot to pdf: chromosome names with stripped chr
# *correct chromosome sort order, sample name, maploc: need to modify DNAcopy.plot function
if (nrow(size_df) > max_chr_plot) {
  draw_plot <- FALSE
  cat('Number of chromosomes in plot exceeds ', max_chr_plot, '. Plotting will be disabled.', sep = "")
}
if (draw_plot) {
  CNA.object <- CNA(flt_pos_df$log10.l.dist,flt_pos_df$num.chrom,flt_pos_df$start,
                    data.type = 'logratio',sampleid = basename)
  if (smooth) {CNA.object <- smooth.CNA(CNA.object)}
  segment.CNA.object <- segment(CNA.object, verbose=1)
  pdf(paste(basename,'.reg.pdf',sep=''), width = pdf_width, height = pdf_height)
  plot(segment.CNA.object, plot.type=plot.type, xmaploc=T,  
       ylim=c(0,max(CNA.object[,3]))) # 'w' type for all chromosomes
  dev.off()
}

# Call regions: chromosome names as initially
CNA.object <- CNA(flt_pos_df$log10.l.dist,flt_pos_df$chrom,flt_pos_df$start,
                  data.type = 'logratio',sampleid = basename)
if (smooth) {CNA.object <- smooth.CNA(CNA.object)}
segment.CNA.object <- segment(CNA.object, verbose=1)
out_regions <- segment.CNA.object$output

# Clean-up
# rename cols
names(out_regions)[3] <- "reg.start"
names(out_regions)[4] <- "reg.end"
#names(out_regions)[5] <- "positions"
#names(out_regions)[6] <- "log10.mean"
# remove unneeded cols: ID, num.mark, seg.mean
out_regions <- out_regions[,-c(1,5,6)] 

# reg.end - change values from beginning to end of position  
out_regions$reg.end <- apply(out_regions[,c('chrom','reg.end')], 1, function(y) { # ugly look-up
  flt_pos_df[flt_pos_df$chrom==y['chrom'] & flt_pos_df$start == as.numeric(y['reg.end']),]$end
})

# Correct for l distance: shift one position to the left, if seg.mean.left>seg.mean.right
# *check for skipped positions between regions
corr_start_end <- c()
for (i in 1:nrow(out_regions)){
  # start of first position in segment
  start_index <- which(flt_pos_df$chrom == as.character(out_regions[i,'chrom']) 
                       & flt_pos_df$start == as.numeric(out_regions[i,'reg.start']))
  # end of last position in segment
  end_index <- which(flt_pos_df$chrom == as.character(out_regions[i,'chrom']) 
                     & flt_pos_df$end == as.numeric(out_regions[i,'reg.end']))
  # correction
  corr_start <- ifelse(i == 1 || out_regions[i-1,'chrom']!=out_regions[i,'chrom'], # start of chr
                       flt_pos_df[start_index,'start'],
                       ifelse(out_regions[i-1,'log10.mean']>out_regions[i,'log10.mean'], # from lower to higher density?
                              flt_pos_df[start_index-1,'start'], # only case to correct
                              flt_pos_df[start_index,'start']))
  corr_end <- ifelse(i == nrow(out_regions) || out_regions[i,'chrom']!=out_regions[i+1,'chrom'], # end of chr
                     flt_pos_df[end_index,'end'],
                     ifelse(out_regions[i,'log10.mean']>out_regions[i+1,'log10.mean'], # from lower to higher density?
                            flt_pos_df[end_index-1,'end'], # only case to correct
                            flt_pos_df[end_index,'end']))
  corr_start_end <- c(corr_start_end, corr_start, corr_end)               
}
corr_se_matrix <- matrix(corr_start_end, ncol=2, byrow=T)
out_regions$reg.start <- corr_se_matrix[,1]
out_regions$reg.end <- corr_se_matrix[,2]

# Calculate additional statistics for each region: number of markers, pd and read coverage mean and sd, pos sizes
add_stats <- apply(out_regions[,c('chrom','reg.start','reg.end')], 1, function(y) {
  # Subset of positions for each region
  flt_pos_sub <- subset(flt_pos_df, chrom==y['chrom']
                        & start>=as.numeric(y['reg.start']) 
                        & end<=as.numeric(y['reg.end']) )
  # distances between positions within region - all but first
  in_dist <- flt_pos_sub$l.dist[2:nrow(flt_pos_sub)]
  # chrom size from .sizes file 
  chr_size <- as.numeric(subset(size_df,V1==y['chrom']))[2] 
  return(c(nrow(flt_pos_sub), # 1 number of positions
           mean(in_dist),sd(in_dist), # 2 mean and 3 sd of pairwise distances between positions inside region
           mean(flt_pos_sub$cov),sd(flt_pos_sub$cov), # 4 mean and 5 sd coverage (actually, number of reads within positions)
           chr_size, # 6 size of chromosomes
           mean(flt_pos_sub$end-flt_pos_sub$start),sum(flt_pos_sub$end-flt_pos_sub$start) # 7 mean and 8 total size of positions within the region
           ))
})
# add to output
out_regions$posiions <- add_stats[1,]
out_regions$pd.mean <- round(add_stats[2,], digits=0)
out_regions$pd.sd <- round(add_stats[3,], digits=0)
out_regions$cov.mean <- round(add_stats[4,], digits=2)
out_regions$cov.sd <- round(add_stats[5,], digits=2)
out_regions$chrom.size <- add_stats[6,]
out_regions$pos.bp.mean <- round(add_stats[7,], digits=0)
out_regions$pos.bp.total <- add_stats[8,]
out_regions$pos.cov <- out_regions$pos.bp.total/(out_regions$reg.end-out_regions$reg.start)

# Write tsv file
write.table(out_regions,file=paste(basename,'.reg.tsv',sep=''),quote=F,sep='\t',
            row.names=F,col.names=T)
