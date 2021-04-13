#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(plyr))

infile <- commandArgs(trailingOnly=TRUE)[1]

#set full chrom length and image width
chrom_length <- as.numeric(commandArgs(trailingOnly=TRUE)[2])
new_length <- as.numeric(commandArgs(trailingOnly=TRUE)[3])

#get filename
pattern <- "(^.+).txt"
pattern_match <- str_match(infile, pattern)

#UNCOMMENT TO GET SEED, S AND POS FROM FILENAME
#pattern <- "^.+/.+_s-(-?\\d.*\\d*)_pos-(\\d+)_seed-(\\d+)_(sample400_)?"
#pattern_match <- str_match(infile, pattern)
#seed <- as.numeric(pattern_match[,4])
#pos <- as.numeric(pattern_match[,3])
#s <- as.numeric(pattern_match[,2])

outfile <- paste0(pattern_match[,1], "ancestry.png")

tracts <- read.table(infile, header = TRUE)

#scale tracts by proportion of chromosome length
tracts$tract_length <- tracts$end_bp - tracts$start_bp
tracts$prop_length <- tracts$tract_length / chrom_length
tracts$scaled_length <- tracts$prop_length*new_length

tracts$colors <- ifelse(tracts$ancID==1, "black", "white")

tracts$childID <- factor(tracts$childID, levels = unique(tracts$childID))
tracts$end_scaled <- ddply(tracts, .(childID), transform, Cumulative.Sum = cumsum(scaled_length))[,"Cumulative.Sum"]

tracts$start_scaled <- ddply(tracts, .(childID), transform, start_scaled = c(0, end_scaled[-length(end_scaled)]))[,"start_scaled"]

#plot tracts & output as png with no whitespace around plot
tractsplot <- ggplot(tracts) +
  geom_segment(mapping=aes(y = childID, yend = childID, x=start_scaled, xend=end_scaled), color = tracts$colors) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  theme(axis.ticks = element_blank(), axis.ticks.length = unit(0, "pt"), axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank(),
        panel.border = element_blank(), aspect.ratio = 1/2, plot.margin = unit(c(0,0,0,0), "null"))

png(file=outfile, width=new_length, height = 200, res=96, type="cairo") #cairo is necessary for running on cluster but may not work on local device
tractsplot
dev.off()


#plot image segmentation mask & output as PNG
#label_pos <- as.integer(ceiling((pos+1) / chrom_length * new_length))

#mask <- c(rep(0, times=label_pos-1), 1, rep(0, times=new_length - label_pos))
#mask <- rep(mask, times=new_length)
#mask_mat <- t(matrix(mask, new_length, new_length, byrow=T))

#png(file=paste0(pattern_match[,1], "ancestry_P.png"), width=new_length, height = new_length, type = "cairo")
#par(mar=c(0, 0, 0, 0))
#image(mask_mat,axes=FALSE, col = grey(c(0,1)))
#dev.off()
