#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(plyr))

infile <- commandArgs(trailingOnly=TRUE)[1]

#set full chrom length and image width
chrom_length <- 25e6
new_length <- 200

#get seed number and selection coeff
pattern <- "^.+/two-strong_seed-(\\d+)_"
pattern_match <- str_match(infile, pattern)
seed <- pattern_match[,2]

outfile <- paste0(pattern_match[,1], "ancestry.png")

#tracts <- read.table("/Users/iman/Desktop/strong_mut_alltracts.txt", header = TRUE)
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
        panel.border = element_blank(), aspect.ratio = 1, plot.margin = unit(c(0,0,0,0), "null"))

png(file=outfile, width=new_length, height = new_length, res=96, type="cairo") #cairo is necessary for running on cluster but may not work on local device
#png(file="/Users/iman/Desktop/strong_mut_out.png", width=new_length, height = new_length, res=96)
tractsplot
dev.off()

#s <- read.table("/Users/iman/Desktop/strong_mut_out.txt", header=TRUE)

s_file <- paste0(pattern_match[,1], "variants.txt")

s <- read.table(s_file, header=TRUE)

#plot image segmentation mask & output as PNG
label_pos <- as.integer(ceiling((s$p1_pos+1) / chrom_length * new_length))
s$p1_label_pos <- label_pos

sum_s <- s %>% dplyr::group_by(p1_label_pos) %>% dplyr::summarise(sum=sum(p1_s))

total_pos <- data.frame(seq(1:new_length))
names(total_pos) <- "pos"

total_pos[sum_s$p1_label_pos, "sum_s"] <- sum_s$sum
total_pos$s <- ifelse(is.na(total_pos$sum_s), 0, total_pos$sum_s)

total_pos$cols <- ifelse(total_pos$s > 0, 1, 0)


mask <- rep(total_pos$cols, times=new_length)
mask_mat <- t(matrix(mask, new_length, new_length, byrow=T))

png(file=paste0(pattern_match[,1], "ancestry_P.png"), width=new_length, height = new_length, type = "cairo")
par(mar=c(0, 0, 0, 0))
image(mask_mat,axes=FALSE, col = grey(c(0,1)))
dev.off()

#png(file=paste0("/Users/iman/Desktop/test_strong_bw", "_P.png"), width=200, height = 200)
#par(mar=c(0, 0, 0, 0))
#image(mask_mat, axes=FALSE, col = grey(c(0,1)))
#dev.off()
