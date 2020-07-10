#! /usr/bin/env Rscript

#outputs .png image file for each simulation
#output image set to 200x200 at 500kb per pixel
#usage: tractsplots_CNN.R /full/path/to/tracts.csv

suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(plyr))

infile <- commandArgs(trailingOnly=TRUE)[1]

#set full chrom length and image width
chrom_length <- 1e8
new_length <- 200

#get seed number and selection coeff
pattern <- "^.+/.+_pos-(\\d+)_.+_s-(\\d.*\\d*)_seed-(\\d+)_alltracts"
pattern_match <- str_match(infile, pattern)
seed <- pattern_match[,4]
s <- pattern_match[,3]
bene_locus <- as.integer(pattern_match[,2])

outfile <- paste0(pattern_match[,1], ".png")

tracts <- read.table(infile, header = TRUE)

#scale tracts by proportion of chromosome length
tracts$tract_length <- tracts$end_bp - tracts$start_bp
tracts$prop_length <- tracts$tract_length / chrom_length
tracts$scaled_length <- tracts$prop_length*new_length

tracts$colors <- ifelse(tracts$ancID==1, "black", "white")

tracts$childID <- factor(tracts$childID, levels = unique(tracts$childID))
tracts$end_scaled <- ddply(tracts, .(childID), transform, Cumulative.Sum = cumsum(scaled_length))[,"Cumulative.Sum"]

tracts$start_scaled <- ddply(tracts, .(childID), transform, start_scaled = c(0, end_scaled[-length(end_scaled)]))[,"start_scaled"]

#order chromosome by tract length spanning variant
newchildorder <- tracts %>%
  group_by(childID) %>%
  summarise(bene_tract = tract_length[start_bp<=bene_locus & end_bp > bene_locus],
            bene_anc = ancID[start_bp<=bene_locus & end_bp > bene_locus],
            childname = childID[start_bp<=bene_locus & end_bp > bene_locus])
childlevels2 <- newchildorder[order(newchildorder$bene_anc, -newchildorder$bene_tract),"childname"]

tracts$childID <- factor(tracts$childID, levels = childlevels2)

#plot tracts & output as png with no whitespace around plot
tractsplot <- ggplot(tracts) +
  geom_segment(mapping=aes(y = childID, yend = childID, x=start_scaled, xend=end_scaled), color = tracts$colors) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  theme(axis.ticks = element_blank(), axis.ticks.length = unit(0, "pt"), axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank(),
        panel.border = element_blank(), aspect.ratio = 1, plot.margin = unit(c(0,0,0,0), "null"))

png(file=outfile, width=200, height = 200, res=96, type = "cairo") #cairo is necessary for running on cluster but may not work on local device
tractsplot
dev.off()

#plot image segmentation mask & output as PNG
label_pos <- as.integer(ceiling((bene_locus+1) / chrom_length * 200))
mask <- c(rep(0, times=label_pos-1), 1, rep(0, times=new_length - label_pos))
mask <- rep(mask, times=200)
mask_mat <- t(matrix(mask, 200, 200, byrow=T))

png(file=paste0(pattern_match[,1], "_P.png"), width=200, height = 200, type = "cairo")
par(mar=c(0, 0, 0, 0))
image(mask_mat, useRaster=TRUE, axes=FALSE, col = grey(c(0,1)))
dev.off()

#print seed, position,  and selection coeff to stdout
cat(seed, bene_locus, label_pos, paste0(s, "\n"), sep="\t")
