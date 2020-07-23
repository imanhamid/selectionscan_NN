#! /usr/bin/env Rscript

#outputs .png image mask file for each simulation
#output image set to 200x200 at 500kb per pixel
#0: neutral, 1: sweep, 2: unlabeled
#usage: unlabeled_masks.R /full/path/to/tracts.txt

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

#plot image segmentation mask & output as PNG
label_pos <- as.integer(ceiling((bene_locus+1) / chrom_length * 200))
mask <- c(rep(0, times=label_pos-1), 1, rep(0, times=new_length - label_pos), rep(2, times = 199*200))
mask_mat <- matrix(mask, 200, 200, byrow=F)

png(file=paste0(pattern_match[,1], "_P.png"), width=200, height = 200, type = "cairo")
par(mar=c(0, 0, 0, 0))
image(mask_mat, useRaster=TRUE, axes=FALSE, col=grey((0:2)/255))
dev.off()

#print seed, position, and selection coeff to stdout
cat(seed, bene_locus, label_pos, paste0(s, "\n"), sep="\t")
