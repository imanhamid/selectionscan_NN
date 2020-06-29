#! /usr/bin/env Rscript

#outputs .png image file for each simulation
#output image set to 500x1200
#usage: tractsplots_CNN.R /full/path/to/tracts.csv

suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(plyr))

infile <- commandArgs(trailingOnly=TRUE)[1]

#set full chrom length and image width
chrom_length <- 1e8
bene_locus <- chrom_length/2 #variant is middle of chrom
new_length <- 1200

#get seed number and selection coeff
pattern <- "^.+/.+_s-(\\d.*\\d*)_seed-(\\d+)_alltracts"
pattern_match <- str_match(infile, pattern)
seed <- pattern_match[,3]
s <- pattern_match[,2]

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
  summarise(bene_tract = tract_length[start_bp<=bene_locus & end_bp >= bene_locus],
            bene_anc = ancID[start_bp<=bene_locus & end_bp >= bene_locus],
            childname = childID[start_bp<=bene_locus & end_bp >= bene_locus])
childlevels2 <- newchildorder[order(newchildorder$bene_anc, -newchildorder$bene_tract),"childname"]

tracts$childID <- factor(tracts$childID, levels = childlevels2)

#plot tracts & output as png
tractsplot <- ggplot(tracts) +
  geom_linerange(mapping=aes(x = childID, ymin=start_scaled, ymax=end_scaled), color = tracts$colors) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1200.1)) +
  theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())

png(file=outfile, width=1200, height = 500, res=96, type = "cairo") #cairo is necessary for running on cluster
tractsplot
dev.off()

#print seed and selection coeff to stdout
cat(seed, paste0(s, "\n"), sep="\t")
