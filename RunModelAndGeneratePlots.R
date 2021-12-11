library(ggplot2)
library(dplyr)
library(data.table)

# File and directory paths
pathToScript <- "/path/to/CHiC-Functions.R"  # contains modeling and helper functions
dataDir <- "/path/to/directory/"  # should contain chicago_input.rds, rmap.csv, baitmap.csv, and mappability.bigWig
outputDir <- "path/to/outputFolder/"  # this is where the model results and plots will output

source(pathToScript)

# Input File: CHiCAGO object
# required columns are: otherEndID, otherEndChr, baitID, baitChr, otherEndLen, 
#                       distSign, N.1, N.2, N.3, N, log.p, log.w, log.q, score
data = readRDS(paste0(dataDir, "chicago_input.rds"))

# Generate additional covariates based on GC and mappability
# rmap file has columns chr, start, end, rID for each restriction fragment
# baitmap file has columns chr, start, end, baitID, gene for each promoter bait
# mappability file is from UCSC portal (such as wgEncodeCrgMapabilityAlign100mer)
restrictionFeatures <- generateFeatureMap(rMapFile=paste0(dataDir, 'rmap.csv'),
                                          baitMapFile=paste0(dataDir, 'baitmap.csv'),
                                          mapFile=paste0(dataDir, "mappability.bigWig"))

#Preprocess into suitable form
processed <- ChiCAGOPreProcessing(data, restrictionFeatures, Dmax=2e6)

# Run model over selected baits (or all baits if NA)
# Specify data file, then baits to be modeled and cores
# Cores should be set at 1 for Windows
results <- ChiCModelRun(processed, baitIDs=NA, cores=60)

saveRDS(results, paste0(dataDir, "modelResults.rds"))  # save results to a file (can comment out if not desired)

# -- generate contact plots --

# read in the baitmap and iterate over it to generate the plots
baitmap <- fread(paste0(dataDir, "baitmap.csv"))
names(baitmap) <- c("chr", "start", "end", "id", "gene")
baitmap <- as.data.frame(baitmap)

for (i in 1:nrow(baitmap)) {
  # filter the results to only the specific ID
  relevantResults <- results %>% filter(baitID == baitmap[i,]$id) %>%
    select(baitChr, start, end, N, padj_fish) %>% arrange(start)
  if (nrow(relevantResults) == 0) { next }  # if no results, nothing to plot
  # create a plot with the counts and p-values on two separate scales for readability
  # use -log(pval) to make low p-values appear visually as spikes, and scale it so the spikes are actually visible; second axis is placed on the right
  ggplot(relevantResults, aes(x=relevantResults$start, y=relevantResults$N)) +
    geom_line(aes(x=relevantResults$start, y=-log10(padj_fish)*(max(relevantResults$N)/(2*max(-log10(relevantResults$padj_fish)))), color="blue"), alpha=0.75, show.legend=FALSE) +
    geom_line(show.legend=FALSE) +
    labs(x="Contact locus", y="Number of reads", title=paste0("Contacts and p-values for bait ", baitmap[i,]$gene)) +
    scale_y_continuous("Number of reads", sec.axis=sec_axis(~ . / (max(relevantResults$N)/(2*max(-log10(relevantResults$padj_fish)))), name="-log p-values"))
  ggsave(paste0(output, "ContactPlots/", baitmap[i,]$gene, ".jpg"), width=7, height=7)  # save the plot as an image
}