library(saps)

#load("data/Breast_small.RData")
load("/Users/Daniel/Desktop/breast_full.RData")
#load("/Users/Daniel/Desktop/enriched_genes.RData")

#results <- saps(geneSets[1:100,], dat, time, event, cpus=4)

#testset <- results$genesets[["random1"]]

plotKM(testset, time/365, event, x.label="Overall survival (years)")

plotRandomDensity(testset)

plotEnrichment(testset, results$rankedGenes)
