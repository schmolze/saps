library(saps)

#load("data/Breast_small.RData")
load("/Users/Daniel/Desktop/breast_full.RData")

test_genesets <- c("SEMBA_FHIT_TARGETS_DN",
                  "CELL_DIVISION",
                  "3_5_CYCLIC_NUCLEOTIDE_PHOSPHODIESTERASE_ACTIVITY",
                  "AAACCAC,MIR-140")

results <- saps(geneSets[test_genesets,], dat, time, event, 10000, 4)

rankedGenes <- results$rankedGenes

testset <- results$genesets[["SEMBA_FHIT_TARGETS_DN"]]

plotKM(testset, time/365, event, x.label="Overall survival (years)")

plotRandomDensity(testset)

plotEnrichment(testset, results$rankedGenes)
