library(saps)

#load("data/Breast_small.RData")
load("/Users/Daniel/Desktop/saps testing/fake.RData")
#load("/Users/Daniel/Desktop/enriched_genes.RData")

test_genesets <- c("NADERI_BREAST_CANCER_PROGNOSIS_UP",
                   "HAHTOLA_MYCOSIS_FUNGOIDES_DN",
                   "SEMBA_FHIT_TARGETS_DN",
                   "NAKAMURA_CANCER_MICROENVIRONMENT_UP",
                   "WINTER_HYPOXIA_UP")


results <- saps(geneSets, dat, time, event,
                random.samples=1000, cpus=4)

testset <- results$genesets[["NADERI_BREAST_CANCER_PROGNOSIS_UP"]]

plotKM(testset, time/365, event, x.label="Overall survival (years)")

plotRandomDensity(testset)

plotEnrichment(testset, results$rankedGenes)
