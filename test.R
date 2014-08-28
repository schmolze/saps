library(saps)

load("data/Breast_small.RData")

results <- saps(geneSets[1:3,], dat, time, event, 1000, 4)

ci <- results$rankedGenes

testset <- results$genesets[["NAKAMURA_CANCER_MICROENVIRONMENT_DN"]]
genes <- testset[["genes"]]

gene_vals <- ci[genes]

plotKM(testset, time/365, event, x.label="Overall survival (years)")

plotRandomDensity(testset)

#stripchart(wing[type=="apf"], pch=1, method="stack", main="Wing", xlim=range(wing), col="blue")
#stripchart(wing[type=="af"],pch=2,method="stack",add=T, col="red")

stripchart(ci, pch=1, main="Main Title", xlab="title")
stripchart(gene_vals, pch="*", add=T, col="red")
