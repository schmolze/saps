library(saps)

load("data/Breast_small.RData")

results <- saps(geneSets[1:30,], dat, time, event, 25, 4)
