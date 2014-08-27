library(saps)

load("data/Breast_small.RData")

results <- saps(geneSets[1:3,], dat, time, event, 25, 4)

testset <- results$genesets[["NAKAMURA_CANCER_MICROENVIRONMENT_UP"]]

cluster <- testset$cluster
p_pure <- testset$saps_unadjusted["p_pure"]
p_pure_adj <- testset$saps_adjusted["p_pure"]
p_random <- testset$saps_unadjusted["p_random"]

dd <- data.frame("time"=time, "event"=event, "cluster"=cluster)

text <- paste("p_pure = ", round(p_pure, digits=3))

km.coxph.plot(formula.s = Surv(time/365, event) ~ cluster, data.s = dd,
              main.title="KM Plot for P_pure",
              y.label="Probability of survival",
              x.label="Overall Survival (years)",
              leg.text=paste(c("Group 1", "Group 2")),
              .lwd=2,
              o.text=text)
