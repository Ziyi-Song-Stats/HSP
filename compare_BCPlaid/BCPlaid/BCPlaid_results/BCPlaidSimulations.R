######
######
args <- commandArgs(trailingOnly = TRUE)
seedsSet <- as.double(args[1])
file_name <- paste0("BCPlaid_Simulation_Results_Seeds", seedsSet, ".RDS")


source("SimData.R")
source("r1p4_F1Measure.R")
library(biclust)

simdata = SimData(seedsSet)

data = simdata$data
labels.matrix = simdata$labels.matrix
  
J=10
I=30

BCPlaid.res = biclust(data, method=BCPlaid(), cluster="r",
                      fit.model=y~m, 
                      shuffle=100, iter.startup=10, iter.layer=100, verbose=FALSE)
BCPlaid.bicluster.res = biclusternumber(BCPlaid.res)

BCPlaid.results = matrix(0, nrow=I, ncol=J)
for (i in 1:length(BCPlaid.bicluster.res)){
  rows.idx = BCPlaid.bicluster.res[[i]]$Rows
  cols.idx = BCPlaid.bicluster.res[[i]]$Cols
  BCPlaid.results[rows.idx, cols.idx] = i
}

F1Measure = r1p4_F1Measure(est.partition = as.vector(t(BCPlaid.results)),
                           true.partition = as.vector(t(labels.matrix)))

saveRDS(F1Measure, file_name)










