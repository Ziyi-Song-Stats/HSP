######
######
args <- commandArgs(trailingOnly = TRUE)
seedsSet <- as.double(args[1])
file_name <- paste0("BCPlaid_Simulation_Results_Seeds", seedsSet, ".RDS")


library(biclust)

##############################################################
##############################################################
r1p4_F1Measure = function(est.partition, true.partition){
  G = max(est.partition)
  K = max(true.partition)
  
  F1A = rep(NA, G)
  for (g in 1:G){
    A = which(est.partition == g)
    a = rep(NA, K)
    for (k in 1:K){
      B = which(true.partition == k)
      precision = length(intersect(A, B)) / length(A)
      recall = length(intersect(A, B)) / length(B)
      if ( (precision+recall) == 0 ){
        a[k] = 0
      }else{
        a[k] = 2*precision*recall/(precision+recall)
      }
    }
    F1A[g] = max(a)
  }
  
  F1B = rep(NA, K)
  for (k in 1:K){
    B = which(true.partition == k)
    a = rep(NA, G)
    for (g in 1:G){
      A = which(est.partition == g)
      precision = length(intersect(B, A)) / length(B)
      recall = length(intersect(B, A)) / length(A)
      if ( (precision+recall) == 0 ){
        a[g] = 0
      }else{
        a[g] = 2*precision*recall/(precision+recall)
      }
    }
    F1B[k] = max(a)
  }
  
  return( (mean(F1A) + mean(F1B))/2 )
}
##############################################################
##############################################################
SimData = function(seedsSet){
  set.seed(seedsSet)
  J=10
  I=30
  labels.matrix = matrix(NA, nrow=I, ncol=J)
  for (j in 1:J){
    labels.matrix[ , j] = c(rep(1,10), rep(2,10), rep(3,10))
  }
  data = matrix(NA, nrow=I, ncol=J)
  data[labels.matrix == 1] = rnorm(I*J/3, -0.97, 0.16)
  data[labels.matrix == 2] = rnorm(I*J/3, 0.15, 0.16)
  data[labels.matrix == 3] = rnorm(I*J/3, 1.37, 0.16)
  
  return(list("data" = data,
              "labels.matrix" = labels.matrix))
}
##############################################################
##############################################################

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










