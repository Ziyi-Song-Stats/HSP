####
###
# calculate a symmetrized version of the F1-measure

F1Measure = function(est.partition, true.partition){
  est.unique = unique(est.partition)
  G = length(est.unique)
  true.unique = unique(true.partition)
  K = length(true.unique)
  
  F1A = rep(NA, G)
  for (g in 1:G){
    A = which(est.partition == est.unique[g])
    a = rep(NA, K)
    for (k in 1:K){
      B = which(true.partition == true.unique[k])
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
    B = which(true.partition == true.unique[k])
    a = rep(NA, G)
    for (g in 1:G){
      A = which(est.partition == est.unique[g])
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











