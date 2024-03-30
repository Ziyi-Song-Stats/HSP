

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






