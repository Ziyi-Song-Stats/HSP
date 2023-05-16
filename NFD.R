

## compute the Normalized Frobenius Distance
NFD = function(A, B){
  p = nrow(A)
  nfd = 0
  for (i in 1:p){
    for (j in 1:p){
      nfd = nfd + (A[i,j] - B[i,j])^2
    }
  }
  nfd = nfd / (p^2)
  return(nfd)
}







