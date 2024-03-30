#####
SimData = function(seedsSet){
  set.seed(seedsSet)
  J=60
  I=30
  
  labels.matrix = matrix(NA, nrow=I, ncol=J)
  a = c(rep(1,5), rep(1,5), rep(2,5), rep(2,5), rep(3,5), rep(3,5))
  b = c(rep(1,5), rep(2,5), rep(3,5), rep(1,5), rep(3,5), rep(2,5))
  c = c(rep(1,5), rep(2,5), rep(1,5), rep(3,5), rep(2,5), rep(3,5))
  for (j in 1:20){
    labels.matrix[,j] = a
  }
  for (j in 21:40){
    labels.matrix[,j] = b
  }
  for (j in 41:J){
    labels.matrix[,j] = c
  }
  
  turbulence = 0.2 * I
  for (j in 1:J){
    indx = sample(seq(1,I), size=turbulence, replace=FALSE)
    for (i in 1:turbulence){
      labels.matrix[,j][indx[i]] = sample(c(1,2,3), size=1)
    }
    labels.matrix[,j] = relabel(labels.matrix[,j])
  }
  
  datas = list()
  
  v = c(-0.97, 0.15, 1.37)
  for (j in 1:J){
    n.clust = max(labels.matrix[,j])
    th = sample(v, size=n.clust, replace=FALSE, prob=NULL)
    datas[[j]] = rnorm(I, mean=th[labels.matrix[,j]], sd=sqrt(0.16))
  }
  
  return(list("datas" = datas,
              "labels.matrix" = labels.matrix))
}



