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







