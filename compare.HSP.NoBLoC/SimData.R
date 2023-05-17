#####
SimData = function(seedsSet){
  set.seed(seedsSet)
  J=60
  #I=12
  I=30
  
  labels.matrix = matrix(NA, nrow=I, ncol=J)
  a = c(rep(1,5), rep(1,5), rep(2,5), rep(2,5), rep(3,5), rep(3,5))
  b = c(rep(1,5), rep(2,5), rep(3,5), rep(1,5), rep(3,5), rep(2,5))
  c = c(rep(1,5), rep(2,5), rep(1,5), rep(3,5), rep(2,5), rep(3,5))
  #a = c(rep(1,2), rep(1,2), rep(2,2), rep(2,2), rep(3,2), rep(3,2))
  #b = c(rep(1,2), rep(2,2), rep(3,2), rep(1,2), rep(3,2), rep(2,2))
  #c = c(rep(1,2), rep(2,2), rep(1,2), rep(3,2), rep(2,2), rep(3,2))
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
  #v = c(-2.5, 0, 2.5)
  #v = c(-3, 0, 3)
  #v = c(-5, 0, 5)
  v = c(-3.5, 0, 3.5)
  for (j in 1:J){
    n.clust = max(labels.matrix[,j])
    th = sample(v, size=n.clust, replace=FALSE, prob=NULL)
    # datas[[j]] = rnorm(I, mean=th[labels.matrix[,j]], sd=sqrt(0.5))
    datas[[j]] = rnorm(I, mean=th[labels.matrix[,j]], sd=sqrt(1))
  }

  return(datas)
}







