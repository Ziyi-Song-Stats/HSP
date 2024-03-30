#####
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

