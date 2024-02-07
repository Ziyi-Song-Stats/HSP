#####
#####
#####
relabel = function(x){
  y = double(length(x))
  unique.values = unique(x)
  L = length(unique.values)
  for ( j in 1:L){
    index = which(x == unique.values[j])
    y[index] = j
  }
  return(y)
}



