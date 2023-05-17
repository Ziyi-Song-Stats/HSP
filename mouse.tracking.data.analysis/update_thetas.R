
# update the list of thetas, and return a list of updated thetas
# Sample kernel parameters (NIG: mu, sigma)
update_thetas = function(pi_partitions, datas, mean0., lambda0., thetas){
  Ncolumn = length(pi_partitions)
  mean0.vector = mean0.
  lambda0.vector = lambda0.
  for (j in 1:Ncolumn){
    mean0. = mean0.vector[j]
    lambda0. = lambda0.vector[j]
    
    z. = pi_partitions[[j]]
    y. = datas[[j]]
    counts = as.vector(table(z.))
    L. = length(counts)
  
    mu.ast = thetas[[j]][,"mu.ast"]
    sigma.ast = thetas[[j]][,"sigma.ast"]
    
    for (l in 1:L.){
      i_l = which(z. == l)
      n_l = counts[l]
      y_l = y.[i_l]
      sigma2.ast.l = 0.3
      lambda_l = 1/(1/lambda0. + n_l/sigma2.ast.l)
      mu_l = lambda_l * (mean0./lambda0. + sum(y_l)/sigma2.ast.l) 
      mu.ast[i_l] = rnorm(n=1, mean=mu_l, sd=sqrt(lambda_l))
    }
    
    thetas[[j]] = cbind(mu.ast=mu.ast, sigma.ast=sigma.ast)
  }
  return(thetas) # return the list of updated_thetas
}





