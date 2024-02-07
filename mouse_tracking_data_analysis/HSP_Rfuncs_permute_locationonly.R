
# update permutation
# using Metropolis-Hastings
update_permutation = function(partition, baseline, crp_para,
                              shrink_para, permutation){
  num = length(partition)
  num_shuffle = num * 1 #0.5
  index_shuffle = sample(1:num, num_shuffle, replace=FALSE)
  permutation_star = permutation    # propose a new permutation
  permutation_star[index_shuffle] = sample(permutation[index_shuffle], num_shuffle, replace=FALSE)
  ratio = log_sp_prob_permute(partition, baseline, crp_para, shrink_para, permutation_star) - 
    log_sp_prob_permute(partition, baseline, crp_para, shrink_para, permutation) 
  if (log(runif(1)) < ratio){
    permutation_update = permutation_star
  } else{
    permutation_update = permutation
  }
  return(permutation_update)
}





###### update permutation_mus
######
update_permutation_mus = function(mu_partitions, baseline_for_mu, 
                                  crp_para_for_mu, shrink_para_for_mu,
                                  permutation_mus, partition_c){
  c = partition_c
  Nclust_c = length(as.vector(table(c)))
  for (r in 1:Nclust_c){
    loc_j = which(c == r)
    loc = match(r, c)
    permutation_mu = permutation_mus[[loc]]
    
    num = length(permutation_mu)
    num_shuffle = num * 1  #0.5
    index_shuffle = sample(1:num, num_shuffle, replace=FALSE)
    permutation_mu_star = permutation_mu    # propose a new permutation
    permutation_mu_star[index_shuffle] = sample(permutation_mu[index_shuffle], num_shuffle, replace=FALSE)
    
    ratio = log_sp_prob_permute(partition=mu_partitions[[loc]], baseline=baseline_for_mu, 
                                crp_para=crp_para_for_mu, shrink_para=shrink_para_for_mu, 
                                permutation=permutation_mu_star) - 
      log_sp_prob_permute(partition=mu_partitions[[loc]], baseline=baseline_for_mu, 
                          crp_para=crp_para_for_mu, shrink_para=shrink_para_for_mu, 
                          permutation=permutation_mu)
    
    if (log(runif(1)) < ratio){
      permutation_mu_update = permutation_mu_star
    } else{
      permutation_mu_update = permutation_mu
    }
    
    for (j in loc_j){
      permutation_mus[[j]] = permutation_mu_update
    }
  }
  return(permutation_mus)
}


######
#####


# update the vector c, within a single iteration
# using Neal 8
update_c = function(current_c, baseline_for_c, shrink_para_for_c, crp_para_for_c,
                    mu_partitions, pi_partitions, shrink_para_for_pi, crp_para_for_pi,
                    baseline_for_mu, shrink_para_for_mu, crp_para_for_mu,
                    permutation_c, permutation_pis){
  z = current_c
  counts = as.vector(table(z))
  Nclust = length(counts)
  Ncolumn = length(baseline_for_c)
  
  for (j in 1:Ncolumn){
    pi_j = pi_partitions[[j]]
    c = z[j]
    counts[c] = counts[c] - 1
    if (counts[c] == 0){
      counts[c] = counts[Nclust]
      loc_z = (z == Nclust)
      z[loc_z] = c
      counts = counts[-Nclust]
      Nclust = Nclust - 1
    }
    z[j] = -1
    
    log_weight = rep(NA, Nclust+1)
    for (c in 1:Nclust){
      loc = match(c, z)
      baseline_mu = mu_partitions[[loc]]
      log_likelihood = log_sp_prob_permute(partition=pi_j, baseline=baseline_mu, 
                                   crp_para=crp_para_for_pi, shrink_para=shrink_para_for_pi,
                                   permutation = permutation_pis[[j]])
      z[j] = c
      log_prior_pmf = log_sp_prob_permute(partition=z, baseline=baseline_for_c, 
                                  crp_para=crp_para_for_c, shrink_para=shrink_para_for_c,
                                  permutation = permutation_c)
      log_weight[c] = log_prior_pmf + log_likelihood
      z[j] = -1
    }
    
    a = prior_simulate_sp_partition(baseline=baseline_for_mu, shrink_para=shrink_para_for_mu,
                                    crp_para=crp_para_for_mu)
    z[j] = Nclust + 1
    log_weight[Nclust+1] = log_sp_prob_permute(partition=z, baseline=baseline_for_c, crp_para=crp_para_for_c,
                                       shrink_para=shrink_para_for_c,permutation=permutation_c) + 
      log_sp_prob_permute(partition=pi_j, baseline=a, crp_para=crp_para_for_pi, 
                          shrink_para=shrink_para_for_pi,
                          permutation = permutation_pis[[j]])
    z[j] = -1
    
    max_weight = max(log_weight)
    log_weight = log_weight - max_weight
    loc_probs = exp(log_weight)
    loc_probs = loc_probs / sum(loc_probs)
    
    newz = sample(1:(Nclust+1), 1, replace=TRUE, prob=loc_probs)
    
    if (newz == Nclust+1){
      mu_partitions[[j]] = a
    } else{
      mu_partitions[[j]] = mu_partitions[[match(newz, z)]]
    }
    
    if (newz == Nclust+1){
      counts = c(counts, 0)
      Nclust = Nclust + 1
    }
    z[j] = newz
    counts[newz] = counts[newz] + 1
    
  }
  return(list("updated_partition_c"=z, "mu_partitions"=mu_partitions))
}


######
######


# return a list of updated mu partitions after one single iteration
# no LSP 
update_mu_partitions = function(current_mu_partitions, partition_c, pi_partitions,
                                baseline_for_mu, shrink_para_for_mu, crp_para_for_mu,
                                shrink_para_for_pi, crp_para_for_pi,
                                permutation_mus, permutation_pis){
  mu_partitions = current_mu_partitions
  c = partition_c
  counts_c = as.vector(table(c))
  Nclust_c = length(counts_c)
  for (r in 1:Nclust_c){
    loc_j = which(c == r)
    loc = match(r, c)
    z = mu_partitions[[loc]]
    
    permutation_mu = permutation_mus[[loc]]
    
    counts = as.vector(table(z))
    Nclust = length(counts)
    Nrow = length(z)
    for (n in 1:Nrow){
      l = z[n]
      counts[l] = counts[l] - 1
      if (counts[l] == 0){
        counts[l] = counts[Nclust]
        loc_z = (z == Nclust)
        z[loc_z] = l
        counts = counts[-Nclust]
        Nclust = Nclust - 1
      }
      z[n] = -1
      
      log_weight = rep(NA, Nclust+1)
      for (l in 1:(Nclust+1)){
        z[n] = l
        log_weight[l] = log_sp_prob_permute(partition=z, baseline=baseline_for_mu, 
                                            crp_para=crp_para_for_mu, 
                                            shrink_para=shrink_para_for_mu,
                                            permutation = permutation_mu)
        for (j in loc_j){
          log_weight[l] = log_weight[l] + log_sp_prob_permute(partition=pi_partitions[[j]], 
                                                              baseline=z, 
                                                              crp_para=crp_para_for_pi, 
                                                              shrink_para=shrink_para_for_pi,
                                                              permutation = permutation_pis[[j]])
        }
        z[n] = -1
      }
      
      max_weight = max(log_weight)
      log_weight = log_weight - max_weight
      loc_probs = exp(log_weight)
      loc_probs = loc_probs / sum(loc_probs)
      
      newz = sample(1:(Nclust+1), 1, replace=TRUE, prob=loc_probs)
      
      if (newz == Nclust+1){
        counts = c(counts, 0)
        Nclust = Nclust + 1
      }
      z[n] = newz
      counts[newz] = counts[newz] + 1
    }
    for (j in loc_j){
      mu_partitions[[j]] = z
    }
  }
  return(mu_partitions)
}


######
######
# update the list of pi partitions within one single iteration, return the updated list of pi partitions 
# and also return the list of thetas, thetas being changed but not updated, like modified Neal's Algorithm 8
update_pi_partitions = function(current_pi_partitions, mu_partitions, thetas,
                                shrink_para_for_pi, crp_para_for_pi, datas,
                                mean0., lambda0.,
                                permutation_pis){
  # current_pi_partitions is a list
  # mu_partitions is a list
  # thetas is a list
  # datas is a list
  
  pi_partitions = current_pi_partitions
  Ncolumn = length(pi_partitions)
  mean0.vector = mean0.
  lambda0.vector = lambda0.
  for (j in 1:Ncolumn){
    mean0. = mean0.vector[j]
    lambda0. = lambda0.vector[j]
    
    z. = pi_partitions[[j]]
    Nrow = length(z.)
    counts = as.vector(table(z.))
    Nclust = length(counts)
    y. = datas[[j]]
    mu.partition. = mu_partitions[[j]]
    theta. = thetas[[j]]
    mu. = theta.[,"mu.ast"]
    sigma. = theta.[,"sigma.ast"]
    for (n in 1:Nrow){
      c = z.[n]
      counts[c] = counts[c] - 1
      if (counts[c] == 0){
        counts[c] = counts[Nclust]
        loc_z = (z. == Nclust)
        z.[loc_z] = c
        counts = counts[-Nclust]
        Nclust = Nclust - 1
      }
      z.[n] = -1
      
      log_weights = rep(NA, Nclust+1)
      
      for (c in 1:Nclust){
        loc = match(c, z.)
        log_likelihood = dnorm(y.[n], mean=mu.[loc], sd=sigma.[loc], log=TRUE)
        z.[n] = c
        log_prior_pmf = log_sp_prob_permute(partition=z., baseline=mu.partition., 
                                            crp_para=crp_para_for_pi,
                                            shrink_para=shrink_para_for_pi,
                                            permutation = permutation_pis[[j]])
        log_weights[c] = log_prior_pmf + log_likelihood
        z.[n] = -1
      }
      
      mu0. = rnorm(1, mean=mean0., sd=sqrt(lambda0.))
      sigma0. = 0.3
      z.[n] = Nclust + 1
      log_weights[Nclust+1] = log_sp_prob_permute(partition=z., baseline=mu.partition., 
                                                  crp_para=crp_para_for_pi,
                                                  shrink_para=shrink_para_for_pi,
                                                  permutation=permutation_pis[[j]]) + 
        dnorm(y.[n], mean=mu0., sd=sigma0., log=TRUE)
      z.[n] = -1
      
      max_weight = max(log_weights)
      log_weights = log_weights - max_weight
      loc_probs = exp(log_weights)
      loc_probs = loc_probs / sum(loc_probs)
      
      newz = sample(1:(Nclust+1), 1, replace=TRUE, prob=loc_probs)
      
      if (newz == Nclust+1){
        mu.[n] = mu0.
        sigma.[n] = sigma0. 
      } else{
        mu.[n] = mu.[match(newz, z.)]
        sigma.[n] = sigma.[match(newz, z.)]
      }
      
      if (newz == Nclust+1){
        counts = c(counts, 0)
        Nclust = Nclust + 1
      }
      z.[n] = newz
      counts[newz] = counts[newz] + 1
      
    }
    pi_partitions[[j]] = z.
    thetas[[j]][,"mu.ast"] = mu.
    thetas[[j]][,"sigma.ast"] = sigma.
  }
  return(list("pi_partitions" = pi_partitions,
              "thetas" = thetas))
}



########
########
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



###########
###########

hierarchical_shrinkage_partition_permute = function(datas, init_partition_c, init_mu_partitions, init_pi_partitions,
                                            init_thetas, maxIters, baseline_for_c, shrink_para_for_c,
                                            crp_para_for_c, baseline_for_mu, shrink_para_for_mu, crp_para_for_mu,
                                            shrink_para_for_pi, crp_para_for_pi,
                                            init_permutation_c, init_permutation_mus, init_permutation_pis){
  ptm = proc.time()
  # datas is a list; init_partition_c is a vector; init_mu_partitions is a list
  # init_pi_partitions is a list; init_thetas is a list
  partition_c = init_partition_c
  mu_partitions = init_mu_partitions
  pi_partitions = init_pi_partitions
  thetas = init_thetas
  
  permutation_c = init_permutation_c
  permutation_mus = init_permutation_mus
  permutation_pis = init_permutation_pis
  
  J = length(partition_c)
  I = length(datas[[1]])
  
  # mu0 is a vector for J columns now
  mu0 = rep(NA, J)
  for (j in 1:J){
    mu0[j] = mean(datas[[j]])
  }
  # lambda0 is vector for J columns now
  lambda0 = rep(NA, J)
  for (j in 1:J){
    lambda0[j] = var(datas[[j]])
  }
  
  
  
  partition_c_iterations = matrix(NA, nrow=maxIters, ncol=J)
  
  mu_partitions_iterations = list()
  for (j in 1:J){
    mu_partitions_iterations = append(mu_partitions_iterations, list(matrix(NA, nrow=maxIters, ncol=I)))
  }
  pi_partitions_iterations = list()
  for (j in 1:J){
    pi_partitions_iterations = append(pi_partitions_iterations, list(matrix(NA, nrow=maxIters, ncol=I)))
  }
  mu.ast_iterations = list()
  for (j in 1:J){
    mu.ast_iterations = append(mu.ast_iterations, list(matrix(NA, nrow=maxIters, ncol=I)))
  }
  sigma.ast_iterations = list()
  for (j in 1:J){
    sigma.ast_iterations = append(sigma.ast_iterations, list(matrix(NA, nrow=maxIters, ncol=I)))
  }
  
  
  permutation_c_iterations = matrix(NA, nrow=maxIters, ncol=J)
  
  permutation_mus_iterations = list()
  for (j in 1:J){
    permutation_mus_iterations = append(permutation_mus_iterations, list(matrix(NA, nrow=maxIters, ncol=I)))
  }
  permutation_pis_iterations = list()
  for (j in 1:J){
    permutation_pis_iterations = append(permutation_pis_iterations, list(matrix(NA, nrow=maxIters, ncol=I)))
  }
  
  
  pb = txtProgressBar(min=0, max=maxIters, style=3)
  for (iter in 1:maxIters){
    step1 = update_c(current_c=partition_c, baseline_for_c=baseline_for_c, shrink_para_for_c=shrink_para_for_c,
                     crp_para_for_c=crp_para_for_c, mu_partitions=mu_partitions, pi_partitions=pi_partitions,
                     shrink_para_for_pi=shrink_para_for_pi, crp_para_for_pi=crp_para_for_pi,
                     baseline_for_mu=baseline_for_mu, shrink_para_for_mu=shrink_para_for_mu, 
                     crp_para_for_mu=crp_para_for_mu,
                     permutation_c=permutation_c, permutation_pis=permutation_pis)
    partition_c = step1$updated_partition_c
    mu_partitions = step1$mu_partitions
    step2 = update_mu_partitions(current_mu_partitions=mu_partitions, partition_c=partition_c,
                                 pi_partitions=pi_partitions, baseline_for_mu=baseline_for_mu,
                                 shrink_para_for_mu=shrink_para_for_mu, crp_para_for_mu=crp_para_for_mu,
                                 shrink_para_for_pi=shrink_para_for_pi, crp_para_for_pi=crp_para_for_pi,
                                 permutation_mus=permutation_mus, permutation_pis=permutation_pis)
    mu_partitions = step2
    
    # step 3
    step3 = update_pi_partitions(current_pi_partitions=pi_partitions, mu_partitions=mu_partitions, thetas=thetas,
                                 shrink_para_for_pi=shrink_para_for_pi, crp_para_for_pi=crp_para_for_pi,
                                 datas=datas, mean0.=mu0, lambda0.=lambda0,
                                 permutation_pis=permutation_pis)
    pi_partitions = step3$pi_partitions
    thetas = step3$thetas
    
    # step 4
    thetas = update_thetas(pi_partitions=pi_partitions, datas=datas, mean0.=mu0, lambda0.=lambda0,
                           thetas=thetas)
    
    # step 5: update permutation_c
    permutation_c = update_permutation(partition=partition_c, baseline=baseline_for_c, 
                                       crp_para=crp_para_for_c, shrink_para=shrink_para_for_c,
                                       permutation=permutation_c)
    
    # step 6: update permutation_mus
    permutation_mus = update_permutation_mus(mu_partitions, baseline_for_mu,
                                             crp_para_for_mu, shrink_para_for_mu,
                                             permutation_mus, partition_c)
    
    # step 7: update permutation_pis
    for (j in 1:J){
      permutation_pis[[j]] = update_permutation(partition=pi_partitions[[j]],
                                                baseline=mu_partitions[[j]],
                                                crp_para=crp_para_for_pi, shrink_para=shrink_para_for_pi,
                                                permutation=permutation_pis[[j]])
    }
    
    
    partition_c_iterations[iter, ] = partition_c
    for (j in 1:J){
      mu_partitions_iterations[[j]][iter, ] = mu_partitions[[j]]
    }
    for (j in 1:J){
      pi_partitions_iterations[[j]][iter, ] = pi_partitions[[j]]
    }
    for (j in 1:J){
      mu.ast_iterations[[j]][iter, ] = thetas[[j]][ , "mu.ast"]
    }
    for (j in 1:J){
      sigma.ast_iterations[[j]][iter, ] = thetas[[j]][ , "sigma.ast"]
    }
    
    
    permutation_c_iterations[iter, ] = permutation_c
    for (j in 1:J){
      permutation_mus_iterations[[j]][iter, ] = permutation_mus[[j]]
    }
    for (j in 1:J){
      permutation_pis_iterations[[j]][iter, ] = permutation_pis[[j]]
    }
    
    setTxtProgressBar(pb, iter)
  } # end for (iter in 1:maxIters) loop
  
  close(pb)
  proc.time() - ptm
  return(list("partition_c_iterations"=partition_c_iterations, 
              "mu_partitions_iterations"=mu_partitions_iterations,
              "pi_partitions_iterations"=pi_partitions_iterations,
              "mu.ast_iterations"=mu.ast_iterations,
              "sigma.ast_iterations"=sigma.ast_iterations,
              "permutation_c_iterations"=permutation_c_iterations,
              "permutation_mus_iterations"=permutation_mus_iterations,
              "permutation_pis_iterations"=permutation_pis_iterations
  ))
}














