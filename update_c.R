# update the vector c, within a single iteration
# using Neal 8
update_c = function(current_c, baseline_for_c, shrink_para_for_c, crp_para_for_c,
                    mu_partitions, pi_partitions, shrink_para_for_pi, crp_para_for_pi,
                    baseline_for_mu, shrink_para_for_mu, crp_para_for_mu){
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
      log_likelihood = log_sp_prob(partition=pi_j, baseline=baseline_mu, 
                           crp_para=crp_para_for_pi, shrink_para=shrink_para_for_pi)
      z[j] = c
      log_prior_pmf = log_sp_prob(partition=z, baseline=baseline_for_c, 
                          crp_para=crp_para_for_c, shrink_para=shrink_para_for_c)
      log_weight[c] = log_prior_pmf + log_likelihood
      z[j] = -1
    }
    
    a = prior_simulate_sp_partition(baseline=baseline_for_mu, shrink_para=shrink_para_for_mu,
                                    crp_para=crp_para_for_mu)
    z[j] = Nclust + 1
    log_weight[Nclust+1] = log_sp_prob(partition=z, baseline=baseline_for_c, crp_para=crp_para_for_c,
                               shrink_para=shrink_para_for_c) + 
      log_sp_prob(partition=pi_j, baseline=a, crp_para=crp_para_for_pi, shrink_para=shrink_para_for_pi)
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




















