# return a list of updated mu partitions after one single iteration
# no LSP 
update_mu_partitions = function(current_mu_partitions, partition_c, pi_partitions,
                                baseline_for_mu, shrink_para_for_mu, crp_para_for_mu,
                                shrink_para_for_pi, crp_para_for_pi){
  mu_partitions = current_mu_partitions
  c = partition_c
  counts_c = as.vector(table(c))
  Nclust_c = length(counts_c)
  for (r in 1:Nclust_c){
    loc_j = which(c == r)
    loc = match(r, c)
    z = mu_partitions[[loc]]
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
        log_weight[l] = log_sp_prob(partition=z, baseline=baseline_for_mu, 
                                    crp_para=crp_para_for_mu, shrink_para=shrink_para_for_mu)
        for (j in loc_j){
          log_weight[l] = log_weight[l] + log_sp_prob(partition=pi_partitions[[j]], baseline=z, 
                                                      crp_para=crp_para_for_pi, shrink_para=shrink_para_for_pi)
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



















