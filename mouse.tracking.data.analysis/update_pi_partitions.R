
# update the list of pi partitions within one single iteration, return the updated list of pi partitions 
# and also return the list of thetas, thetas being changed but not updated, like modified Neal's Algorithm 8
update_pi_partitions = function(current_pi_partitions, mu_partitions, thetas,
                                shrink_para_for_pi, crp_para_for_pi, datas,
                                mean0., lambda0.){
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
        log_prior_pmf = log_sp_prob(partition=z., baseline=mu.partition., crp_para=crp_para_for_pi,
                                    shrink_para=shrink_para_for_pi)
        log_weights[c] = log_prior_pmf + log_likelihood
        z.[n] = -1
      }

      mu0. = rnorm(1, mean=mean0., sd=sqrt(lambda0.))
      sigma0. = 0.3
      z.[n] = Nclust + 1
      log_weights[Nclust+1] = log_sp_prob(partition=z., baseline=mu.partition., crp_para=crp_para_for_pi,
                                          shrink_para=shrink_para_for_pi) + dnorm(y.[n], mean=mu0., 
                                                                                  sd=sigma0., log=TRUE)
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
  







