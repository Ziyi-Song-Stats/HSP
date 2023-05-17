hierarchical_shrinkage_partition = function(datas, init_partition_c, init_mu_partitions, init_pi_partitions,
                                            init_thetas, maxIters, baseline_for_c, shrink_para_for_c,
                                            crp_para_for_c, baseline_for_mu, shrink_para_for_mu, crp_para_for_mu,
                                            shrink_para_for_pi, crp_para_for_pi){
  ptm = proc.time()
  # datas is a list; init_partition_c is a vector; init_mu_partitions is a list
  # init_pi_partitions is a list; init_thetas is a list
  partition_c = init_partition_c
  mu_partitions = init_mu_partitions
  pi_partitions = init_pi_partitions
  thetas = init_thetas
  
  J = length(partition_c)
  I = length(datas[[1]])
  
  S0 = 2 #3 # shape parameter
  s0 = 1 # rate parameter
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
  
  pb = txtProgressBar(min=0, max=maxIters, style=3)
  for (iter in 1:maxIters){
    step1 = update_c(current_c=partition_c, baseline_for_c=baseline_for_c, shrink_para_for_c=shrink_para_for_c,
                     crp_para_for_c=crp_para_for_c, mu_partitions=mu_partitions, pi_partitions=pi_partitions,
                     shrink_para_for_pi=shrink_para_for_pi, crp_para_for_pi=crp_para_for_pi,
                     baseline_for_mu=baseline_for_mu, shrink_para_for_mu=shrink_para_for_mu, 
                     crp_para_for_mu=crp_para_for_mu)
    partition_c = step1$updated_partition_c
    mu_partitions = step1$mu_partitions
    step2 = update_mu_partitions(current_mu_partitions=mu_partitions, partition_c=partition_c,
                                 pi_partitions=pi_partitions, baseline_for_mu=baseline_for_mu,
                                 shrink_para_for_mu=shrink_para_for_mu, crp_para_for_mu=crp_para_for_mu,
                                 shrink_para_for_pi=shrink_para_for_pi, crp_para_for_pi=crp_para_for_pi)
    mu_partitions = step2
    
    # step 3
    step3 = update_pi_partitions(current_pi_partitions=pi_partitions, mu_partitions=mu_partitions, thetas=thetas,
                                         shrink_para_for_pi=shrink_para_for_pi, crp_para_for_pi=crp_para_for_pi,
                                         datas=datas, mean0.=mu0, lambda0.=lambda0, s0.=s0, S0.=S0)
    pi_partitions = step3$pi_partitions
    thetas = step3$thetas
    
    # step 4
    thetas = update_thetas(pi_partitions=pi_partitions, datas=datas, mean0.=mu0, lambda0.=lambda0, s0.=s0, S0.=S0,
                           thetas=thetas)

    
    partition_c_iterations[iter, ] = partition_c
    for (j in 1:J){
      mu_partitions_iterations[[j]][iter, ] = mu_partitions[[j]]
    }
    for (j in 1:J){
      pi_partitions_iterations[[j]][iter, ] = pi_partitions[[j]]
    }
    
    setTxtProgressBar(pb, iter)
  } # end for (iter in 1:maxIters) loop
  
  close(pb)
  proc.time() - ptm
  return(list("partition_c_iterations"=partition_c_iterations, 
              "mu_partitions_iterations"=mu_partitions_iterations,
              "pi_partitions_iterations"=pi_partitions_iterations
              ))
}












