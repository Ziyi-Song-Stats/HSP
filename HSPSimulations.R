######
######
args <- commandArgs(trailingOnly = TRUE)
seedsSet <- as.double(args[1])
file_name <- paste0("HSP_Simulation_Results_Seeds", seedsSet, ".RDS")


source("SimData.R")
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("log_sp_prob.cpp")
Rcpp::sourceCpp("prior_simulate_sp_partition.cpp")
source("update_c.R")
source("update_mu_partitions.R")
source("update_pi_partitions.R")
source("update_thetas.R")
source("hierarchical_shrinkage_partition.R")
library(MASS)
library(mnormt)
library(MCMCprecision)
library(invgamma)


datas = SimData(seedsSet)
J=60
I=30

init_partition_c = rep(1,J)
init_mu_partitions = init_pi_partitions = init_thetas = list()
for (j in 1:J){
  init_mu_partitions[[j]] = rep(1,I)
  init_pi_partitions[[j]] = rep(1,I)
  init_thetas[[j]] = cbind(mu.ast=rep(mean(datas[[j]]),I), 
                           sigma.ast=rep(sd(datas[[j]]),I) )
}


Iters = 10000

baseline_for_c = c(rep(1,20), rep(2,20), rep(3,20))
shrink_para_for_c = rep(0, J)
crp_para_for_c = 1

baseline_for_mu = rep(1:6, each=5)
shrink_para_for_mu = rep(3, I) 
crp_para_for_mu = 1 

shrink_para_for_pi = rep(3, I)
crp_para_for_pi = 1


hsp = hierarchical_shrinkage_partition(datas=datas, init_partition_c=init_partition_c, 
                                       init_mu_partitions=init_mu_partitions,
                                       init_pi_partitions=init_pi_partitions,
                                       init_thetas=init_thetas, maxIters=Iters,
                                       baseline_for_c=baseline_for_c,
                                       shrink_para_for_c=shrink_para_for_c, crp_para_for_c=crp_para_for_c,
                                       baseline_for_mu=baseline_for_mu, shrink_para_for_mu=shrink_para_for_mu,
                                       crp_para_for_mu=crp_para_for_mu, shrink_para_for_pi=shrink_para_for_pi,
                                       crp_para_for_pi=crp_para_for_pi)

saveRDS(hsp, file_name)


