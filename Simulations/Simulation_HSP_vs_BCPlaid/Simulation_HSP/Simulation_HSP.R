######
######
args <- commandArgs(trailingOnly = TRUE)
seedsSet <- as.double(args[1])
file_name <- paste0("HSP_Simulation_Results_Seeds", seedsSet, ".RDS")

############################################################
r1p4_F1Measure = function(est.partition, true.partition){
  G = max(est.partition)
  K = max(true.partition)
  
  F1A = rep(NA, G)
  for (g in 1:G){
    A = which(est.partition == g)
    a = rep(NA, K)
    for (k in 1:K){
      B = which(true.partition == k)
      precision = length(intersect(A, B)) / length(A)
      recall = length(intersect(A, B)) / length(B)
      if ( (precision+recall) == 0 ){
        a[k] = 0
      }else{
        a[k] = 2*precision*recall/(precision+recall)
      }
    }
    F1A[g] = max(a)
  }
  
  F1B = rep(NA, K)
  for (k in 1:K){
    B = which(true.partition == k)
    a = rep(NA, G)
    for (g in 1:G){
      A = which(est.partition == g)
      precision = length(intersect(B, A)) / length(B)
      recall = length(intersect(B, A)) / length(A)
      if ( (precision+recall) == 0 ){
        a[g] = 0
      }else{
        a[g] = 2*precision*recall/(precision+recall)
      }
    }
    F1B[k] = max(a)
  }
  
  return( (mean(F1A) + mean(F1B))/2 )
}
############################################################

library(mcclust)
library(mcclust.ext)
library(MASS)
library(mnormt)
library(MCMCprecision)
library(invgamma)
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("HSP_Rcpp_Functions.cpp")
source("HSP_MCMC_R_Functions.R")
source("Simulation_HSP_vs_BCPlaid_Data_Generate.R")


simdata = SimData(seedsSet)
J=10
I=30

sim.data.matrix = simdata$data
labels.matrix = simdata$labels.matrix
datas = list()
for (j in 1:J){
  datas[[j]] = sim.data.matrix[ , j]
} 

init_partition_c = rep(1,J)
init_mu_partitions = init_pi_partitions = init_thetas = list()
for (j in 1:J){
  init_mu_partitions[[j]] = rep(1,I)
  init_pi_partitions[[j]] = rep(1,I)
  init_thetas[[j]] = cbind(mu.ast=rep(mean(datas[[j]]),I), 
                           sigma.ast=rep(sd(datas[[j]]),I) )
}

init_permutation_c = seq(1, J)
init_permutation_mus = init_permutation_pis = list()
for (j in 1:J){
  init_permutation_mus[[j]] = seq(1, I)
  init_permutation_pis[[j]] = seq(1, I)
}


Iters = 10000

baseline_for_c = rep(1,10)
shrink_para_for_c = rep(0, J)
crp_para_for_c = 1

baseline_for_mu = rep(1:6, each=5)
shrink_para_for_mu = rep(0, I) 
crp_para_for_mu = 1 

shrink_para_for_pi = rep(3.5, I)
crp_para_for_pi = 1


hsp = hierarchical_shrinkage_partition_permute(datas=datas, init_partition_c=init_partition_c, 
                                               init_mu_partitions=init_mu_partitions,
                                               init_pi_partitions=init_pi_partitions,
                                               init_thetas=init_thetas, maxIters=Iters,
                                               baseline_for_c=baseline_for_c,
                                               shrink_para_for_c=shrink_para_for_c, crp_para_for_c=crp_para_for_c,
                                               baseline_for_mu=baseline_for_mu, shrink_para_for_mu=shrink_para_for_mu,
                                               crp_para_for_mu=crp_para_for_mu, shrink_para_for_pi=shrink_para_for_pi,
                                               crp_para_for_pi=crp_para_for_pi,
                                               init_permutation_c=init_permutation_c, 
                                               init_permutation_mus=init_permutation_mus,
                                               init_permutation_pis=init_permutation_pis)


# Thinning
R = 10000
seq_thin = seq(0.2*R, R, by=1)
library(mcclust)
library(mcclust.ext)

HSP.results = matrix(NA, nrow=I, ncol=J)
for (j in 1:J){
  MCMClabels.pi.j = hsp$pi_partitions_iterations[[j]][seq_thin, ]
  psm.pi.j = comp.psm(MCMClabels.pi.j)
  avg.j = minVI(psm.pi.j, MCMClabels.pi.j, method=("all"), include.greedy=TRUE)
  HSP.results[ , j] = relabel(avg.j$cl[1,])
}

F1Measure = r1p4_F1Measure(est.partition = as.vector(t(HSP.results)),
                           true.partition = as.vector(t(labels.matrix)))

saveRDS(F1Measure, file_name)










