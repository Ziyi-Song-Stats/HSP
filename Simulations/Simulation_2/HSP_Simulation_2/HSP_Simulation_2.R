######
######
args <- commandArgs(trailingOnly = TRUE)
seedsSet <- as.double(args[1])
file_name <- paste0("HSP_Simulation_Results_Seeds", seedsSet, ".RDS")

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
source("Simulation_2_Data_Generate.R")

simdata = SimData(seedsSet)

datas = simdata$datas
true.partition.subjects = simdata$labels.matrix

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

init_permutation_c = seq(1, J)
init_permutation_mus = init_permutation_pis = list()
for (j in 1:J){
  init_permutation_mus[[j]] = seq(1, I)
  init_permutation_pis[[j]] = seq(1, I)
}


Iters = 10000

baseline_for_c = c(rep(1,20), rep(2,20), rep(3,20))
shrink_para_for_c = rep(0, J)
crp_para_for_c = 1

baseline_for_mu = rep(1:6, each=5)
shrink_para_for_mu = rep(4, I) 
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

true.col.partition = c( rep(1,20), rep(2,20), rep(3,20) )
MCMClabels.c = hsp$partition_c_iterations[seq_thin, ]       
psm.c = comp.psm(MCMClabels.c)
avg = minVI(psm.c, MCMClabels.c, method=("all"), include.greedy=TRUE)
subject_groups_number = max(avg$cl[1,])
ari_HSP = arandi(avg$cl[1,], true.col.partition)
f1measure_HSP = F1Measure(est.partition=avg$cl[1,], true.partition=true.col.partition)


######
true.partition.subjects = simdata$labels.matrix
ari_HSP.subjects = rep(NA, J)
f1measure_HSP.subjects = rep(NA, J)
for (j in 1:J){
  MCMClabels.pi.j = hsp$pi_partitions_iterations[[j]][seq_thin, ]
  psm.pi.j = comp.psm(MCMClabels.pi.j)
  avg.j = minVI(psm.pi.j, MCMClabels.pi.j, method=("all"), include.greedy=TRUE)
  ari_HSP.subjects[j] = arandi(avg.j$cl[1,], true.partition.subjects[,j])
  f1measure_HSP.subjects[j] = F1Measure(est.partition=avg.j$cl[1,], 
                                        true.partition = true.partition.subjects[,j])
}
ari_HSP.subjects.avg = mean(ari_HSP.subjects)
f1measure_HSP.subjects.avg = mean(f1measure_HSP.subjects)


results = list("subject_groups_number" = subject_groups_number,
               "ari_HSP" = ari_HSP,
               "f1measure_HSP" = f1measure_HSP,
               "ari_HSP.subjects.avg" = ari_HSP.subjects.avg,
               "f1measure_HSP.subjects.avg" = f1measure_HSP.subjects.avg)
saveRDS(results, file_name)








