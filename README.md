# Hierarchical Shrinkage Priors
R code for the paper:
"Clustering Computer Mouse Tracking Data with Informed Hierarchical Shrinkage Partition Priors".

Feel free to contact me via ziyis9@uci.edu if you are interested in it.

## Code 

The R scripts in the folder `Functions` are functions needed to run our HSP method. 

* The script `HSP_Rcpp_Functions.cpp` calculates the probability of a given partition under Shrinkage Partition (SP) prior and simulate a new partition via SP prior given a base partition of the same vector of items.

* The script `HSP_MCMC_R_Functions.R` contain all the other functions needed for the posterior inference using our Hierarchical Shrinkage Partition (HSP) prior.

Libraries: Rcpp, RcppArmaillo, MASS, mnormt, MCMCprecision, invgamma, mcclust, mcclust.ext. 

# Example: instructions for use

We now explain how to run our HSP method via a simple simulation example. 

We first set up all the required files and libararies.
```{r}
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
```

## Simulating data

Then we start with simulating data in the format needed to run the HSP method. For simplicity, we generate a small matrix of true labels with 12 columns and 18 rows, `labels.matrix`. We focus on the nested clusterings of rows within each column. The labels are not shared across columns. 
``` r
J=12; I=18

labels.matrix = matrix(NA, nrow=I, ncol=J)
for (j in 1:4){
  labels.matrix[,j] = c(rep(1,3), rep(1,3), rep(2,3), rep(2,3), rep(3,3), rep(3,3))
}
for (j in 5:8){
  labels.matrix[,j] = c(rep(1,3), rep(2,3), rep(3,3), rep(1,3), rep(3,3), rep(2,3))
  }
for (j in 9:J){
  labels.matrix[,j] = c(rep(1,3), rep(2,3), rep(1,3), rep(3,3), rep(2,3), rep(3,3))
}
turbulence = 0.1 * I
for (j in 1:J){
  indx = sample(seq(1,I), size=turbulence, replace=FALSE)
  for (i in 1:turbulence){
    labels.matrix[,j][indx[i]] = sample(c(1,2,3), size=1)
    }
  labels.matrix[,j] = relabel(labels.matrix[,j])
}

labels.matrix
#>
#>        [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
#> [1,]    1    1    1    1    1    1    1    1    1     1     1     1
#> [2,]    1    1    1    1    1    1    1    1    1     2     1     1
#> [3,]    2    1    1    1    1    1    1    1    1     1     1     1
#> [4,]    1    1    1    1    2    2    2    2    2     2     2     1
#> [5,]    1    1    1    1    2    2    2    2    2     2     3     2
#> [6,]    1    1    1    1    2    2    2    2    2     2     3     2
#> [7,]    3    2    2    2    3    3    3    3    1     1     1     1
#> [8,]    3    1    2    2    3    3    3    3    1     1     1     1
#> [9,]    3    2    2    2    3    3    3    3    1     1     1     1
#> [10,]   3    2    2    2    1    1    1    1    3     3     2     3
#> [11,]   3    2    2    2    1    1    1    1    3     3     2     3
#> [12,]   3    2    2    2    1    1    1    1    1     3     2     3
#> [13,]   2    3    3    3    3    3    3    3    2     2     3     2
#> [14,]   2    3    3    3    3    3    3    3    2     2     3     2
#> [15,]   2    3    3    3    2    3    3    3    2     2     3     2
#> [16,]   2    3    3    3    2    3    2    2    3     3     2     3
#> [17,]   2    3    3    3    2    2    2    2    3     3     2     3
#> [18,]   2    3    3    3    2    2    2    2    3     3     2     3
```
We suppose that the above 12 columns are clustered into 3 groups, i.e., { column 1-4, column 5-8, column 9-12 }, since the nested clusterings of rows are more similar (not necessarily identical) within each of these column groups rather than across the column groups. 

According to the matrix of true labels, we simulate data points correspondingly stored in `datas` in the format of `List` in `R`.

```{r}
datas = list()
v = c(-0.97, 0.15, 1.37)
for (j in 1:J){
  n.clust = max(labels.matrix[,j])
  th = sample(v, size=n.clust, replace=FALSE, prob=NULL)
  datas[[j]] = rnorm(I, mean=th[labels.matrix[,j]], sd=sqrt(0.16))
}
```


## Running HSP

Before fitting our method, we need to assign initial values for the partition of columns, base partitions for rows, partitions of rows, all the kernel parameters, and all the permutattion parameters, for the sake of MCMC samplings. Users can choose their own initial values based on the domain knowledge about their real data. 
```{r}
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
```

We also need to specify base partition of columns, base partition of rows, and values for the three shrinkage parameters. See details in the paper. 

Here we let $\boldsymbol{\tau} = \boldsymbol{\rho} = 0$ and $\boldsymbol{\lambda} = 3.5$, without leveraging any base partitions of columns or rows provided by prior knowledge. 
```{r}
shrink_c = 0; shrink_mu = 0; shrink_pi = 3.5

baseline_for_c = c(rep(1,4), rep(2,4), rep(3,4))
shrink_para_for_c = rep(shrink_c, J)
crp_para_for_c = 1

baseline_for_mu = rep(1:6, each=3)
shrink_para_for_mu = rep(shrink_mu, I) 
crp_para_for_mu = 1 

shrink_para_for_pi = rep(shrink_pi, I)
crp_para_for_pi = 1
```

We fit our method and MCMC posterior samplings are stored in `hsp'.
```{r}
Iters = 10000
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
```

## Post-processing and MCMC summary

```{r}
# Thinning
seq_thin = seq(0.2*Iters, Iters, by=1)

true.col.partition = c( rep(1,4), rep(2,4), rep(3,4) )
MCMClabels.c = hsp$partition_c_iterations[seq_thin, ]       
psm.c = comp.psm(MCMClabels.c)
avg = minVI(psm.c, MCMClabels.c, method=("all"), include.greedy=TRUE)
subject_groups_number = max(avg$cl[1,])
ari_HSP = arandi(avg$cl[1,], true.col.partition)
f1measure_HSP = F1Measure(est.partition=avg$cl[1,], true.partition=true.col.partition)

ari_HSP.subjects = rep(NA, J)
f1measure_HSP.subjects = rep(NA, J)
for (j in 1:J){
  MCMClabels.pi.j = hsp$pi_partitions_iterations[[j]][seq_thin, ]
  psm.pi.j = comp.psm(MCMClabels.pi.j)
  avg.j = minVI(psm.pi.j, MCMClabels.pi.j, method=("all"), include.greedy=TRUE)
  ari_HSP.subjects[j] = arandi(avg.j$cl[1,], labels.matrix[,j])
  f1measure_HSP.subjects[j] = F1Measure(est.partition=avg.j$cl[1,], 
                                        true.partition = labels.matrix[,j])
}
ari_HSP.subjects.avg = mean(ari_HSP.subjects)
f1measure_HSP.subjects.avg = mean(f1measure_HSP.subjects)


results = list("subject_groups_number" = subject_groups_number,
               "ari_HSP" = ari_HSP,
               "f1measure_HSP" = f1measure_HSP,
               "ari_HSP.subjects.avg" = ari_HSP.subjects.avg,
               "f1measure_HSP.subjects.avg" = f1measure_HSP.subjects.avg)
```



The folders Simulation 1a, Simulation 1b, and Simulation 2 correspond to the simulation scenarios in the paper. Each folder contains the needed codes and simulated data and results. You can also see a Description text file in each folder that explains details about the folder. 

The folder compare_BCPlaid provides codes to compare with the bi-clustering method BCPlaid, which is implemented in biclust R package. See details in our Supplementary Material.



"HSPSimulations.R": Given simulated data, informed base partitions of subjects and conditions, and shrinkage parameter values and others, it runs our method and put mcmc results in a list, called hsp. Then, it calculate ARI and F1 measurements between the estimated partitions and ground truth or the given base partitions. You can change the shrinkage parameter values by yourself when you run it. See details in our paper.
