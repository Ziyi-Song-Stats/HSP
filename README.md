# Hierarchical Shrinkage Priors
R code for the paper:
"Clustering Computer Mouse Tracking Data with Informed Hierarchical Shrinkage Partition Priors".

Feel free to contact me via ziyis9@uci.edu if you are interested in it.

## Code 

In the folder `"Functions"`:

* The script `HSP_Rcpp_Functions.cpp` calculates the probability of a given partition under Shrinkage Partition (SP) prior and simulate a new partition via SP prior given a base partition of the same vector of items.

* The script `HSP_MCMC_R_Functions.R` contains all the other functions needed for the posterior inference using our Hierarchical Shrinkage Partition (HSP) prior.

* The script `HSP_Demo.Rmd` explains how to run the method on a simple example, which is illustrated in the followings. 

Libraries we need: 

* `Rcpp`, `RcppArmaillo`, `MASS`, `mnormt`, `MCMCprecision`, `invgamma`, `mcclust`, `mcclust.ext`.

In the folder `"Simulations"`:

* Folder `"Simulation_1a"`: "HSP_Simulation_1a" contains our codes and results on 50 replications by running the HSP method under Simulation scenario 1(a); "NoBLoC_Simulation_1a" contains results by the NoB-LoC method under scenario 1(a); "Simulation_1a_plots.Rmd" plots Figure 1 in the supplementary file.

* Folder `"Simulation_1b"`: "Simulation_1b_HSP_level10" contains our codes and replication results by running HSP under Simulation scenario 1(b) at contamination level of 10%, while "Simulation_1b_NoBLoC_level10_Results" are replication results by NoB-LoC under Simulation scenario 1(b) at contamination level of 10%. Similarly, "Simulation_1b_HSP_level20" and "Simulation_1b_NoBLoC_level20_Results" at contamination level of 20%, and "Simulation_1b_HSP_level30" and "Simulation_1b_NoBLoC_level30_Results" at contamination level of 30%. "Simulation_1b_plots.Rmd" plots Figure 2 in the supplementary file.   

* Folder `"Simulation_2"`: "HSP_Simulation_2" are our codes and 50 replications results by running the HSP under Simulation scenario 2; "HHDP_Simulation_2" are results by the alternative HHDP method under scenario 2; "Simulation_2_plots.Rmd" plots Figure 3 in the supplementary file.

* Folder `"Simulation_HSP_vs_BCPlaid"`: For the comparison of HSP and BCPlaid, "Simulation_HSP" contains our codes and replication results by running the HSP, while "Simulation_BCPlaid" are codes and results by running the BCPlaid; "Simulation_HSP_vs_BCPlaid_plots.Rmd" plots Figure 11 in the supplementary file. The bi-clustering method `BCPlaid` is implemented in R `biclust` package.


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

Then we start with simulating data in the format needed to run the HSP method. Here we generate a matrix of true labels with 60 columns and 30 rows, `labels.matrix`. We focus on the nested clusterings of rows within each column. The labels are not shared across columns. 
```{r}
set.seed(1)

J=60; I=30

labels.matrix = matrix(NA, nrow=I, ncol=J)
for (j in 1:20){
  labels.matrix[,j] = c(rep(1,5), rep(1,5), rep(2,5), rep(2,5), rep(3,5), rep(3,5))
}
for (j in 21:40){
  labels.matrix[,j] = c(rep(1,5), rep(2,5), rep(3,5), rep(1,5), rep(3,5), rep(2,5))
  }
for (j in 41:J){
  labels.matrix[,j] = c(rep(1,5), rep(2,5), rep(1,5), rep(3,5), rep(2,5), rep(3,5))
}
turbulence = 0.1 * I
for (j in 1:J){
  indx = sample(seq(1,I), size=turbulence, replace=FALSE)
  for (i in 1:turbulence){
    labels.matrix[,j][indx[i]] = sample(c(1,2,3), size=1)
    }
  labels.matrix[,j] = relabel(labels.matrix[,j])
}
```
We suppose that these 60 columns are clustered into 3 groups, i.e., { column 1-20, column 21-40, column 41-60 }, since the nested clusterings of rows are more similar (not necessarily identical) within each of these column groups rather than across the column groups. 

We print out the whole `labels.matrix` for bettering understanding.
```{r}
#labels.matrix
options(max.print = I*J)
prmatrix(labels.matrix, rowlab=rep("",I), collab=rep("",J))
#> 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 2 1 2 1
#> 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 2 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 2 1 1 1
#> 2 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1
#> 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 2 1 1 3 1 1 1 2 1 1 1
#> 1 1 2 1 1 1 2 1 2 2 1 2 1 1 1 2 1 1 2 1 2 2 1 3 2 2 2 3 2 2 3 2 3 2 2 1 2 2 2 2 2 2 2 1 2 3 2 2 2 3 2 2 3 2 2 2 3 2 3 1
#> 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 2 2 2 3 2 2 3 3 2 2 3 2 3 2 2 2 2 2 2 2 2 1 2 3 2 3 2 2 2 3 2 3 3 2 2 2 3 1 3 2
#> 1 1 2 1 1 1 1 1 1 1 1 3 1 1 1 1 1 1 2 1 2 2 2 3 2 2 2 3 2 2 3 2 3 2 2 2 2 2 2 2 1 1 2 3 2 3 2 2 2 3 2 2 3 2 2 2 3 2 3 2
#> 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 2 2 2 3 2 2 2 3 1 2 3 2 3 2 2 3 3 1 2 2 2 2 3 3 1 3 2 2 2 3 2 2 3 2 3 2 3 2 3 2
#> 1 1 2 2 1 2 3 1 1 2 1 1 1 1 1 1 1 2 3 1 2 2 2 3 2 2 2 3 2 2 3 2 3 2 2 2 2 2 1 2 2 2 2 3 2 3 3 2 2 3 2 2 3 2 2 2 3 2 3 3
#> 2 2 1 3 2 3 2 2 3 2 2 2 1 2 2 3 2 3 3 2 3 2 3 2 3 3 3 2 3 3 2 3 2 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1
#> 2 3 1 3 2 2 2 2 3 2 2 2 2 2 2 3 2 3 1 2 3 3 3 2 3 3 3 2 3 3 2 3 2 2 3 3 1 3 3 3 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 2 1 1 1
#> 2 3 1 3 2 3 2 2 3 2 2 2 2 2 1 3 2 3 1 2 3 3 3 2 3 3 3 2 3 3 2 3 2 3 3 3 3 3 3 3 1 1 1 1 3 1 1 1 1 1 1 1 1 1 1 1 3 1 1 1
#> 2 2 1 3 2 3 2 2 3 2 3 2 2 2 2 3 2 3 1 2 3 3 3 3 1 3 3 2 3 3 2 3 2 3 3 3 3 3 1 3 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 2 1 1 1
#> 2 3 1 3 2 3 2 2 3 2 2 2 2 2 2 3 2 3 1 2 3 3 3 2 3 3 1 2 3 3 2 3 1 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1
#> 2 3 1 3 2 3 2 2 3 2 2 2 2 2 2 3 2 3 1 2 1 1 1 1 1 1 1 3 1 1 1 2 1 1 1 1 1 1 1 1 3 3 3 2 3 2 3 2 3 2 3 3 2 3 3 3 1 3 1 3
#> 2 3 1 3 2 3 2 2 3 2 2 2 2 2 2 2 2 3 1 2 1 1 1 2 1 3 1 1 1 3 1 1 1 1 1 1 1 1 1 1 3 3 3 2 3 2 3 3 1 2 1 3 2 3 3 3 1 3 2 3
#> 2 3 1 3 2 3 2 2 3 2 2 2 2 2 2 3 2 3 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 3 3 3 2 3 2 3 3 3 2 2 3 2 3 3 3 1 3 2 3
#> 2 3 1 3 2 3 2 2 3 2 2 3 2 2 2 3 2 3 1 2 1 1 1 1 2 1 1 1 1 3 1 1 1 1 1 1 1 1 1 1 3 3 3 2 3 2 3 3 3 2 3 3 2 3 3 3 2 3 2 3
#> 2 3 1 3 2 3 2 2 3 2 2 2 2 2 2 3 2 3 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 3 3 3 2 3 2 1 3 3 2 3 3 1 3 3 3 1 3 1 3
#> 3 2 1 2 3 2 3 3 2 3 3 3 3 3 3 1 3 2 3 3 3 3 3 2 3 3 3 2 3 3 2 3 2 3 3 3 3 2 3 3 2 2 2 3 2 3 2 2 2 3 2 2 3 2 2 2 3 2 3 2
#> 3 2 3 2 3 2 3 3 2 3 1 3 3 3 1 2 3 2 3 3 3 3 1 2 3 3 3 2 3 3 2 3 2 3 3 3 3 3 3 3 2 2 2 3 2 3 2 2 2 3 2 2 3 2 2 2 3 2 3 2
#> 3 3 3 2 3 2 3 1 3 3 3 3 3 1 3 2 3 2 3 3 3 3 3 2 3 3 3 2 3 3 2 3 2 3 3 3 3 3 3 3 2 2 2 3 2 3 2 2 2 3 2 2 3 2 2 2 3 2 3 2
#> 3 2 3 2 3 2 3 3 2 1 3 3 3 3 3 2 3 2 3 3 1 3 3 2 3 2 3 2 3 3 2 3 2 3 3 1 3 3 3 3 2 2 2 3 2 3 2 2 2 3 2 2 3 2 2 2 3 2 3 2
#> 1 2 3 2 1 1 3 1 2 3 3 3 3 3 3 2 3 2 3 3 3 3 3 2 3 3 3 2 3 3 2 3 2 3 3 3 3 3 2 3 2 2 2 3 3 3 2 2 2 2 2 2 3 2 2 2 3 2 3 2
#> 3 2 3 2 3 2 3 3 2 3 3 3 3 3 3 2 3 2 3 2 2 1 2 3 2 2 2 3 2 2 3 2 3 2 2 2 2 2 2 2 3 3 3 2 3 2 3 3 3 2 3 3 2 3 3 3 1 3 2 3
#> 3 2 3 2 3 2 3 3 2 3 3 3 3 3 3 2 3 2 3 3 2 2 2 3 2 2 1 3 2 2 3 2 3 2 2 2 1 2 2 2 3 3 1 2 3 2 3 3 3 2 3 1 2 2 3 3 1 3 2 3
#> 3 2 3 2 3 2 3 3 2 3 3 3 3 3 1 2 3 3 3 3 2 2 1 3 2 2 2 3 2 2 3 2 3 2 2 2 2 2 2 2 3 3 3 2 3 3 3 3 2 2 3 3 2 3 2 3 1 3 2 1
#> 3 2 1 2 3 2 2 3 2 3 3 3 3 3 3 2 3 2 3 1 2 2 2 3 2 2 2 3 2 2 3 2 3 2 2 2 2 1 2 2 3 3 3 2 3 2 2 2 3 2 3 3 2 3 3 2 1 3 2 3
#> 3 2 3 2 3 2 3 3 2 3 3 3 3 3 3 2 3 2 3 3 2 2 2 3 2 2 2 3 2 2 1 2 3 2 2 2 2 2 2 2 3 3 2 2 3 2 3 3 3 2 3 2 2 3 3 3 1 3 2 3
```



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

baseline_for_c = c(rep(1,20), rep(2,20), rep(3,20))
shrink_para_for_c = rep(shrink_c, J)
crp_para_for_c = 1

baseline_for_mu = rep(1:6, each=5)
shrink_para_for_mu = rep(shrink_mu, I) 
crp_para_for_mu = 1 

shrink_para_for_pi = rep(shrink_pi, I)
crp_para_for_pi = 1
```

We fit our method and MCMC posterior samplings are stored in `hsp'.
```{r}
Iters = 5000
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

The model provides a `list` as output, containing the following elements:

+ ``partition_c_iterations``: MCMC chain for the partition of subjects (columns)
+ `mu_partitions_iterations`: MCMC chain for the base partition of conditions (rows) within each subject
+ `pi_partitions_iterations`: MCMC chain for the partition of conditions within each subject
+ `permutation_c_iterations`: MCMC chain for the permutation parameter that permutes order of subjects
+ `permutation_mus_iterations`: MCMC chain for the permutation parameter that permutes order of base partition of conditions within each subject
+ `permutation_pis_iterations`: MCMC chain for the permutation parameter that permutes otder of conditions within each subject   

## Post-processing and MCMC summary

```{r}
# Thinning
seq_thin = seq(0.2*Iters, Iters, by=1)
```

```{r}
true.col.partition = c( rep(1,20), rep(2,20), rep(3,20) )
MCMClabels.c = hsp$partition_c_iterations[seq_thin, ]       
psm.c = comp.psm(MCMClabels.c)
avg = minVI(psm.c, MCMClabels.c, method=("all"), include.greedy=TRUE)
subject_groups_number = max(avg$cl[1,])
ari_HSP = arandi(avg$cl[1,], true.col.partition)
f1measure_HSP = F1Measure(est.partition=avg$cl[1,], true.partition=true.col.partition)
```
We use `mcclust.ext::minVI` to obtain an estimated partition of columns, see it in `avg$cl[1,]`. Then we can calculate ARI and F1 measurements between the estimated partition of columns and the ground truth partition of columns, see `ari_HSP` and `f1measure_HSP`. 
```{r}
avg$cl[1,]
#> 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
ari_HSP
#> 1
f1measure_HSP
#> 1
```

Likewise, in the following, we can also use `mcclust.ext::minVI` to obtain an estimated partition of rows within each column. Then we calculate ARI's and F1's of the estimated partition of rows under each column, and average them across all the columns to use as evaluation indexes.    
```{r}
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

See `ari_HSP.subjects` and `f1measure_HSP.subjects` for the ARI's and F1's of the estimated partition of rows within each column.
```{r}
ari_HSP.subjects
#> [1] 0.7146016 0.8930258 1.0000000 0.6630845 0.8930258 0.6448297 0.7885475 0.8139524 0.8188699 0.7652962 0.7157276 0.7288461
#> [13] 0.7179530 0.7351772 0.8896193 0.7967490 0.6167401 0.7967490 1.0000000 0.6980006 1.0000000 0.8930258 0.7908955 0.6135394
#> [25] 0.6740984 0.8981704 0.8896193 1.0000000 0.6853364 0.7721718 0.9037232 0.8157371 0.5879901 0.7902071 0.8136622 0.6600213
#> [37] 0.6600213 0.7134387 0.6478329 0.7442869 0.8066667 0.7075791 1.0000000 0.8077870 0.8887752 0.7975604 0.7114575 1.0000000
#> [49] 1.0000000 0.8535956 1.0000000 0.8010671 0.9454190 0.8934529 0.8981704 0.8989819 0.8934529 0.8981704 0.8074464 0.8887752
```

```{r}
f1measure_HSP.subjects
#> [1] 0.8989821 0.9682540 1.0000000 0.8653199 0.9682540 0.8703704 0.9388889 0.9267399 0.9296296 0.8857143 0.8958548 0.8997494
#> [13] 0.9011765 0.8933333 0.9691228 0.9363636 0.8222222 0.9363636 1.0000000 0.8376068 1.0000000 0.9682540 0.9351433 0.7277679
#> [25] 0.8023504 0.9665831 0.9691228 1.0000000 0.8211153 0.8206522 0.9649123 0.8492325 0.7893248 0.8511346 0.8414855 0.8189084
#> [37] 0.8189084 0.9047619 0.8591228 0.8935574 0.9333333 0.8994709 1.0000000 0.8544212 0.9696342 0.9342161 0.8227208 1.0000000
#> [49] 1.0000000 0.8549114 1.0000000 0.8532491 0.8819444 0.9679634 0.9665831 0.9658994 0.9679634 0.9665831 0.9333333 0.9696342
```

### Codes for alternative methods
If one is interested in codes for the alternative methods, NoB-LoC and HHDP, please contact the authors Dr. Juhee Lee at UC Santa Cruz and Dr. Giovanni Rebaudo at University of Torino.


