# Hierarchical Shrinkage Priors
R code for the paper:
"Clustering Computer Mouse Tracking Data with Informed Hierarchical Shrinkage Partition Priors".

## Code 

The R scripts in the folder "Functions" are functions needed to run our HSP method. 

* The script "HSP_Rcpp_Functions.cpp" calculates the probability of a given partition under Shrinkage Partition (SP) prior and simulate a new partition via SP prior given a base partition of the same vector of items.

* The script "HSP_MCMC_R_Functions.R" contain all the other functions needed for the posterior inference using our Hierarchical Shrinkage Partition (HSP) prior.

Libraries: Rcpp, RcppArmaillo, MASS, mnormt, MCMCprecision, invgamma. 

# Example: instructions for use

We now explain how to run our HSP method via a simple simulation example. 

We first set up all the required files and libararies.
```{r}
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

Then we start with simulating data in the format needed to run the HSP method. For simplicity, we generate a small matrix of true labels with 12 columns and 18 rows. We focus on the nested clusterings of rows within each column. The labels are not shared across columns. 
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
We suppose that the 12 columns are clustered into 3 groups, i.e., { column 1-4, column 5-8, column 9-12 }, such that nested clusterings of rows are more similar (not necessarily identical) within each column group rather than across the column groups.  


## Running HSP

sdf

## Post-processing

sdf

Feel free to download the zip file. Codes used for the paper should be in it.

The folders Simulation 1a, Simulation 1b, and Simulation 2 correspond to the simulation scenarios in the paper. Each folder contains the needed codes and simulated data and results. You can also see a Description text file in each folder that explains details about the folder. 

The folder compare_BCPlaid provides codes to compare with the bi-clustering method BCPlaid, which is implemented in biclust R package. See details in our Supplementary Material.


The folder "Functions" contains the main functions implemented for our HSP method. Followings are brief explanations about the functions.

"SimData.R": Data simulation. It might be different across simulation scenarios. See scenarios in the paper. 

"relabel.R": Given labels of a vector of items, it relabel the items using incremental label numbers, while the clustering unchanged. For example, {3,3,1,2,2} --> {1,1,2,3,3}.

"F1Measure.R": Calculate a symmetrized version of the F1-measure between two clusterings on the same items.

"HSP_Rcppfuncs_permute.cpp": log_sp_prob_permute() computes the log probability density of a partition given the base partition and the baseline CRP distribution under Shrinkage Partition distribution. prior_simulate_sp_partition() generates a new partition given a base partition under Shrinkage Partition distribution. See details in the paper. This file is written in Rcpp for speeding up. 

"HSP_Rfuncs_permute.R": Functions for MCMC samplings. See MCMC details in our Supplementary Material. The function hierarchical_shrinkage_partition_permute() outputs a list of MCMC iteration results, which include subject partition, condition partition within each subject, permutation parameters, etc., for each iteration. 

"HSPSimulations.R": Given simulated data, informed base partitions of subjects and conditions, and shrinkage parameter values and others, it runs our method and put mcmc results in a list, called hsp. Then, it calculate ARI and F1 measurements between the estimated partitions and ground truth or the given base partitions. You can change the shrinkage parameter values by yourself when you run it. See details in our paper.



Feel free to contact me via ziyis9@uci.edu if you are interested in it.




