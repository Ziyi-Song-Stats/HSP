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

Then we start with simulating data in the format needed to run the HSP method. We generate a matrix of true labels with 30 columns and 18 rows. We focus on the nested clusterings of rows within each column. The labels are not shared across columns. 

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




