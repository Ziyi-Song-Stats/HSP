#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;


// [[Rcpp::export]]
double log_sp_prob(arma::colvec partition, arma::colvec baseline, double crp_para, arma::colvec shrink_para){
  int label, base, Nclust, denom;
  arma::colvec pre_labels, pre_bases, unique_labels, d;
  
  int N = partition.size();
  arma::colvec logpmf(N-1);
  
  for (int i=1; i<N; i++){
    label = partition(i);
    pre_labels = partition(arma::span(0,i-1));
    base = baseline(i);
    pre_bases = baseline(arma::span(0,i-1));
    
    unique_labels = unique(pre_labels);
    Nclust = unique_labels.size();
    denom = accu(pre_bases == base);
    
    arma::colvec b(Nclust+1);
    arma::colvec a(Nclust+1);
    
    for (int j=0; j<Nclust; j++){
      a(j) = accu(pre_labels == unique_labels(j));
    }
    a(Nclust) = crp_para;
    a = a / (crp_para + i);
    
    if (denom == 0){
      b(Nclust) = 1;
    }
    else{
      for (int c=0; c<Nclust; c++){
        b(c) = accu((pre_labels == unique_labels(c)) % (pre_bases == base)) / denom;
      }
    }
    b = exp(shrink_para(i) * b);
    
    d = a % b;
    d = d / accu(d);
    
    arma::uvec loc = find(unique_labels == label);
    if (loc.n_elem == 0){
      logpmf(i-1) = log(d(Nclust));
    } else{
      logpmf(i-1) = log(d(loc(0)));
    }
    
  }
  return accu(logpmf);
}


////////////
/////////////


// [[Rcpp::export]]
double log_sp_prob_permute(arma::colvec partition, 
                           arma::colvec baseline, 
                           double crp_para, 
                           arma::colvec shrink_para,
                           arma::colvec permutation){
  // Convert the permutation to 0-based indexing
  arma::uvec zero_based_permu = arma::conv_to<arma::uvec>::from(permutation - 1);  //very important
  //   // Apply the permutation
  partition = partition(zero_based_permu);  //very important
  baseline = baseline(zero_based_permu);  //very important
  
  int label, base, Nclust, denom;
  arma::colvec pre_labels, pre_bases, unique_labels, d;
  
  int N = partition.size();
  arma::colvec logpmf(N-1);
  
  for (int i=1; i<N; i++){
    label = partition(i);
    pre_labels = partition(arma::span(0,i-1));
    base = baseline(i);
    pre_bases = baseline(arma::span(0,i-1));
    
    unique_labels = unique(pre_labels);
    Nclust = unique_labels.size();
    denom = accu(pre_bases == base);
    
    arma::colvec b(Nclust+1);
    arma::colvec a(Nclust+1);
    
    for (int j=0; j<Nclust; j++){
      a(j) = accu(pre_labels == unique_labels(j));
    }
    a(Nclust) = crp_para;
    a = a / (crp_para + i);
    
    if (denom == 0){
      b(Nclust) = 1;
    } 
    else{
      for (int c=0; c<Nclust; c++){
        b(c) = accu((pre_labels == unique_labels(c)) % (pre_bases == base)) / denom;
      }
    } 
    b = exp(shrink_para(i) * b);
    
    d = a % b;
    d = d / accu(d);
    
    arma::uvec loc = find(unique_labels == label);
    if (loc.n_elem == 0){
      logpmf(i-1) = log(d(Nclust));
    } else{ 
      logpmf(i-1) = log(d(loc(0)));
    } 
    
  } 
  return accu(logpmf);
} 




/////////////
////////////

// // [[Rcpp::export]]
// double log_sp_prob_permute(arma::colvec partition,
//                            arma::colvec baseline, 
//                            double crp_para, 
//                            arma::colvec shrink_para,
//                            arma::colvec permutation){
//   
//   int N = partition.size();
//   arma::colvec logpmf(N-1);
//   
//   // Convert the permutation to 0-based indexing
//   arma::uvec zero_based_permu = arma::conv_to<arma::uvec>::from(permutation - 1);  //very important
//   // Apply the permutation
//   partition = partition(zero_based_permu);  //very important
//   baseline = baseline(zero_based_permu);  //very important
//   
//   // Use arma::uvec for storing indices
//   arma::uvec pre_labels_idx, pre_bases_idx, unique_labels_idx;
//   arma::colvec unique_labels, a, b, d;
//   int label, base, Nclust, denom;
//   
//   // Efficiently create a vector filled with crp_para
//   arma::colvec crp_vec(N);
//   crp_vec.fill(crp_para);
//   
//   a.set_size(N);  // max size allocation outside loop
//   b.set_size(N);  // max size allocation outside loop
//   
//   for (int i = 1; i < N; ++i){
//     label = partition(i);
//     arma::colvec current_partition_head = partition.head(i);
//     base = baseline(i);
//     arma::colvec current_baseline_head = baseline.head(i);
//     
//     unique_labels = unique(current_partition_head);
//     Nclust = unique_labels.size();
//     pre_bases_idx = find(current_baseline_head == base);
//     denom = pre_bases_idx.size();
//     
//     a.head(Nclust + 1).fill(0);  // reset values
//     b.head(Nclust + 1).fill(0);  // reset values
//     
//     for (size_t j = 0; j < Nclust; ++j){
//       a(j) = accu(current_partition_head == unique_labels(j));
//     }      
//     
//     a(Nclust) = crp_para;
//     a.head(Nclust + 1) /= (crp_para + i);
//      
//     if (denom == 0){
//       b(Nclust) = 1;
//     } else{      
//       for (size_t c = 0; c < Nclust; ++c){
//         b(c) = accu((current_partition_head == unique_labels(c)) % (current_baseline_head == base)) / static_cast<double>(denom);
//       } 
//     }    
//     
//     b.head(Nclust + 1) = exp(shrink_para(i) * b.head(Nclust + 1));
//      
//     d = a.head(Nclust + 1) % b.head(Nclust + 1);
//     d /= accu(d);
//      
//     unique_labels_idx = find(unique_labels == label);
//     if (unique_labels_idx.is_empty()){
//       logpmf(i - 1) = log(d(Nclust));
//     } else {      
//       logpmf(i - 1) = log(d(unique_labels_idx(0)));
//     }      
//     
//   }      
//   
//   return accu(logpmf);
// }    






/////
///////

// [[Rcpp::export]]
arma::colvec prior_simulate_sp_partition(arma::colvec baseline, 
                                         arma::colvec shrink_para, double crp_para){
  arma::colvec pre_labels, pre_bases, unique_labels, pp;
  int base, Nclust, denom, newz;
  
  int N = baseline.size();
  arma::colvec z(N);
  z[0] = 1;
  
  for (int i=1; i<N; i++){
    pre_labels = z(arma::span(0,i-1));
    base = baseline(i);
    pre_bases = baseline(arma::span(0,i-1));
    
    unique_labels = unique(pre_labels);
    Nclust = unique_labels.size();
    denom = accu(pre_bases == base);
    
    arma::colvec b(Nclust+1);
    arma::colvec a(Nclust+1);
    arma::colvec possible_label = arma::linspace<arma::vec>(1, Nclust+1, Nclust+1);
    
    for (int j=0; j<Nclust; j++){
      a(j) = accu(pre_labels == unique_labels(j));
    }
    a(Nclust) = crp_para;
    a = a / (crp_para + i);
    
    if (denom == 0){
      b(Nclust) = 1;
    } else{
      for (int c=0; c<Nclust; c++){
        b(c) = accu((pre_labels == unique_labels(c)) % (pre_bases == base)) / denom;
      }
    }
    b = exp(shrink_para(i) * b);
    
    pp = a % b;
    newz = Rcpp::RcppArmadillo::sample(possible_label, 1, 1, pp)[0];
    z(i) = newz;
  }
  return z;
}











