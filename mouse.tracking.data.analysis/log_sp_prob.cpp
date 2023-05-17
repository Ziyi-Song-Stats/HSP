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








