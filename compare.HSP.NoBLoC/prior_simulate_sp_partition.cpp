#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;


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







