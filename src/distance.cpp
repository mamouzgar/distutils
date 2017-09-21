// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// C++ code taken from here : https://www.r-bloggers.com/pairwise-distances-in-r/
//

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix PartialDistance(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(), 
    n = Br.nrow(),
    k = Ar.ncol();
  
  arma::mat A = arma::mat(Ar.begin(), m, k, false); 
  arma::mat B = arma::mat(Br.begin(), n, k, false); 
  
  arma::colvec An =  sum(square(A),1);
  arma::colvec Bn =  sum(square(B),1);
  
  arma::mat C = -2 * (A * B.t());
  C.each_col() += An;
  C.each_row() += Bn.t();
  
  return wrap(sqrt(C)); 
}


// [[Rcpp::export]]
List Partition(NumericMatrix X, NumericMatrix NP, NumericVector SquaredX) {
  
  // Get dimension of the matrix
  int n = X.nrow(), 
    m = NP.nrow(),
    k = NP.ncol();
  
  // Copy matrices to internal structures
  arma::mat A = arma::mat(X.begin(), n, k, false); 
  arma::mat B = arma::mat(NP.begin(), m, k, false); 
  
  arma::colvec An = SquaredX;
  arma::colvec Bn = sum(square(B),1);
  
  arma::mat C = -2 * (A * B.t());
  C.each_col() += An;
  C.each_row() += Bn.t();
  
  arma::uvec IdxVect = arma::index_min(C, 1);
  arma::urowvec IdxVect2(n);
    
  for(int i = 0; i< n; i++){
    IdxVect2[i] = i;
  }
  
  arma::umat locs(2, n);
  locs.row(0) = IdxVect2;
  locs.row(1) = IdxVect.t();
  
  arma::uvec eids = sub2ind(size(C), locs);
  arma::vec dist = C.elem(eids);
  
  std::vector<int> Partition = arma::conv_to<std::vector<int> >::from(IdxVect+1);
  std::vector<double>Disttance = arma::conv_to<std::vector<double> >::from(dist);
  
  List RetList = List::create(Named("Patition") = wrap(Partition),
                              Named("Dist") = wrap(Disttance));
  
  return RetList;
  
}

