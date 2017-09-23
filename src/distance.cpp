// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;

//' Compute the partial distance between matrices
//' 
//' This function computes the distance between a matrix of data points and a
//' matrix of reference points. Both matrices need to have the same number of
//' dimensions (columns). The code was taked from here: // C++ code taken from here:
//' https://www.r-bloggers.com/pairwise-distances-in-r/
//' 
//' @param Ar A numeric n-by-m matrix containing the position of n data points m-dimensional points
//' @param Br A numeric k-by-m matrix containing the position of k reference m-dimensional points
//' 
//' @return A numeric n-by-k matrix reporting the euclidean distance between the n data points and 
//' the k reference points
//' 
//' @export
//' 
//' @examples 
//' 
//' A <- matrix(runif(10000*100), nrow = 10000)
//' B <- matrix(runif(100*100), nrow = 100)
//' 
//' library(distutils)
//' 
//' print(system.time(C1 <- PartialDistance(A, B)))
//' print(system.time(C2 <- as.matrix(dist(rbind(A,B)))[1:10000, 10001:10100]))
//' 
//' summary(as.vector(C1 - C2))
//' 
//' 
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




//' Parition a set of data ponints into groups based on the distance from a set of reference points
//' 
//' This function computes the distance between a matrix of data points and a
//' matrix of reference points. Then associate each point with a id inicating the closet
//' reference point
//' 
//' @param Ar A numeric n-by-m matrix containing the position of n data points m-dimensional points
//' @param Br A numeric k-by-m matrix containing the position of k reference m-dimensional points
//' @param SquaredAr A numeric vector reporting the sum, by row, of Ar (i.e. SquaredX = rowSums(Ar^2)).
//' This is done to speed up the calculation when the same set o data points is clustered against
//' a different set of reference point (e.g. in k-means)
//' 
//' @return A list with two elements:
//' * Patition is an integer vector indicating, for each point, the index of the closest reference point
//' * Distance is a numeric vector indicating, for each point, the squared distance to the closest reference point
//' @md
//' 
//' @export
//' 
//' @examples 
//' 
//' A <- matrix(runif(10000*100), nrow = 10000)
//' B <- matrix(runif(100*100), nrow = 100)
//'   
//'   library(distutils)
//'   
//'   print(
//'     system.time(
//'       Part <- Partition(Ar = A, Br = B, SquaredAr = rowSums(A^2))
//'     )
//'   )
//'   
//'   print(
//'     system.time({
//'       C2 <- as.matrix(dist(rbind(A,B)))[1:10000, 10001:10100]
//'       Partition <- apply(C2, 1, which.min)
//'       Dist <- apply(C2, 1, min)^2
//'     })
//'   )
//'   
//'   summary(Partition - Part$Patition)
//'   summary(Dist - Part$Dist)
//' 
// [[Rcpp::export]]
List Partition(NumericMatrix Ar, NumericMatrix Br, NumericVector SquaredAr) {
  
  // Get dimension of the matrix
  int n = Ar.nrow(), 
    m = Br.nrow(),
    k = Br.ncol();
  
  // Copy matrices to internal structures
  arma::mat A = arma::mat(Ar.begin(), n, k, false); 
  arma::mat B = arma::mat(Br.begin(), m, k, false); 
  
  arma::colvec An = arma::vec(SquaredAr.begin(), SquaredAr.size(), false);
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
  
  std::vector<uint> Partition = arma::conv_to<std::vector<uint> >::from(IdxVect+1);
  std::vector<double> Distance = arma::conv_to<std::vector<double> >::from(dist);
  
  List RetList = List::create(Named("Patition") = wrap(Partition),
                              Named("Dist") = wrap(Distance));
  
  return RetList;
  
}

