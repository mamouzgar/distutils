// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Compute the number of points within a given distance
//' 
//' This function computes the distance between a matrix of data points and a
//' single reference points. A sorted distance vector define the reference distances
//' to consider
//' 
//' @param Ar A numeric n-by-m matrix containing the position of n data points m-dimensional points
//' @param P A m dimensional vector matrix containing the position of the reference point
//' @param DVect A k dimensional vector matrix containing the k distances to consider
//' 
//' @return A k-dimensional vector containing the number of points with a distance lower than k
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
List RadialCount(NumericMatrix Ar, NumericVector Pr, NumericVector SquaredAr, NumericVector DVect) {
  
  // Get dimension of the matrix
  int n = Ar.nrow(), 
    k = Ar.ncol();
  
  // Copy matrices to internal structures
  // std::cout << "Init A" << std::endl;
  arma::mat A = arma::mat(Ar.begin(), n, k, false); 
  
  // std::cout << "Init B" << std::endl;
  arma::colvec B = arma::vec(Pr.begin(), Pr.size(), false); 
  
  // std::cout << "Init An";
  arma::colvec An = arma::vec(SquaredAr.begin(), SquaredAr.size(), false);
  
  // std::cout << "Init Bn";
  arma::colvec Bn = sum(square(B),0);
  
  arma::mat C = -2 * (A * B);
  C.each_col() += An;
  C.each_row() += Bn.t();
  
  
  arma::colvec Dn = arma::vec(DVect.begin(), DVect.size(), false);
  arma::uvec PCount = arma::uvec(Dn.size());
  PCount.zeros();
  
  arma::colvec CVect_Org = sqrt(C.col(0));
  CVect_Org.replace(arma::datum::nan, 0);
  
  arma::uvec ValidPoints = arma::uvec(CVect_Org.size());
  
  double MaxDist = CVect_Org.max();
  double MinDist = CVect_Org.min();
  double MaxVal = Dn.max();
  arma::uvec ToUse = arma::uvec(2);
  
  // std::cout << MinDist;
  
  // std::cout << CVect_Org;
  
  arma::colvec CVect;
  
  if(MaxDist > MaxVal){
    ValidPoints = arma::find(CVect_Org < MaxVal);
    // std::cout << ValidPoints;
    CVect = CVect_Org.elem(ValidPoints);
    // std::cout << CVect;
    MaxDist = CVect.max();
  } else {
    CVect = CVect_Org;
    for(int j=0; j<ValidPoints.size(); j++){
      ValidPoints[j] = j;
    }
  }
  
  CVect = arma::sort(CVect);
  
  int p_idx = 0;
  
  ToUse[0] = 0;
  ToUse[1] = PCount.size()-1;
  
  // std::cout << Dn;
  
  // std::cout << MaxDist;
  
  while(Dn[ ToUse[1] ] > MaxDist){
    // std::cout << i << Dn[i] << " " << MaxDist << std::endl;
    ToUse[1]--;
  }
  if(ToUse[1] < PCount.size()-1){
    ToUse[1]++;
  }
  
  // std::cout << MinDist;
  
  while(Dn[ ToUse[0] ] < MinDist){
    ToUse[0]++;
  }
  if(ToUse[0]>0){
    ToUse[0]--;
  }
    
  // std::cout << ToUse;
  
  // std::cout << j << " " << i << " " << p_idx << std::endl;
  
  for(int k = ToUse[0]; k < ToUse[1]; k++){
    // std::cout << CVect[p_idx] << " " << Dn[k] << std::endl;
    while( CVect[p_idx] <= Dn[k] ){
    // while( (CVect[p_idx] <= Dn[k]) & (p_idx < CVect.size()) ){
      p_idx++;
      // std::cout << CVect[p_idx] << " " << Dn[k] << std::endl;
    }
    PCount[k] = p_idx;
  }
  
  p_idx = CVect.size();
  
  for(int k = ToUse[1]; k < PCount.size(); k++){
    PCount[k] = p_idx;
  }
  
  ValidPoints = ValidPoints + 1;
  
  std::vector<uint> PCountVect = arma::conv_to<std::vector<uint> >::from(PCount);
  std::vector<uint> VPointsVect = arma::conv_to<std::vector<uint> >::from(ValidPoints);
  
  List RetList = List::create(Named("PCount") = wrap(PCountVect),
                              Named("Idxs") = wrap(VPointsVect));
  
  //List RetList = List::create(Named("PCount") = wrap(PCountVect));
  
  return RetList;
}