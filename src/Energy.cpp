// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Compute the elastic energy associated with a particular configuration 
//' 
//' This function computes the elastic energy associate to a set of points and graph embedded
//' into them. See XXX for reference
//' 
//' @param X A numeric n-by-m matrix containing the position of n data points m-dimensional points
//' @param NodePositions A numeric k-by-m matrix containing the position of the k nodes of the embedded graph
//' @param ElasticMatrix A numeric l-by-l matrix containing the elastic parameters associates with the edge
//' of the embedded graph
//' @param Dists A numeric vector containind the squared distance of the data points to the closest node of the graph
//' @param BranchingFee a numeric value currently unused
//' 
//' @return A list with four elements:
//' * ElasticEnergy is the total energy
//' * EP is the EP component of the energy
//' * RS is the RS component of the energy
//' * MSE is the MSE component of the energy
//' @md
//' 
//' @export
//' 
//' @examples 
//' 
// [[Rcpp::export]]
List ElasticEnergy(NumericMatrix X, NumericMatrix NodePositions, NumericMatrix ElasticMatrix, NumericVector Dists, double BranchingFee) {
  
  int k = NodePositions.nrow(),
    n = NodePositions.ncol(),
    i,
    j;
    
  double MSE = sum(Dists)/X.nrow(),
    EP = 0,
    RP = 0,
    TotEnergy = 0;
  
  arma::mat Lambda = arma::mat(ElasticMatrix.begin(), k, k, true);
  arma::mat NP = arma::mat(NodePositions.begin(), k, n, false);
  arma::vec Mu = Lambda.diag(0);
  Lambda.diag(0).zeros();
  
  arma::uvec StarCenterIndices = arma::find(Mu);
  arma::uvec PosLambdaIdxs = find(Lambda);
  arma::umat PosLambdaIdxsMat = ind2sub( size(Lambda), PosLambdaIdxs);
  
  arma::rowvec dev;
  double l;
  
  for(i=0; i<PosLambdaIdxs.size(); i++){
    dev = NP.row(PosLambdaIdxsMat(0,i)) - NP.row(PosLambdaIdxsMat(1,i));
    l = Lambda(PosLambdaIdxsMat(0,i),PosLambdaIdxsMat(1,i));
    EP += l*arma::dot(dev,dev);
  }
  
  arma::uvec leafs;
  
  for(i=0; i<StarCenterIndices.size(); i++){
    leafs = find(Lambda.col(StarCenterIndices[i]));
    j = leafs.size();
    dev = NP.row(StarCenterIndices[i]);
    for(k=0; k<j; k++){
      dev -= 1/((double)j)*NP.row(leafs[k]);
    }
    RP += Mu[StarCenterIndices[i]]*arma::dot(dev,dev);
  }
  
  
  TotEnergy = EP + RP + MSE;
  
  List RetList = List::create(Named("ElasticEnergy") = TotEnergy,
                              Named("EP") = EP,
                              Named("RP") = RP,
                              Named("MSE") = MSE);

  return RetList;
  
}

