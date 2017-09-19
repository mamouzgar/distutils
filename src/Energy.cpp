// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
List ComputeElasticEnergy(NumericMatrix X, NumericMatrix NodePositions, NumericMatrix ElasticMatrix, NumericVector Dists, double BranchingFee) {
  
  int k = NodePositions.nrow(),
    // n = NodePositions.ncol(),
    i,
    j;
    
  double MSE = sum(Dists)/X.nrow(),
    EP = 0,
    RP = 0,
    TotEnergy = 0;
  
  arma::mat Lambda = arma::mat(ElasticMatrix.begin(), k, k, false);
  // arma::mat NP = arma::mat(NodePositions.begin(), k, n, false);
  arma::vec Mu = Lambda.diag(0);
  Lambda.diag(0).zeros();
  
  arma::uvec StarCenterIndices = arma::find(Mu);
  arma::uvec PosLambdaIdxs = find(Lambda);
  arma::umat PosLambdaIdxsMat = ind2sub( size(Mu), PosLambdaIdxs);
  
  arma::vec dev;
  double l;
  
  for(i=0; i<k; i++){
    dev = NodePositions.row(PosLambdaIdxsMat(i,1)) - NodePositions.row(PosLambdaIdxsMat(i,2));
    l = Lambda(PosLambdaIdxsMat(i,1),PosLambdaIdxsMat(i,2));
    EP =+ l*dot(dev,dev);
  }
  
  arma::uvec leafs;
  arma::rowvec tVect;
  
  for(i=0; i<StarCenterIndices.size(); i++){
    leafs = find(Lambda.col(StarCenterIndices[i]));
    j = leafs.size();
    tVect = NodePositions.row(StarCenterIndices[i]);
    for(k=0; k<j; k++){
      tVect =- (1/k+1)*NodePositions.row(leafs[k]);
    }
    RP =+ Mu[StarCenterIndices[i]]*dot(dev,dev);
  }
  
  TotEnergy = EP + RP + MSE;
  
  List RetList = List::create(Named("ElasticEnergy") = TotEnergy,
                              Named("EP") = EP,
                              Named("RP") = RP,
                              Named("MSE") = MSE);

  return RetList;
  
}

