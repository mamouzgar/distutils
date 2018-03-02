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
List ElasticEnergy(NumericMatrix X,
                   NumericMatrix NodePositions,
                   NumericMatrix ElasticMatrix,
                   NumericVector Dists) {
  
  int k = NodePositions.nrow(),
    n = NodePositions.ncol(),
    i,
    j;
  
  double MSE = sum(Dists)/X.nrow(),
    EP = 0,
    RP = 0,
    TotEnergy = 0;
  
  arma::mat dev;
  
  arma::mat EM = arma::mat(ElasticMatrix.begin(), k, k, true);
  arma::mat NP = arma::mat(NodePositions.begin(), k, n, false);
  arma::vec Mu = EM.diag(0);
  
  arma::mat Lambda = trimatu(EM,  1);
  arma::uvec StarCenterIndices = arma::find(Mu > 0);
  
  EM.diag(0).zeros(); 
  
  arma::uvec PosLambdaIdxs = find(Lambda);
  arma::umat PosLambdaIdxsMat = ind2sub(size(Lambda), PosLambdaIdxs);
  
  dev = NP.rows(PosLambdaIdxsMat.row(0)) - NP.rows(PosLambdaIdxsMat.row(1));
  
  arma::mat l = Lambda(find(Lambda > 0));
  
  // Rcpp::Rcout << "lpenalized" << std::endl;
  // Rcpp::Rcout << lpenalized << std::endl;
  
  arma::mat tEP = l.t() * sum(pow(dev, 2), 1);
  
  EP = tEP(0,0);
  
  // Rcpp::Rcout << EP << std::endl;
  
  arma::uvec leafs;
  // double K = 0;
  
  for(i=0; i<StarCenterIndices.size(); i++){
    
    leafs = find(EM.col(StarCenterIndices(i)));

    j = leafs.size();
    
    dev = NP.row(StarCenterIndices(i)) - sum(NP.rows(leafs)) / j;

    RP += Mu(StarCenterIndices(i))*arma::dot(dev,dev);

  }
  
  TotEnergy = EP + RP + MSE;
  
  List RetList = List::create(Named("ElasticEnergy") = TotEnergy,
                              Named("EP") = EP,
                              Named("RP") = RP,
                              Named("MSE") = MSE);
  
  
  return RetList;
  
}

















//' Compute the elastic energy associated with a particular configuration (Old version)
//' 
//' This function computes the elastic energy associate to a set of points and graph embedded
//' into them. See XXX for reference
//' 
//' @param X A numeric n-by-m matrix containing the position of n data points m-dimensional points
//' @param NodePositions A numeric k-by-m matrix containing the position of the k nodes of the embedded graph
//' @param ElasticMatrix A numeric l-by-l matrix containing the elastic parameters associates with the edge
//' of the embedded graph
//' @param Dists A numeric vector containind the squared distance of the data points to the closest node of the graph
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
List ElasticEnergy_V0(NumericMatrix X,
                   NumericMatrix NodePositions,
                   NumericMatrix ElasticMatrix,
                   NumericVector Dists) {
  
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
  arma::umat PosLambdaIdxsMat = ind2sub(size(Lambda), PosLambdaIdxs);
  
  arma::rowvec dev;
  double l;
  
  for(i=0; i<PosLambdaIdxs.size(); i++){
    dev = NP.row(PosLambdaIdxsMat(0,i)) - NP.row(PosLambdaIdxsMat(1,i));
    l = Lambda(PosLambdaIdxsMat(0,i),PosLambdaIdxsMat(1,i));
    EP += l*arma::dot(dev,dev);
  }
  
  EP /= 2;
  
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




























//' Compute the penalized elastic energy associated with a particular configuration 
//' 
//' This function computes the elastic energy associate to a set of points and graph embedded
//' into them. See XXX for reference
//' 
//' @param X A numeric n-by-m matrix containing the position of n data points m-dimensional points
//' @param NodePositions A numeric k-by-m matrix containing the position of the k nodes of the embedded graph
//' @param ElasticMatrix A numeric l-by-l matrix containing the elastic parameters associates with the edge
//' of the embedded graph
//' @param Dists A numeric vector containind the squared distance of the data points to the closest node of the graph
//' @param alpha 
//' @param beta
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
List PenalizedElasticEnergy(NumericMatrix X,
                   NumericMatrix NodePositions,
                   NumericMatrix ElasticMatrix,
                   NumericVector Dists,
                   double alpha,
                   double beta) {
  
  int k = NodePositions.nrow(),
    n = NodePositions.ncol(),
    i,
    j;
  
  double MSE = sum(Dists)/X.nrow(),
    EP = 0,
    RP = 0,
    TotEnergy = 0;
  
  arma::mat dev;
  
  arma::mat EM = arma::mat(ElasticMatrix.begin(), k, k, true);
  arma::mat NP = arma::mat(NodePositions.begin(), k, n, false);
  arma::vec Mu = EM.diag(0);
  
  arma::mat Lambda = trimatu(EM,  1);
  arma::uvec StarCenterIndices = arma::find(Mu > 0);
  
  EM.diag(0).zeros(); 
  
  // Rcpp::Rcout << "Lambda" << std::endl;
  // Rcpp::Rcout << Lambda << std::endl;
  
  arma::uvec PosLambdaIdxs = find(Lambda);
  arma::umat PosLambdaIdxsMat = ind2sub(size(Lambda), PosLambdaIdxs);
  
  // Rcpp::Rcout << "PosLambdaIdxsMat" << std::endl;
  // Rcpp::Rcout << PosLambdaIdxsMat << std::endl;  
  
  dev = NP.rows(PosLambdaIdxsMat.row(0)) - NP.rows(PosLambdaIdxsMat.row(1));
  
  // Rcpp::Rcout << "dev" << std::endl;
  // Rcpp::Rcout << dev << std::endl;
  
  arma::mat l = Lambda(find(Lambda > 0));
  
  // Rcpp::Rcout << "l" << std::endl;
  // Rcpp::Rcout << l << std::endl;
  
  arma::mat BinEM = EM;
  BinEM(find(BinEM >0)).ones();
  arma::rowvec Ks = arma::sum(BinEM);
  
  // Rcpp::Rcout << "Ks" << std::endl;
  // Rcpp::Rcout << Ks << std::endl;
  
  // Rcpp::Rcout << Ks(PosLambdaIdxsMat.row(0)) << std::endl;
  
  // Rcpp::Rcout << Ks(PosLambdaIdxsMat.row(1)) << std::endl;
  
  // Rcpp::Rcout << arma::join_cols(Ks(PosLambdaIdxsMat.row(0)).t(), Ks(PosLambdaIdxsMat.row(1)).t()) << std::endl;
  
  arma::mat lp = arma::max(arma::join_cols(Ks(PosLambdaIdxsMat.row(0)).t(), Ks(PosLambdaIdxsMat.row(1)).t()), 0);
  
  // Rcpp::Rcout << lp.t() << std::endl;
  
  lp = lp-2;
  
  // Rcpp::Rcout << lp.t() << std::endl;
  
  lp(find(lp<0)).zeros();
  
  // Rcpp::Rcout << lp.t() << std::endl;
  
  // Rcpp::Rcout << "lp" << std::endl;
  // Rcpp::Rcout << lp << std::endl;
  
  arma::mat lpenalized = l + alpha*lp.t();
  
  // Rcpp::Rcout << "lpenalized" << std::endl;
  // Rcpp::Rcout << lpenalized << std::endl;
  
  arma::mat tEP = lpenalized.t() * sum(pow(dev, 2), 1);
  
  EP = tEP(0,0);
  
  // Rcpp::Rcout << EP << std::endl;
  
  arma::uvec leafs;
  // double K = 0;
  
  for(i=0; i<StarCenterIndices.size(); i++){
    
    leafs = find(EM.col(StarCenterIndices(i)));
    
    // Rcpp::Rcout << leafs << std::endl;
    
    j = leafs.size();
    
    // Rcpp::Rcout << j << std::endl;
    
    // Rcpp::Rcout << NP.row(StarCenterIndices(i)) << std::endl;
    
    // Rcpp::Rcout << sum(NP.rows(leafs)) << std::endl;
    
    // Rcpp::Rcout << leafs << std::endl;
    
    dev = NP.row(StarCenterIndices(i)) - sum(NP.rows(leafs)) / j;
    
    // for(k=0; k<j; k++){
    //  dev -= 1/((double)j)*NP.row(leaves[k]);
    //}
    
    // Rcpp::Rcout << dev << std::endl;
    
    RP += Mu(StarCenterIndices(i))*pow(j, beta)*arma::dot(dev,dev);
    
    // Rcpp::Rcout << dev << std::endl;
      
    // Rcpp::Rcout << Mu(StarCenterIndices(i)) << std::endl;
    
    // Rcpp::Rcout << arma::dot(dev,dev) << std::endl;
    
    // Rcpp::Rcout << pow(j, beta) << std::endl;

    
  }
  
  TotEnergy = EP + RP + MSE;
  
  List RetList = List::create(Named("ElasticEnergy") = TotEnergy,
                              Named("EP") = EP,
                              Named("RP") = RP,
                              Named("MSE") = MSE);
  
  
  return RetList;
  
}






























// [[Rcpp::export]]
List ComputeWeightedAverage(NumericMatrix X,
                            IntegerVector partition,
                            NumericVector PointWeights,
                            unsigned int NumberOfNodes) {
  
  unsigned int n = X.nrow(),
    m = X.ncol(),
    i, j;
  
  arma::mat Data = arma::mat(X.begin(), n, m, false);
  // Rcout << "Data created" << std::endl;
  
  arma::irowvec Part = arma::irowvec(partition.begin(), n, false);
  // Rcout << "Part created" << std::endl;
  
  arma::rowvec PWei = arma::rowvec(PointWeights.begin(), n, false);
  // Rcout << "PWei created" << std::endl;
  
  arma::mat NodeClusterCenters = arma::mat(NumberOfNodes, m, arma::fill::zeros);
  // Rcout << "NodeClusterCenters created" << std::endl;
  
  arma::vec NodeClusterRelativeSize = arma::vec(NumberOfNodes, arma::fill::zeros);
  // Rcout << "NodeClusterRelativeSize created" << std::endl;
  
  double TotalWeight = sum(PointWeights);
  arma::uvec Indices;
  
  for(i = 1; i <= NumberOfNodes; i++){
    
    Indices = find(Part == i);
    
    // Rcout << "i= " << i << "size= " << Indices.size() << std::endl;
    
    NodeClusterRelativeSize(i-1) = sum(PWei(Indices));
    
    for(j = 0; j < Indices.size(); j++){
      NodeClusterCenters.row(i-1) += Data.row(Indices(j)) * PWei(Indices(j))/NodeClusterRelativeSize(i-1);
    }
    
  }
  
  std::vector<double> RelCluSize = arma::conv_to<std::vector<double> >::from(NodeClusterRelativeSize/TotalWeight);
  
  List RetList = List::create(Named("NodeClusterCenters") = wrap(NodeClusterCenters),
                              Named("NodeClusterRelativeSize") = wrap(RelCluSize)
                              );
  
  return RetList;

}











//' 
// [[Rcpp::export]]
arma::mat FitGraph2DataGivenPartition(
    NumericMatrix X,
    NumericVector PointWeights,
    NumericMatrix NodePositions,
    NumericMatrix SpringLaplacianMatrix,
    IntegerVector partition,
    bool FastSolve) {
  
  unsigned int n = X.nrow(),
    m = X.ncol(),
    NumberOfNodes = SpringLaplacianMatrix.nrow(),
    i, j;
  
  arma::mat Data = arma::mat(X.begin(), n, m, false);
  arma::irowvec Part = arma::irowvec(partition.begin(), n, false);
  arma::rowvec PWei = arma::rowvec(PointWeights.begin(), n, false);
  arma::mat LMat = arma::mat(SpringLaplacianMatrix.begin(), NumberOfNodes, NumberOfNodes, true);
  
  
  arma::mat NodeClusterCenters = arma::mat(NumberOfNodes, m, arma::fill::zeros);
  arma::vec NodeClusterRelativeSize = arma::vec(NumberOfNodes, arma::fill::zeros);
  
  double TotalWeight = sum(PointWeights);
  arma::uvec Indices;
  
  Indices = find(Part == 0);
  
  // Rcout << "Ignored points=" << Indices.size() << std::endl;
  
  for(i = 1; i <= NumberOfNodes; i++){
    
    Indices = find(Part == i);
    
    // Rcout << "i= " << i << "size= " << Indices.size() << std::endl;
    
    NodeClusterRelativeSize(i-1) = sum(PWei(Indices));
    
    // if(Indices.size() == 0){
    //   Rcout << "Warning: point " << i << " is not associated with any data point" << std::endl;
    // }
    
    for(j = 0; j < Indices.size(); j++){
      NodeClusterCenters.row(i-1) += Data.row(Indices(j)) * PWei(Indices(j))/NodeClusterRelativeSize(i-1);
    }
    
  }
  
  // std::vector<double> RelCluSize = arma::conv_to<std::vector<double> >::from(NodeClusterRelativeSize/TotalWeight);
  NodeClusterRelativeSize = NodeClusterRelativeSize/TotalWeight;
  
  LMat.diag(0) += NodeClusterRelativeSize;
  
  NodeClusterCenters.each_col() %= NodeClusterRelativeSize;
  
  arma::mat FinalVal;
  
  if(FastSolve){
    FinalVal = arma::solve(LMat, NodeClusterCenters, arma::solve_opts::fast);
  } else {
    FinalVal = arma::solve(LMat, NodeClusterCenters);
  }
  
  return FinalVal;
}










