X = read.table('../rpgraph2/examples/tree23/tree23.data');

TT <- NULL

for(inflationFactor in round(10^seq(0, 3.4, by=.1))){
  
  print(inflationFactor)
  
  X1 <- X
  
  # inflate the number of points
  for(i in 1:inflationFactor){
    X1 <- rbind(X1, X)
  }
  
  Npoints = nrow(X1)
  
  Edges = matrix(c(1,2,2,3,3,4,4,5,5,6,6,7,7,8,4,9,9,10,10,11,11,12), byrow = TRUE, ncol = 2)
  
  NNodes = max(Edges)
  NodePositions = X1[sample(x = 1:nrow(X1), size = NNodes, replace = TRUE),]
  
  #ElasticMatrix = MakeUniformElasticMatrix(Edges,0.01,1);
  #ElasticMatrix = MakeUniformElasticMatrix(Edges,0,0);
  
  X1 <- data.matrix(X1)
  NodePositions <- data.matrix(NodePositions)
  SquaredX = rowSums(X1^2)
  
  TT <- rbind(TT,
              c(
                system.time(
                  Dists <- RFastDistance::fastPdist2List(X1, NodePositions, SquaredX = SquaredX, BlockSize = 1)
              ), nrow(X1)))
  
}

plot(TT[,6], TT[,3])

