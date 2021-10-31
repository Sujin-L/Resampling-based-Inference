library(readxl)
library(tidyverse)

#### transform X ####
min_pos <- function(x) min(x[x > 0])
x <- X_rumen_o
trans_X <- function(x){
  
  X <- data.matrix(x)
  n <- nrow(X); p <- ncol(X)
  
  # calculate row sum
  X_sum <- apply(X, 1, sum)
  
  # scaling
  X1 <- matrix(NA, n, p)
  for(i in 1:nrow(X)){
    X1[i,] <- X[i,]/X_sum[i]
  }
  colnames(X1) <- colnames(X)
  
  # substitute 0 value with minimum value
  X1[X1==0] <- min_pos(X1)
  
  # log-transformation
  log_X1 <- log(X1)
  
  # centering log-transformed X
  for(j in 1:p){
    log_X1[,j] <- log_X1[,j] - apply(log_X1,2,mean)[j]
  }
  
  return(log_X1)
}


#### transform Y ####
trans_Y <- function(y){
  
  #centering Y
  y1 <- scale(y, center = TRUE, scale = FALSE)
  y1 <- as.vector(y1)
  
  return(y1)
}
