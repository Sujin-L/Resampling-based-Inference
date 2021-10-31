#### exhaustive leave-2-out cross-validation GLM for our data ####
source("range_lambda.R")

cv.glm <- function(K, x, y){
  
  range_lambda <- range_lambda(x, y)
  n <- nrow(x)
  weights_cv <- c(rep(1, n*(K-1)/K), 1000)
  
  collect <- NULL
  
  s <- combn(n, n/K)
  for(i in 1:ncol(s)){
    
    ss <- s[,i]
    
    ## data transformation
    x_cv <- x[-ss,]; y_cv <- y[-ss]
    
    for(j in 1:ncol(x_cv)){
      x_cv[,j] <- x_cv[,j] - apply(x_cv,2,mean)[j]
    }
    y_cv <- scale(y[-ss], center = TRUE, scale = FALSE)
    
    
    ## bind the data and glmnet
    train <- list(X = as.matrix(rbind(x_cv, rep(1,ncol(x_cv)))), y=c(y_cv,0)) 
    test <- list(X = as.matrix(x[ss,]), y = y[ss])
    
    seq_lambda <- seq(log(range_lambda$min), log(range_lambda$max), length=100)
    
    gfit <- glmnet(x = train$X, y = train$y, weights = weights_cv, intercept = FALSE,
                   lambda = exp(seq_lambda))
    
    y_hat <- predict(gfit, newx = test$X, s = gfit$lambda, type="response")
    mse <- apply(y_hat, 2, function(x) mean((x - test$y)^2))
    
    collect <- rbind(collect, cbind(gfit$lambda, mse))
  }
  
  colnames(collect) <- c("lambda","mse")
  result <- as.data.frame(collect)
  return(result)
}



