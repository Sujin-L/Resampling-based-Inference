
#### Caculate value of minimum and maximum lambda ####
range_lambda <- function(x, y){
  
  min_lambda <- max_lambda <- NULL
  
  n <- nrow(x)
  weights <- c(rep(1,n), 1000)
  
  gfit <- glmnet(x = rbind(x, rep(1,ncol(x))), y = c(y,0), weights = weights, intercept = FALSE)
  
  min_lambda <- min(gfit$lambda) 
  max_lambda <- max(gfit$lambda)
  
  if(max_lambda<1){
    max1 <- 1
  }
  else{
    max1 <- max_lambda
  }
  
  result = list(min = max(min_lambda), max = max(max_lambda), max1 = max1)
  return(result)
}  

