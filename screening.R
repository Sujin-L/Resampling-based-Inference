library(MASS)

install.packages("usethis")
library(usethis)
use_git_config(user.name = "Sujin-L", user.email = "sujennt77@gmail.com")
git_vaccinate()

create_github_token()


stat <- function(X,y){
  
  n <- nrow(X); p <- ncol(X)
  D <- t(X) %*% X
  inv <- ginv(D)
  Q <- rep(1, ncol(X))
  A <- t(Q)%*%inv%*%Q
  
  #### estimate in original data ####
  beta_ols <- inv%*%t(X)%*%y
  M <- diag(ncol(X))-inv%*%Q%*%solve(A)%*%t(Q)
  
  ## estimate beta
  beta_hat <- M%*%beta_ols; beta_hat <- as.vector(t(beta_hat))
  names(beta_hat) <- colnames(X) 
  
  ## t-statistic
  SSE_o <- t(y)%*%(diag(nrow(X))-X%*%inv%*%t(X))%*%y + t(beta_ols)%*%Q%*%solve(A)%*%t(Q)%*%beta_ols
  var_beta_o <- as.vector(SSE_o)*( M %*% inv) /(n-p+1)
  t_ref <- beta_hat/sqrt(diag(var_beta_o))
  
  value <- list(beta = beta_hat, t= t_ref, M = M, inv = inv, Q=Q, A=A, se=sqrt(diag(var_beta_o)))
  return(value)
}

