
#### re-estimating the regression parameters ####
re_estimation <- function(sig_var, ori_X, ori_y){
  
  sig_index <- which(colnames(ori_X)%in%rownames(sig_var))
  s <- length(sig_index)
  
  ori_X <- ori_X/100 #scaling for real data
  ori_X[ori_X==0] <- min_pos(ori_X)
  
  if(s==1){
    
    new_X <- cbind(ori_X[,sig_index], apply(ori_X[,-sig_index],1,sum))
    trans_X <- log(new_X[,1]/new_X[,2])
    
    lfit <- lm(ori_y ~ trans_X)
    beta_hat <- lfit$coefficients
    names(beta_hat) <- c("intercept", rownames(sig_var))
    return(beta_hat)
  }
  
  if(s>1){
    
    log_X <- log(ori_X)
    log_X <- new_X <- log_X[, sig_index] 
    
    for(j in 1:ncol(new_X)){
      new_X[,j] <- new_X[,j] - apply(new_X,2,mean)[j]
    }
    new_X <- as.matrix(new_X)
    new_y <- as.matrix(trans_Y(ori_y))
    
    new_beta_hat <- stat(new_X,new_y)$beta
    mu_hat <- mean(ori_y) - apply(log_X,2,mean)%*%new_beta_hat
    
    beta_hat <- c(intercept=mu_hat, new_beta_hat)
    return(beta_hat)
  }
}



