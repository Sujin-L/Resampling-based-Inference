permute_test <- function(X_screen, y, rep, level){
  
  X <- X_screen
  
  #### reference value #### 
  ref <- stat(X,y)
  beta_hat <- ref$beta ; t_ref <- ref$t; M <- ref$M; inv <- ref$inv; Q <- ref$Q; A <- ref$A
  
  #### permutation ####
  beta_per <- t_per <- NULL 
  
  for(i in 1:rep){
    
    s <- sample(1:length(y))
    per_y <- y[s] ## permute only y
    stat_per <- stat(X, per_y)
    
    ## permuted beta
    sol <- stat_per$beta
    
    ## permuted t-statistic
    t_star <- stat_per$t
    beta_per <- cbind(beta_per, sol); t_per <- cbind(t_per, t_star)
  }
  
  p_value_t <- sig_var_t <- p_value_beta <- sig_var_beta <-  NULL
  
  for(i in 1:length(t_ref)){
    
    t_per_i <- t_per[i,]; t_ref_i <- t_ref[i]
    pt <- mean( abs(t_per_i) >= abs(t_ref_i) )
    
    if(pt <= level){
      p_value_t <- rbind(p_value_t, c(pt, "***"))
    }
    else{
      p_value_t <- rbind(p_value_t, c(pt, ""))
    }
    
  }
  
  p_value_t <- cbind(round(beta_hat,4), p_value_t)
  
  
  
  for(i in 1:length(beta_hat)){
    b_per_i <- beta_per[i,]; beta_hat_i <- beta_hat[i]
    pbeta <- mean( abs(b_per_i) >= abs(beta_hat_i) )
    
    if(pbeta <= level){
      p_value_beta <- rbind(p_value_beta, c(pbeta, "***"))
    }
    else{
      p_value_beta <- rbind(p_value_beta, c(pbeta, ""))
    }
    
  }
  
  p_value_beta <- cbind(round(beta_hat,4), p_value_beta)
  
  
  colnames(p_value_t) <- colnames(p_value_beta) <- c("beta_hat","p_value","significance")
  rownames(p_value_t) <- rownames(p_value_beta) <- names(beta_hat)
  p_value_t <- as.data.frame(p_value_t); p_value_beta <- as.data.frame(p_value_beta)
  
  ind_sig_t <- rownames(p_value_t[which(as.numeric(as.vector(p_value_t$p_value)) <= level),])
  sig_var_t <- as.data.frame(p_value_t[ind_sig_t,]$beta_hat)
  colnames(sig_var_t) <- ""; rownames(sig_var_t) <- ind_sig_t
  
  ind_sig_beta <- rownames(p_value_beta[which(as.numeric(as.vector(p_value_beta$p_value)) <= level),])
  sig_var_beta <- as.data.frame(p_value_beta[ind_sig_beta,]$beta_hat)
  colnames(sig_var_beta) <- ""; rownames(sig_var_beta) <- ind_sig_beta
  
  
  result <- list(beta_hat= beta_hat, p_value_t = p_value_t, sig_var_t = sig_var_t, p_value_beta = p_value_beta, sig_var_beta = sig_var_beta)
  return(result)
}









bootstrap_test <- function(X_screen, y, rep, level){
  
  X <- X_screen
  
  #### reference value #### 
  ref <- stat(X,y)
  beta_hat <- ref$beta ; t_ref <- ref$t; se <- ref$se
  
  data <- cbind(y, X); colnames(data)[1] <- "y"
  
  #### bootstrap testing - residual resampling ####
  per_boot <- stud_boot <- NULL 
  
  
  for(i in 1:rep){
    
    s <- sample(1:length(y),size=length(y), replace = TRUE)
    boot <- data[s,]
    boot_y <- boot[,1]; boot_X <- boot[,2:ncol(boot)]
    
    stat_boot <- stat(boot_X, boot_y)
    
    per_boot <- cbind(per_boot, as.vector(stat_boot$beta))
    t_boot <- (stat_boot$beta-beta_hat)/stat_boot$se
    stud_boot <- cbind(stud_boot, t_boot)
  }
  
  per_quant <- stud_quant <-  per_quant1 <- stud_int1 <- NULL
  
  for(j in 1:length(t_ref)){
    per_quant <- rbind(per_quant, quantile(per_boot[j,], c(0.025,0.975),na.rm=TRUE))
    stud_quant <- rbind(stud_quant, quantile(stud_boot[j,], c(0.025,0.975), na.rm = TRUE))
  }
  
  per_quant <- cbind(round(beta_hat,3), round(per_quant,3))
  stud_quant <- cbind(round(beta_hat,3), round(stud_quant,3))
  
  for(j in 1:nrow(per_quant)){
    
    if(sign(as.numeric(per_quant[j,2]))==sign(as.numeric(per_quant[j,3]))){
      per_quant1 <- rbind(per_quant1, c(per_quant[j,],"***"))
    }
    else{
      per_quant1 <- rbind(per_quant1, c(per_quant[j,],""))
    }
  }
  
  #per_quant1 <- per_quant
  rownames(per_quant1) <- names(beta_hat)
  colnames(per_quant1)[1] <- colnames(stud_quant)[1]<- "beta_hat"
  colnames(per_quant1)[4] <- "significance"
  per_quant1 <- as.data.frame(per_quant1); stud_quant <- as.data.frame(stud_quant)
  
  
  ind_sig_per <- rownames(per_quant1[which(sign(as.numeric(as.vector(per_quant1[,2])))==sign(as.numeric(as.vector(per_quant1[,3])))),])
  sig_var_per <- as.data.frame(per_quant1[ind_sig_per,]$beta_hat)
  colnames(sig_var_per) <- ""; rownames(sig_var_per) <- ind_sig_per
  
  stud_int <- matrix(0,nrow(stud_quant),2) 
  for(i in 1:nrow(stud_quant)){
    stud_int[i,1] <- beta_hat[i] - se[i]*stud_quant[i,3]
    stud_int[i,2] <- beta_hat[i] - se[i]*stud_quant[i,2]
  }
  
  rownames(stud_int) <- names(beta_hat)
  stud_int <- cbind(round(beta_hat,3), round(stud_int,3))
  
  for(j in 1:nrow(stud_int)){
    
    if(sign(as.numeric(stud_int[j,2]))==sign(as.numeric(stud_int[j,3]))){
      stud_int1 <- rbind(stud_int1, c(stud_int[j,],"***"))
    }
    else{
      stud_int1 <- rbind(stud_int1, c(stud_int[j,],""))
    }
  }
  
  stud_int1 <- as.data.frame(stud_int1); rownames(stud_int1) <- names(beta_hat)
  colnames(stud_int1) <- c("beta_hat","2.5%","97.5%","significance")
  
  ind_sig_stud <- rownames(stud_int1[which(sign(as.numeric(as.vector(stud_int1[,2])))==sign(as.numeric(as.vector(stud_int1[,3])))),])
  sig_var_stud <- as.data.frame(stud_int1[ind_sig_stud,]$beta_hat)
  colnames(sig_var_stud) <- ""; rownames(sig_var_stud) <- ind_sig_stud
  
  result <- list(beta_hat= beta_hat, t_ref = t_ref, per_int = per_quant1, stud_int = stud_int1, 
                 sig_var_per = sig_var_per, sig_var_stud = sig_var_stud)
  return(result)
  
}








t_test <- function(X_screen, X, y, level){
  
  ref <- stat(X_screen,y)
  
  beta_hat <- ref$beta ; t_ref <- ref$t; M <- ref$M; inv <- ref$inv; Q <- ref$Q; A <- ref$A
  tt <- cbind(round(beta_hat,4),round(t_ref,4),2*pt(abs(t_ref),df=nrow(X)-length(screen),lower.tail = FALSE))
  
  colnames(tt) <- c("beta_hat","t value", "p_value"); tt <- as.data.frame(tt)
  
  
  ind_sig <- rownames(tt[which(as.numeric(as.vector(tt$p_value)) <= level),])
  sig_var <- as.data.frame(tt[ind_sig,]$beta_hat)
  colnames(sig_var) <- ""; rownames(sig_var) <- ind_sig
  
  beta_hat_index <- which(colnames(X)%in%rownames(sig_var))
  beta_hat <- rep(0,ncol(X))
  beta_hat[beta_hat_index] <- as.numeric(as.vector(unlist(sig_var)))
  
  p_t <- NULL
  for(j in 1:nrow(tt)){
    if(tt$p_value[j] <= level){
      p_t <- rbind(p_t, c(round(tt[j,],3), "***"))
    }
    else{
      p_t <- rbind(p_t, c(round(tt[j,],3) ,""))
    }
  }
  
  p_t <- as.data.frame(p_t)
  colnames(p_t)[4] <- "significance"; rownames(p_t) <- names(beta_hat)
  
  result <- list(p_value = p_t, sig_var = sig_var)
  return(result)
}

