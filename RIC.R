## Resampling-based inferences for compositional regression when sample sizes are limited 
# option = "1se", "min"

source("preprocess.R")
source("cv_glm.R")
source("stat.R")
source("screening.R")
source("re_estimation.R")
source("testing.R")


RIC <- function(X, y, option , level = 0.05){
  
  ## data preprocessing 
  X_o <- X; y_o <- y # save original value 
  X <- trans_X(X)
  y <- trans_Y(y)
  
  #### screening ####
  # output: graph of lambda
  
  # 10 fold cross-validation -> repetition for all combinations 
  collect <- cv.glm(K = 10, X, y)
  
  #### screening covariates #### 
  screen <- beta_screen(collect, X, y, option)
  X_screen <- as.matrix(X[, screen$beta_list]) ## screened X
  
  #### model fit ####
  
  
  # output: five suggesting models -> p-value
  
  permute <- permute_test(X_screen, y, 1000, level)
  boot <- bootstrap_test(X_screen, y, 1000, level)
  ttest <- t_test(X_screen, X, y, level)
  
  p_value <- list(permute_t = permute$p_value_t, permute_beta = permute$p_value_beta,
                  percentile_bootstrap = boot$per_int, studentized_bootstrap = boot$stud_int,
                  ttest = ttest$p_value)
  
  
  
  #### final models ####
  # output: five final models after re-estimation 
  beta_hat1 <- re_estimation(permute$sig_var_t, X_o, y_o)
  beta_hat2 <- re_estimation(permute$sig_var_beta, X_o, y_o)  
  beta_hat3 <- re_estimation(boot$sig_var_per, X_o, y_o)    
  beta_hat4 <- re_estimation(boot$sig_var_stud, X_o, y_o)
  beta_hat5 <- re_estimation(ttest$sig_var, X_o, y_o)
  
  
  final_model <- list(permute_t = beta_hat1 , permute_beta = beta_hat2, 
                      percentile_bootstrap = beta_hat3, studentized_bootstrap = beta_hat4,
                      ttest = beta_hat5)
  
  
  result <- list(lambda = screen$lambda_star, plot = screen$plot, 
                 p_value = p_value, final_model = final_model)
  
  return(result)
  
}

