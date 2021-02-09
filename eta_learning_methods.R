eta_learning_original <- function(MCMC_sample_size, learning_iter, etak, g_obs){
  # Normalizing eta won't change probability (see normalization term)
  etak <- etak / sqrt(sum(etak^2))
  
  # Step 2 of algorithm: generate MCMC sample using eta0
  eta0 <- etak
  MCMC_list <- simulate(network(18) ~ edges + mutual, coef=eta0, nsim=MCMC_sample_size)
  MCMC_list_features <- summary(MCMC_list ~ edges + mutual);rm(MCMC_list)
  p <- dim(MCMC_list_features)[2]
  
  for(iter in 1:learning_iter){    
    # initialize place holder
    wi_vector <- numeric(MCMC_sample_size)
    wi_ZiZi <- I_hat <- outer(numeric(p),numeric(p))
    wi_Zi <- numeric(p)
    
    # trick: avoid wi to explode
    for(i in 1:MCMC_sample_size){
      wi_vector[i] <- c((etak-eta0) %*% MCMC_list_features[i,])
    }
    wi_vector <- exp(wi_vector - max(wi_vector))
    wi_vector <- wi_vector/sum(wi_vector)
    
    # construct input for equation (3.5)
    for(i in 1:MCMC_sample_size){
      Zi <- MCMC_list_features[i,]
      wi_Zi <- wi_Zi + wi_vector[i]*Zi
    }
    
    # equation (3.5)
    for(i in 1:MCMC_sample_size){
      Zi <- MCMC_list_features[i,]
      # NOTE: outer(wi_Zi,wi_Zi) is a constant here and being applied
      #       MCMC_sample_size number of times
      I_hat <- I_hat + wi_vector[i]*outer(Zi,Zi) - outer(wi_Zi,wi_Zi)
    }
    
    # equation (3.4)
    etak <- c(etak + solve(I_hat) %*% (g_obs - wi_Zi))
    etak <- etak / sqrt(sum(etak^2))
    
  } # end of eta learning
  
  cat('eta =', etak)
  
  return(etak)
}


eta_learning_compact <- function(MCMC_sample_size, learning_iter, etak, g_obs){
  # Normalizing eta won't change probability (see normalization term)
  etak <- etak / sqrt(sum(etak^2))
  
  # Step 2 of algorithm: generate MCMC sample using eta0
  eta0 <- etak
  MCMC_list <- simulate(network(18) ~ edges + mutual, coef=eta0, nsim=MCMC_sample_size)
  Z <- summary(MCMC_list ~ edges + mutual);rm(MCMC_list)
  
  for(iter in 1:learning_iter){
    w <- Z %*% (etak-eta0)
    w <- exp(w - max(w))
    w <- w/sum(w)

    wZ <- array(w, dim(Z)) * Z
    wZ_colsum <- colSums(wZ)
    wZZ_sum <- t(wZ) %*% Z

    I_hat <- wZZ_sum - (matrix(wZ_colsum, ncol=1) %*% wZ_colsum)
      
    etak <- c(etak + solve(I_hat) %*% (g_obs - wZ_colsum))
    etak <- etak / sqrt(sum(etak^2))
    
  } # end of eta learning
  
  return(etak)
}


eta_learning_distributional <- function(learning_iter, MCMC_size, g_obs, par_gamma){
    eta <- c(0,0)
    llr_sum <- -choose(18,2) * log(2)
    for(l_iter in 1:learning_iter){
        MCMC_list <- simulate(network(18) ~ edges + mutual, coef=eta, nsim=MCMC_size)
        MCMC_list_features <- summary(MCMC_list ~ edges + mutual)
        #MCMC_list <- gen_MCMC(MCMC_size, eta)
        #MCMC_list_features <- gen_features_MCMC(MCMC_list)
        xi_sample_mean <- colMeans(MCMC_list_features)
        cov_matrix <- cov(MCMC_list_features)
        xi_hat <- par_gamma * g_obs + (1-par_gamma) * xi_sample_mean
        eta_old <- eta
        eta <- c(eta + solve(cov_matrix) %*% (xi_hat - xi_sample_mean))
        cat('iter =',l_iter,', eta = (',eta,')','\n')
    }
    eta_se <- diag(solve(var(MCMC_list_features)))
    cat('\n')
    cat('eta:',eta,'\n')
    cat('standard error:',eta_se,'\n')
    return(MCMC_list_features)
}