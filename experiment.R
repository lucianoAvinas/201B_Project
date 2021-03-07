library(network)
library(ergm.count)
library(pkgcond)

load('movie_data.RData')

nets_in_genre <- function(genre){
    genre.info <- movies_by_genre[[genre]]
    genre.nets <- list()
    i <- 1
    for(net.info in genre.info){
        net <- as.network(net.info[['adj_mat']], directed=F, matrix.type='a',
                          ignore.eval=F, names.eval='e_weights')
        net %v% 'ranks' <- net.info[['ranks']]
        net %v% 'names' <- net.info[['names']]

        genre.nets[[i]] <- net
        i <- i + 1
    }
    genre.nets
}


simulate_obs <- function(nets, sample_per_obs, feat_formula,
                         ref_formula, eta) {
    summ_list <- list()
    for (i in 1:length(nets)) {
        sims <- simulate(feat_formula, coef=eta, nsim=sample_per_obs,
                         response='e_weights', reference=ref_formula,
                         basis=nets[[i]])
        summ_list[[i]] <- summary(update(sims~., feat_formula), 
                                  response='e_weights')
    }
    do.call(rbind, summ_list)
}


mom_sgd_ergm <- function(nets, feat_formula, ref_formula, eta, 
                         l_iter, beta, lr_func, use_print, expect_NA){
    #eta <- eta / norm(eta, type='2')

    n <- length(nets)
    step <- 0
    sp <- paste(rep(' ', 2), collapse='')
    #trust <- 75

    if (use_print){
        cat('(iter, step, net) eta = (',eta,')','\n')
    }

    mom <- 0 # will be vector later by broadcasting
    for (i in 1:l_iter) {
        inds <- sample(1:n)
        for (j in inds) {
            eta_prev <- eta

            step <- step + 1 
            if (length(expect_NA) > 0) {
                padded_eta <- add_infs(eta, expect_NA)
            } else {
                padded_eta <- padded_eta
            }
            sgd <- suppress_messages(ergm(update(nets[[j]]~.,feat_formula), 
                          coef=padded_eta, response='e_weights', 
                          reference=ref_formula, control=control.ergm(CD.maxit=1,
                          MCMLE.maxit=1, MCMLE.termination='none')))$gradient
        
            if (length(expect_NA) > 0) {
                sgd <- sgd[-expect_NA]
            }
            more_NA <- is.na(sgd)
            if (any(more_NA)) {
                NA_inds <- which(more_NA)
                sgd[NA_inds] <- 0
                cat('Some NA at indices:', NA_inds, '\n')
            }
            #if (norm(sgd, type='2') > trust) {
            #    cat('Sgd over trust with norm', norm(sgd, type='2'), '\n')
            #    sgd <- 0
            #}
            mom <- beta*mom + sgd
            eta <- eta + lr_func(i)*mom

            eta <- eta #/ norm(eta, type='2')
            if (use_print) {
                #cat(sgd,'\n')
                cat(sp,' (',i,step,j,')',sp,' eta = (',eta,')','\n')
                cat(sp,sp,'eta diff =',norm(eta-eta_prev, type='2'),'\n')
            }
        }
    }

    list(eta=eta, eta_prev=eta_prev)
}


add_infs <- function(vec, inds) {
    temp <- rep(0, length(vec)+length(inds))
    temp[-inds] <- vec
    temp[inds] <- -Inf
    
    temp
}


eta_distributional <- function(nets, feat_formula, ref_formula, eta, 
                               l_iter, gamma, sample_per_obs, use_print,
                               expect_NA){
    #eta <- eta / norm(eta, type='2')
    cat('initial eta = (',eta,')','\n')

    # Product of multiples obsv. in log turns into a sum
    # Going through the math will ultimately have division by 
    # n at cov matrix step, so might as well take average now
    g_obs <- rep(0,length(eta))
    for(i in 1:length(nets)){
        g_obs <- g_obs + summary(update(nets[[i]]~., feat_formula),
                                        response='e_weights')
    }
    if (length(expect_NA) > 0) {
        g_obs <- g_obs[-expect_NA]
    }
    g_obs <- g_obs/length(nets)

    # New sampling done for each iteration
    for(i in 1:l_iter){
        if (length(expect_NA) > 0) {
                padded_eta <- add_infs(eta, expect_NA)
        } else {
            padded_eta = padded_eta
        }
        Z <- simulate_obs(nets, sample_per_obs, feat_formula, ref_formula, 
                          padded_eta)
        if (length(expect_NA) > 0) {
            Z <- Z[,-expect_NA]
        }

        eta_prev <- eta

        xi_sample_mean <- colMeans(Z)
        cov_matrix <- cov(Z)
        xi_hat <- gamma * g_obs + (1-gamma) * xi_sample_mean

        eta <- c(eta + solve(cov_matrix) %*% (xi_hat - xi_sample_mean))
        eta <- eta #/ norm(eta, type='2')

        if (use_print) {
            cat('iter =',i,', eta = (',eta,')','\n')
            cat('eta diff',norm(eta-eta_prev, type='2'),'\n')
        } 
    }
    list(eta=eta, eta_prev=eta_prev)
}


robust.llk.rel <- function(eta, eta_prev, Z_nets, Z_sims){
    # Suggested in Hummel's Improving Simulation.. end of Section 2
    eta_diff <- eta - eta_prev
    etaZ <- Z_sims %*% eta_diff

    mu_hat <- median(etaZ)
    var_hat <- 1.483 * median(abs(etaZ - mu_hat))
    llk_avg <- mean(Z_nets %*% eta_diff) - mu_hat - var_hat/2

    llk_avg
}


ergm_MCMCtranslate <- function(eta, Z_sims) {
    av <- apply(Z_sims, 2, mean)
    Z_sims <- sweep(Z_sims, 2, av, "-")

    basepred <- Z_sims %*% eta
    prob <- max(basepred)
    prob <- exp(basepred - prob)
    prob <- prob/sum(prob)
    g_sim <- Z_sims
    E <- apply(sweep(g_sim, 1, prob, "*"), 2, sum)

    htmp <- sweep(sweep(g_sim, 2, E, "-"), 1, sqrt(prob), "*")
    H <- crossprod(htmp, htmp)
    cat('H', diag(H),'\n')

    cov.zbar <- spectrum0.mvar(g_sim) * sum(prob^2)
    cat('cov.zbar', diag(cov.zbar),'\n')
    mc.cov <- solve(H, cov.zbar, tol=1e-20)
    mc.cov <- solve(H, t(mc.cov), tol=1e-20)
    
    cat('mc.cov',diag(mc.cov),'\n\n')
    mc.cov
} 


eta_MCMC_Handcock <- function(eta, Z_nets, Z_sims) {
    m <- nrow(Z_sims)

    w0 <- Z_sims %*% eta
    w <- exp(w0 - max(w0))
    w <- w/sum(w)

    wZ <- array(w, dim(Z_sims)) * Z_sims
    wZ_colsum <- colSums(wZ)
    wZZ_sum <- t(wZ) %*% Z_sims

    I_hat <- wZZ_sum - (matrix(wZ_colsum, ncol=1) %*% wZ_colsum)
    cat('I_hat', diag(I_hat), '\n')

    W_i <- -sweep(Z_sims, 2, Z_nets, '-')
    #W_i <- array(w, dim(Z_sims)) * W_i

    #V_tilde <- colSums(acf(W_i, type='covariance', plot=F)$acf)/m^2
    V_tilde <- spectrum0.mvar(W_i)/m

    cat('V_tilde', diag(V_tilde), '\n\n')
    #cat('V_oth', diag(V_oth), '\n')
    #mc.cov <- solve(I_hat, V_tilde, tol=1e-20)
    #mc.cov <- solve(I_hat, t(mc.cov), tol=1e-20)/m
    #cat('mc.cov', diag(mc.cov),'\n')
    #mc.cov
    list(I_hat=I_hat, V_tilde=V_tilde)
}


ergm_MCMCcov <- function(m_sims, eta, Z_sims, Z_nets) {
    n <- as.integer(nrow(Z_sims)/m_sims)
    Zsims_av <- 0
    Znets_av <- 0
    #I_tot <- 0
    #V_tot <- 0
    for (i in 1:n) {
        Zsims_av <- Zsims_av + Z_sims[(1+(i-1)*m_sims):(i*m_sims),]
        Znets_av <- Znets_av + Z_nets[i,]
        #res <- eta_MCMC_Handcock(eta, Z_nets[i,], 
        #                         Z_sims[(1+(i-1)*m_sims):(i*m_sims),])
        #I_tot <- I_tot + res$I_hat
        #V_tot <- V_tot + res$V_tilde
    }
    res <- eta_MCMC_Handcock(eta, Znets_av/n, Zsims_av/n)
    I_tot <- res$I_hat
    V_tot <- res$V_tilde
    cat('I_tot', diag(I_tot), '\n')
    cat('V_tot', diag(V_tot), '\n')
    var_tilde <- solve(I_tot, V_tot, tol=1e-20)
    var_tilde <- solve(I_tot, t(var_tilde), tol=1e-20)
    var_tilde
}


fit_eta <- function(train_nets, net_initializer, method, method_args) {
    train_nets <- net_initializer(train_nets)

    if (method == 'momentum') {
        eta_pair <- mom_sgd_ergm(train_nets, method_args$feat_formula, 
                                 method_args$ref_formula, method_args$eta, 
                                 method_args$l_iter, method_args$beta, 
                                 method_args$lr_func, method_args$use_print,
                                 method_args$expect_NA)
    } else {
        eta_pair <- eta_distributional(train_nets, method_args$feat_formula, 
                                       method_args$ref_formula, method_args$eta, 
                                       method_args$l_iter, method_args$gamma, 
                                       method_args$sample_per_obs, 
                                       method_args$use_print)
    } 
    cat('eta:', eta_pair$eta, '\n')
    eta_pair
}


run_test <- function(test_nets, net_initializer, eta_pair, feat_formula, ref_formula,
                     sims_per_val, expect_NA) {
    # Compare to zero eta
    n <- length(test_nets)
    test_nets <- net_initializer(test_nets)
    eta <- eta_pair$eta

    Z_nets <- matrix(rep(0, n*(length(eta_pair[[1]])+length(expect_NA))), nrow=n)
    for(i in 1:n){
        Z_nets[i,] <- summary(update(test_nets[[i]]~., feat_formula),
                                 response='e_weights')
    }


    Z_sims <- simulate_obs(test_nets, sims_per_val, feat_formula, 
                           ref_formula, add_infs(eta, expect_NA))
    Z_sims <- Z_sims[,-expect_NA]
    Z_nets <- Z_nets[,-expect_NA]
    
    # llk increases with number of samples, better to look at average
    llk_avg <- robust.llk.rel(eta, rep(0, length(eta)), Z_nets, Z_sims)
    cat('llk_ratio_rel0:', llk_avg, '\n')
    eta_var <- ergm_MCMCcov(sims_per_val, eta, Z_sims, Z_nets)
    SE <- sqrt(diag(eta_var))
    pvals <- 2*pnorm(-abs(eta/SE))
    #SE <- rep(NA, length(eta))
    #pvals <- rep(NA, length(eta))

    cat('pvals:', pvals, '\n')

    list(llk_avg=llk_avg, SE=SE, pvals=pvals)
}

fit_and_val <- function(genre, save_name, n_train, n_val, sims_per_val,
                        method, net_initializer, method_args) {
    nets <- nets_in_genre(genre)
    set.seed(5)
    nets <- nets[sample(length(nets))]

    train_nets <- nets[1:n_train]
    val_nets <- nets[(n_train+1):(n_train+n_val)]

    eta_pair <- fit_eta(train_nets, net_initializer, method, method_args)
    val_res <- run_test(val_nets, net_initializer, eta_pair, 
                        method_args$feat_formula, method_args$ref_formula, 
                        sims_per_val, method_args$expect_NA)

    all_eta_save_info <- list(eta_pair=eta_pair, val_res=val_res, 
          genre=genre, n_train=n_train, n_val=n_val, method=method, 
          method_args=method_args, sims_per_val=sims_per_val,
          net_initializer=net_initializer)
    save(all_eta_save_info, file=file.path('eta_bank', 
                                 paste0(save_name, '.RData')))

    list(eta_pair=eta_pair, val_res=val_res)
}

test_eta <- function(save_name, n_test, sims_per_test){
    load(file.path('eta_bank', paste0(save_name, '.RData')))
    genre <- all_eta_save_info$genre
    n_train <- all_eta_save_info$n_train
    n_val <- all_eta_save_info$n_val
    net_initializer <- all_eta_save_info$net_initializer

    nets <- nets_in_genre(genre)
    set.seed(5)
    nets <- nets[sample(length(nets))]

    test_nets <- nets[(n_train+n_val+1):(n_train+n_val+n_test)]

    
    feat_formula <- all_eta_save_info$method_args$feat_formula
    ref_formula <- all_eta_save_info$method_args$ref_formula
    eta_pair <- all_eta_save_info$eta_pair
    test_res <- run_test(test_nets, net_initializer, eta_pair, 
                         feat_formula, ref_formula, sims_per_test,
                         all_eta_save_info$method_args$expect_NA)

    test_res
}

view_val_sims <- function(save_name) {
    ## Jumps in MCMC correspond to different net sims
    ## We want correct behavior on average
    ## abline is this average
    load(file.path('eta_bank', paste0(save_name, '.RData')))

    eta <- all_eta_save_info$eta_pair$eta
    genre <- all_eta_save_info$genre
    feat_formula <- all_eta_save_info$method_args$feat_formula
    ref_formula <- all_eta_save_info$method_args$ref_formula
    n_train <- all_eta_save_info$n_train
    n_val <- all_eta_save_info$n_val
    sims_per_val <- all_eta_save_info$sims_per_val
    net_initializer <- all_eta_save_info$net_initializer
    expect_NA <- all_eta_save_info$method_args$expect_NA

    nets <- nets_in_genre(genre)
    set.seed(5)
    nets <- nets[sample(length(nets))]

    nets <- net_initializer(nets[(n_train+1):(n_train+n_val)])
    if (length(expect_NA) > 0) {
                padded_eta <- add_infs(eta, expect_NA)
    } else {
        padded_eta = padded_eta
    }
    Z_sims <- simulate_obs(nets, sims_per_val, feat_formula, 
                           ref_formula, padded_eta)

    Z_avg <- 0 
    for(i in 1:length(nets)){
        Z_avg <- Z_avg + summary(update(nets[[i]]~., feat_formula),
                                        response='e_weights')
    }
    Z_avg <- Z_avg/length(nets)

    if (length(expect_NA) > 0) {
        Z_sims <- Z_sims[,-expect_NA]
        Z_avg <- Z_avg[-expect_NA]
    }

    #variable_names <- labels(terms(feat_formula))
    #k <- length(variable_names)
    k <- length(Z_avg)
    height <- round(sqrt(k))
    par(mfrow=c(height, 2*ceiling(k/height)))
    print('red line is average from val observations')
    print('black lines are simulate MCMC. we want similar avg behavior')
    for (i in 1:k) {
        Z_diff <- Z_sims[,i]-Z_avg[i]
        plot(1:nrow(Z_sims), Z_diff,type='l', xlab=NULL, main=paste(
             'MCMC Convergence: eta coef', i))
        abline(h=mean(Z_diff),col='red',lwd=2)
        plot(density(Z_diff), main=paste('Density plot: eta coef', 
                                         i))
    }
}


add_groups <- function(nets) {
    for (i in 1:length(nets)) {
        char_groups <- cut(nets[[i]] %v% 'ranks', c(1,2,6,Inf), 
                           right=F, ordered.result=T)
        nets[[i]] %v% 'char_groups' <- as.integer(char_groups)
        main_v_misc <- cut(nets[[i]] %v% 'ranks', c(1,6,Inf), 
                           right=F, ordered.result=T)
        nets[[i]] %v% 'main_v_misc' <- as.integer(main_v_misc)
    }
    nets
}


set.seed(3)
#fit_and_val('Drama', 'test', 50, 20, 1000, 'momentum', add_groups,
#            list(feat_formula=~nodemix('char_groups')+
#                               atleast(threshold=3)+
#                               nodesqrtcovar(center=TRUE), 
#                 ref_formula=~Poisson, eta=rep(1,7), l_iter=3, 
#                 beta=0, lr_func=function(x) 1e-2*exp(-(x-1)), 
#                 use_print=TRUE, expect_NA=c(1)))
#comedy_fit used 9e-3*exp(-(x-1)) learning rate and eta=c(1,1,1,0.5,1,1,1)

#view_val_sims('drama_fit')

test_eta('test', 20, 2000)