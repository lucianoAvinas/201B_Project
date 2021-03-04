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


fisher_info <- function(eta, eta0, Z_sims) {
    w <- Z_sims %*% (eta-eta0)
    w <- exp(w - max(w))
    w <- w/sum(w)

    wZ <- array(w, dim(Z_sims)) * Z_sims
    wZ_colsum <- colSums(wZ)
    wZZ_sum <- t(wZ) %*% Z_sims

    I_hat <- wZZ_sum - (matrix(wZ_colsum, ncol=1) %*% wZ_colsum)
    f_info <- solve(I_hat)

    f_info
}


mom_sgd_ergm <- function(nets, feat_formula, ref_formula, eta, 
                         l_iter, beta, lr_func, use_print, expect_NA){
    eta <- eta / norm(eta, type='2')

    n <- length(nets)
    step <- 0
    sp <- paste(rep(' ', 2), collapse='')

    if (use_print){
        cat('(iter, step, net) eta = (',eta,')','\n')
    }

    mom <- 0 # will be vector later by broadcasting
    for (i in 1:l_iter) {
        inds <- sample(1:n)
        for (j in inds) {
            eta_prev <- eta

            step <- step + 1 
            sgd <- suppress_messages(ergm(update(nets[[j]]~.,feat_formula), 
                          coef=eta, response='e_weights', reference=ref_formula,
                          control=control.ergm(CD.maxit=1, MCMLE.maxit=1, 
                          MCMLE.termination='none')))$gradient
            if (expect_NA > 0) {
                sgd <- sgd[-expect_NA]
            }
            more_NA <- is.na(sgd)
            if (any(more_NA)) {
                NA_inds <- which(more_NA)
                sgd[NA_inds] <- 0
                cat('Some NA at indices', NA_inds, '\n')
            }
            mom <- beta*mom + sgd
            eta <- eta + lr_func(i)*mom

            eta <- eta / norm(eta, type='2')
            if (use_print) {
                #cat(sgd,'\n')
                cat(sp,' (',i,step,j,')',sp,' eta = (',eta,')','\n')
                cat(sp,sp,'eta diff =',norm(eta-eta_prev, type='2'),'\n')
            }
        }
    }
    list(eta=eta, eta_prev=eta_prev)
}


eta_distributional <- function(nets, feat_formula, ref_formula, eta, 
                               l_iter, gamma, sample_per_obs, use_print){
    eta <- eta / norm(eta, type='2')
    cat('initial eta = (',eta,')','\n')

    # Product of multiples obsv. in log turns into a sum
    # Going through the math will ultimately have division by 
    # n at cov matrix step, so might as well take average now
    g_obs <- rep(0,length(eta))
    for(i in 1:length(nets)){
        g_obs <- g_obs + summary(update(nets[[i]]~., feat_formula),
                                        response='e_weights')
    }
    g_obs <- g_obs/length(nets)

    # New sampling done for each iteration
    for(i in 1:l_iter){
        Z <- simulate_obs(nets, sample_per_obs, feat_formula, ref_formula, 
                          eta)
        eta_prev <- eta

        xi_sample_mean <- colMeans(Z)
        cov_matrix <- cov(Z)
        xi_hat <- gamma * g_obs + (1-gamma) * xi_sample_mean

        eta <- c(eta + solve(cov_matrix) %*% (xi_hat - xi_sample_mean))
        eta <- eta / norm(eta, type='2')

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


eta_MCMC_Handcock <- function(eta, eta_prev, Z_sims) {
    # Assumes MCMC draws are uncorrelated
    # Need to consider auto-cov lag if not
    m <- nrow(Z_sims)

    Ihat_inv <- fisher_info(eta, eta_prev, Z_sims)
    expZ <- sum(Z_sims %*% (eta_prev - eta))/m # eta_diff is reverse here
    var_tilde <- expZ^2 * (Ihat_inv %*% Ihat_inv) / m

    var_tilde
}

eta_MCMCvar <- function(m_sims, Z_sims) {
    # Assumes MCMC draws are uncorrelated
    # Need to consider auto-cov lag if not
    n <- as.integer(nrow(Z_sims)/m_sims)
    var_tilde <- 0
    for (i in 1:n) {
        Z_sel <- Z_sims[(1+(i-1)*m_sims):(i*m_sims),]
        var_tilde <- var_tilde + solve(var(Z_sel))
    }

    var_tilde/n
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
                     sims_per_val) {
    # Compare to zero eta
    n <- length(test_nets)
    test_nets <- net_initializer(test_nets)

    Z_nets <- matrix(rep(0, n*length(eta_pair[[1]])), nrow=n)
    for(i in 1:n){
        Z_nets <- summary(update(test_nets[[i]]~., feat_formula),
                                 response='e_weights')
    }


    Z_sims <- simulate_obs(test_nets, sims_per_val, feat_formula, 
                           ref_formula, eta_pair$eta)
    
    # llk increases with number of samples, better to look at average
    llk_avg <- robust.llk.rel(eta_pair$eta, rep(0, length(eta_pair$eta)), Z_nets, Z_sims)
    eta_var <- eta_MCMCvar(sims_per_val, Z_sims)
    SE <- sqrt(diag(eta_var))
    pvals <- 2*pnorm(-abs(eta_pair$eta/SE))

    cat('llk_avg:', llk_avg, '\n')
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
                        sims_per_val)

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
                         feat_formula, ref_formula, sims_per_test)

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

    nets <- nets_in_genre(genre)
    set.seed(5)
    nets <- nets[sample(length(nets))]

    nets <- net_initializer(nets[(n_train+1):(n_train+n_val)])
    Z_sims <- simulate_obs(nets, sims_per_val, feat_formula, 
                           ref_formula, eta)

    Z_avg <- 0 
    for(i in 1:length(nets)){
        Z_avg <- Z_avg + summary(update(nets[[i]]~., feat_formula),
                                        response='e_weights')
    }
    Z_avg <- Z_avg/length(nets)

    variable_names <- labels(terms(feat_formula))
    k <- length(variable_names)
    height <- round(sqrt(k))
    par(mfrow=c(height, 2*round(k/height)))
    for (i in 1:k) {
        Z_diff <- Z_sims[,i]-Z_avg[i]
        plot(1:nrow(Z_sims), Z_diff,type='l', xlab=NULL, main=paste(
             'MCMC Convergence:', variable_names[[i]]))
        abline(h=mean(Z_diff),col='red',lwd=2)
        plot(density(Z_diff), main=paste('Density plot:', 
                                         variable_names[[i]]))
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
fit_and_val('Drama', 'drama_fit', 50, 20, 1000, 'momentum', add_groups,
            list(feat_formula=~nodemix('char_groups')+
                               atmost(threshold=2)+
                               nodesqrtcovar(center=TRUE), 
                 ref_formula=~Poisson, eta=rnorm(7), l_iter=5, 
                 beta=0.99, lr_func=function(x) 1e-2*exp(-(x-1)/2), 
                 use_print=TRUE, expect_NA=1))
#view_val_sims('test')

#fit_and_val('Drama', 'test', 10, 5, 1000, 'dist', add_groups,
#            list(feat_formula=~nodematch('main_v_misc')+
#                               absdiffcat('char_groups',levels=-3)+
#                               atmost(threshold=2)+
#                               nodesqrtcovar(center=TRUE), 
#                 ref_formula=~Poisson, eta=rnorm(5), l_iter=5, 
#                 gamma=0.5, sample_per_obs=1000, use_print=TRUE))