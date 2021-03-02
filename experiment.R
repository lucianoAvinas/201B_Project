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
	                     l_iter, beta, lr_func, use_print){
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
            mom <- beta*mom + sgd
            eta <- eta + lr_func(i)*mom

            eta <- eta / norm(eta, type='2')
            if (use_print) {
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


eta_MCMCvar <- function(eta, eta_prev, Z_sims) {
	# Assumes MCMC draws are uncorrelated
	# Need to consider auto-cov lag if not
	m <- nrow(Z_sims)

	Ihat_inv <- fisher_info(eta, eta_prev, Zsim_sel)
	expZ <- exp(Zsim_sel %*% eta_diff)/m
	var_tilde <- (expZ %*% expZ) * (Ihat_inv %*% Ihat_inv) / m

	var_tilde
}


fit_eta <- function(train_nets, method, method_args) {
	if (method == 'momentum') {
		eta_pair <- mom_sgd_ergm(train_nets, method_args$feat_formula, 
			                     method_args$ref_formula, method_args$eta, 
	                             method_args$l_iter, method_args$beta, 
	                             method_args$lr_func, method_args$use_print)
	} else {
		eta_pair <- eta_distributional(train_nets, method_args$feat_formula, 
			                           method_args$ref_formula, method_args$eta, 
	                                   method_args$l_iter, method_args$gamma, 
	                                   method_args$sample_per_obs, 
	                                   method_args$use_print)
	}
	
	eta_pair
}


run_test <- function(test_nets, eta_pair, feat_formula, ref_formula,
	                 sims_per_val) {
	n <- length(test_nets)
	Z_nets <- matrix(rep(0, n*eta_pair[[1]]), nrow=n)
    for(i in 1:n){
        Z_nets <- summary(update(test_nets[[i]]~., feat_formula),
                                 response='e_weights')
    }

    # llk estimate only holds for eta0 close to eta
    # We have chosen the second last eta to be our eta0
    # If we have convergence or near convergence this should be ok
    Z_sims <- simulate_obs(test_nets, sims_per_val, feat_formula, 
    	                   ref_formula, eta_pair$eta_prev)
    
    # llk increases with number of samples, better to look at average
    llk_avg <- robust.llk.rel(eta_pair$eta, eta_pair$eta_prev, Z_nets, Z_sims)
    eta_var <- eta_MCMCvar(eta_pair$eta, eta_pair$eta_prev, Z_sims)
    pvals <- 2*pnorm(-abs(eta_pair$eta/sqrt(diag(eta_var))))

    print('llk_avg:', llk_avg)
    print('pvals:', pvals)

    list(llk_avg=llk_avg, SE=SE, pvals=pvals)
}

fit_and_val <- function(genre, n_train, n_val, method, method_args, 
	                    sims_per_val, save_name) {
	nets <- nets_in_genre(genre)
	set.seed(5)
	nets <- nets[sample(length(nets))]

	train_nets <- nets[1:n_train]
	val_nets <- nets[(n_train+1):(n_train+n_val)]

	eta_pair <- fit_eta(train_nets, method, save_name, extra_args)
	val_res <- run_test(val_nets, eta_pair, method_args$feat_formula,
		                method_args$ref_formula, sims_per_val)

	save('all_eta_save_info', list(eta_pair=eta_pair, val_res=val_res, 
		  genre=genre, n_train=n_train, n_val=n_val, method=method, 
		  method_args=method_args, sims_per_val=sims_per_val),
	      file=file.path('eta_bank', save_name))

	list(eta_pair=eta_pair, val_res=val_res)
}

test_eta <- function(n_test, sims_per_test, save_name){
    load(file.path('eta_bank', save_name))
    genre <- all_eta_save_info$genre
    n_train <- all_eta_save_info$n_train
    n_val <- all_eta_save_info$n_val

	nets <- nets_in_genre(genre)
	set.seed(5)
	nets <- nets[sample(length(nets))]
	test_nets <- nets[(n_train+n_val+1):(n_train+n_val+n_test)]

    
    feat_formula <- all_eta_save_info$method_args$feat_formula
    ref_formula <- all_eta_save_info$method_args$ref_formula
    eta_pair <- all_eta_save_info$eta_pair
	test_res <- run_test(test_nets, eta_pair, feat_formula,
		                 ref_formula, sims_per_test)

	test_res
}

view_sims <- function(save_name, net_index, n_sims) {
	load(file.path('eta_bank', save_name))
	eta <- all_eta_save_info$eta_pair$eta
	genre <- all_eta_save_info$genre
	feat_formula <- all_eta_save_info$method_args$feat_formula
    ref_formula <- all_eta_save_info$method_args$ref_formula

	nets <- nets_in_genre(genre)
	net_sims <- simulate(feat_formula, coef=eta, nsim=n_sims,
                         response='e_weights', reference=ref_formula,
                         basis=nets[[net_index]])
    Z_sims <- summary(update(net_sims~., feat_formula), response='e_weights')
    Z_net <- summary(update(nets[[net_index]]~., feat_formula), response='e_weights')

    k <- length(eta)
    height <- as.integer(sqrt(k))
    par(mfrom=c(height, 2*round(k/height)))
    variable_names <- labels(terms(feat_formula))
    for (i in 1:k) {
    	plot(1:n_sims, Z_sims[,i]-Z_net[i],type='l', xlab=NULL, main=paste(
    		 'MCMC Convergence:', variable_names[[i]])
        abline(h=0,col='red',lwd=2)
        plot(density(Z_sims[,i]-Z_net[i]), main=paste('Density plot:', 
        	                                          variable_names[[i]]))
    }
}



