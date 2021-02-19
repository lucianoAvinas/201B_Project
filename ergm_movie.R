library(network)
library(ergm.count)

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
                         ref_formula, etak) {
    summ_list <- list()
    for (i in 1:length(nets)) {
        sims <- simulate(feat_formula, coef=etak, nsim=sample_per_obs,
                         response='e_weights', reference=ref_formula,
                         basis=nets[[i]])
        summ_list[[i]] <- summary(update(sims~., feat_formula), 
                                  response='e_weights')
    }
    do.call(rbind, summ_list)
}


eta_traditional <- function(nets, sample_per_obs, feat_formula, ref_formula, 
                            etak, n_samplings, gradient_steps){
    etak <- etak / norm(etak, type='2')

    # Product of multiples obsv. in log turns into a sum
    # Going through the math will ultimately have division by 
    # n at info matrix step, so might as well take average now
    g_obs <- rep(0,length(etak))
    for(i in 1:length(nets)){
        g_obs <- g_obs + summary(update(nets[[i]]~., feat_formula),
                                        response='e_weights')
    }
    g_obs <- g_obs/length(nets)

    for (i in 1:n_samplings) {
        Z <- simulate_obs(nets, sample_per_obs, feat_formula, ref_formula, 
                          etak)
        eta0 <- etak
        for (j in 1:gradient_steps) {
            eta_prev <- etak

            w <- Z %*% (etak-eta0)
            w <- exp(w - max(w))
            w <- w/sum(w)

            wZ <- array(w, dim(Z)) * Z
            wZ_colsum <- colSums(wZ)
            wZZ_sum <- t(wZ) %*% Z

            I_hat <- wZZ_sum - (matrix(wZ_colsum, ncol=1) %*% wZ_colsum)
            etak <- c(etak + solve(I_hat) %*% (g_obs - wZ_colsum))

            etak <- etak / norm(etak, type='2')
        }
        cat('iter =',i,', eta = (',etak,')','\n')
        cat('eta diff',norm(etak-eta0, type='2'),'\n')
    }
    # shuffle Z to remove temporal depedence of observations
    list(eta=etak, last_sim=Z[sample(nrow(Z)),])
}


eta_distributional <- function(nets, sample_per_obs, feat_formula, ref_formula, 
                               etak, learning_iter, gamma){
    etak <- etak / norm(etak, type='2')

    # Product of multiples obsv. in log turns into a sum
    # Going through the math will ultimately have division by 
    # n at cov matrix step, so might as well take average now
    g_obs <- rep(0,length(etak))
    for(i in 1:length(nets)){
        g_obs <- g_obs + summary(update(nets[[i]]~., feat_formula),
                                        response='e_weights')
    }
    g_obs <- g_obs/length(nets)

    # New sampling done for each iteration
    for(i in 1:learning_iter){
        Z <- simulate_obs(nets, sample_per_obs, feat_formula, ref_formula, 
                          etak)
        eta0 <- etak

        xi_sample_mean <- colMeans(Z)
        cov_matrix <- cov(Z)
        xi_hat <- gamma * g_obs + (1-gamma) * xi_sample_mean

        etak <- c(etak + solve(cov_matrix) %*% (xi_hat - xi_sample_mean))
        etak <- etak / norm(etak, type='2')

        cat('iter =',i,', eta = (',etak,')','\n')
        cat('eta diff',norm(etak-eta0, type='2'),'\n')
    }
    # shuffle Z to remove temporal depedence of observations
    list(eta=etak, last_sim=Z[sample(nrow(Z)),])
}


nllratio <- function(eta, nets.holdout, feat_formula, ref_formula, eta0=NULL){
    if (is.null(eta0)) {
        eta0 <- rep(0, length(eta))
    }

    all_obs <- list()
    for(i in 1:length(nets.holdout)){
        all_obs[[i]] <- summary(update(nets[[i]]~., feat_formula),
                           response='e_weights')
    }
    all_obs <- do.call(c, rbind)
    sum_obs <- rowSums(all_obs)

    # negative log likelihood
    r_eta_eta0 <- -sum_obs %*% (eta - eta0) + log(mean(exp(all_obs %*% (eta - eta0))))

    r_eta_eta0
}


example_methods <- function(method) {
    # Method takes int = 0 or 1

    test.nets <- nets_in_genre('Drama')[1:10]
    for (i in 1:length(test.nets)) {
        char_groups <- cut(test.nets[[i]] %v% 'ranks', c(1,4,8,Inf), 
                           include.lowest=T, ordered.result=T)
        test.nets[[i]] %v% 'char_groups' <- as.integer(char_groups)
    }

    # Absdiffcat takes categorical variables of 3 categories, so it needs 3 - 1 etas
    if (method == 0) {
        res <- eta_distributional(test.nets, 1000, 
            ~absdiffcat('char_groups')+smallerthan(3)+nodesqrtcovar(TRUE)+sum+edges,
            ~Poisson, rnorm(6), 10, 0.5)
    } else {
        res <- eta_traditional(test.nets, 1000, 
            ~absdiffcat('char_groups')+smallerthan(3)+nodesqrtcovar(TRUE)+sum+edges,
            ~Poisson, rnorm(6), 10, 500)
    }

    res
}


example_methods(0)