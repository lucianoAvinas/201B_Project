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
                               etak, l_iter, gamma){
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
    for(i in 1:l_iter){
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


sequential_ergm <- function(nets, feat_formula, ref_formula, etak, 
                            l_iter, cd_steps, mle_steps, termination){
    # For stochastic MCMLE set cd_steps to 1 and mle_steps to few
    # ergm will stop if cd if iteration = cd_steps
    etak <- etak / norm(etak, type='2')

    n <- length(nets)
    step <- 0

    cat('(iter, step, net) eta = (',etak,')','\n')
    for (i in 1:l_iter) {
        inds <- sample(1:length(nets))
        eta_prev <- etak

        for (j in inds) {
            step <- step + 1
            etak <- suppress_messages(ergm(update(nets[[j]]~.,feat_formula), 
                        coef=etak, response='e_weights', reference=ref_formula,
                        control=control.ergm(
                        CD.maxit=cd_steps, MCMLE.maxit=mle_steps, 
                        MCMLE.termination=termination)))$coef
            etak <- etak / norm(etak, type='2')
        }      
        cat('    (',i,step,j,')   eta = (',etak,')','\n')
        cat('     eta diff',norm(etak-eta_prev, type='2'),'\n')
    }
}


mom_sgd_ergm <- function(nets, feat_formula, ref_formula, 
                         etak, l_iter, beta, lr_step){
    # For stochastic MCMLE set cd_steps to 1 and mle_steps to few
    # ergm will stop if cd if iteration = cd_steps
    etak <- etak / norm(etak, type='2')

    n <- length(nets)
    step <- 0
    sp <- paste(rep(' ', 2), collapse='')

    cat('(iter, step, net) eta = (',etak,')','\n')
    mom <- 0 # will be vector later by broadcasting
    for (i in 1:l_iter) {
        inds <- sample(1:n)
        eta_prev <- etak
        for (j in inds) {
            step <- step + 1 
            sgd <- suppress_messages(ergm(update(nets[[j]]~.,feat_formula), 
                          coef=etak, response='e_weights', reference=ref_formula,
                          control=control.ergm(CD.maxit=1, MCMLE.maxit=1, 
                          MCMLE.termination='none')))$gradient
            mom <- beta*mom + sgd
            etak <- etak + lr_step*mom

            etak <- etak / norm(etak, type='2')
        }
        cat(sp,' (',i,step,j,')',sp,' eta = (',etak,')','\n')
        cat(sp,sp,'eta diff =',norm(etak-eta_prev, type='2'),'\n')
    }
}


mombatch_gd_ergm <- function(nets, feat_formula, ref_formula, 
                             etak, l_iter, beta, lr_step, batch_sz){
    # For stochastic MCMLE set cd_steps to 1 and mle_steps to few
    # ergm will stop if cd if iteration = cd_steps
    etak <- etak / norm(etak, type='2')

    n <- length(nets)
    step <- 0
    sp <- paste(rep(' ', 2), collapse='')

    cat('(iter, step, net) eta = (',etak,')','\n')
    mom <- 0 # will be vector later by broadcasting
    batched_gd <- 0 # will be vector later by broadcasting    
    for (i in 1:l_iter) {
        inds <- sample(1:n)
        eta_prev <- etak
        for (j in inds) {
            step <- step + 1
            gd <- suppress_messages(ergm(update(nets[[j]]~.,feat_formula), 
                          coef=etak, response='e_weights', reference=ref_formula,
                          control=control.ergm(CD.maxit=1, MCMLE.maxit=1, 
                          MCMLE.termination='none')))$gradient
            batched_gd <- batched_gd + gd

            if (step%%batch_sz == 0) {
                mom <- beta*mom + batched_gd/batch_sz
                etak <- etak + lr_step*mom

                etak <- etak / norm(etak, type='2')
                batched_gd <- 0
            }
        }
        cat(sp,' (',i,step,j,')',sp,' eta = (',etak,')','\n')
        cat(sp,sp,'eta diff =',norm(etak-eta_prev, type='2'),'\n')
    }
}


example_methods <- function(method, extra.args=NULL) {
    # Method takes int \in (0,1,2,3,4)

    test.nets <- nets_in_genre('Drama')[1:32]
    for (i in 1:length(test.nets)) {
        char_groups <- cut(test.nets[[i]] %v% 'ranks', c(1,4,8,Inf), 
                           include.lowest=T, ordered.result=T)
        test.nets[[i]] %v% 'char_groups' <- as.integer(char_groups)
    }

    # Absdiffcat takes categorical variables of 3 categories, so it needs 3 - 1 etas
    if (method == 0) {
        res <- eta_distributional(test.nets, 1024, 
            ~absdiffcat('char_groups')+sum+edges,
            ~Poisson, rnorm(4), 20, 0.5)
    } else if (method == 1) {
        res <- eta_traditional(test.nets, 1024, 
            ~absdiffcat('char_groups')+sum+edges,
            ~Poisson, rnorm(4), 20, 500)
    } else if (method == 2) {
        if (is.null(extra.args)) {
            cd_steps = 60
            mle_steps = 20
            termination = 'Hummel'
        } else {
            cd_steps = extra.args[[1]]
            mle_steps = extra.args[[2]]
            # termination = 'none' to use all mle steps
            # example list(1,3,'none')
            termination = extra.args[[3]]

        }
        res <- sequential_ergm(test.nets,
            ~absdiffcat('char_groups')+sum+edges,
            ~Poisson, rnorm(4), 20, cd_steps, mle_steps, termination)
    } else if (method == 3) {
        res <- mom_sgd_ergm(test.nets,
            ~absdiffcat('char_groups')+sum+edges,
            ~Poisson, rnorm(4), 20, 0.99, 1e-2) 
    } else {
        res <- mombatch_gd_ergm(test.nets,
            ~absdiffcat('char_groups')+sum+edges,
            ~Poisson, rnorm(4), 20, 0.99, 1e-2, 4) 
    }

    res
}

# Try plotting loglik to see effect of mom
#example_methods(3)