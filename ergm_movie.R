library(network)
library(ergm)

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


eta_distributional <- function(g_obs, feature_func, sample_per_obs, ref_formula, 
	                           learning_iter, gamma){
	#### Not Filled Out ###
	# High Level Idea
	# Takes in list of nets 'g_obs', a function that turns obs into features,
	# a number samples to simulate per obs, a reference formula for simulate,
	# number of learning iterations, and gamma for smoothing obs and sample mean
}


eta_traditional <- function(g_obs, feature_func, sample_per_obs, ref_formula, 
	                           learning_iter, par_gamma){
	#### Not Filled Out ###
	# Similar to above
}