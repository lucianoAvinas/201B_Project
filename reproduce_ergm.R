library(networkdata)
library(ergm)

source('eta_learning_methods.R')

data('sampson')
g_obs <- summary(samplike~edges+mutual)

res <- ergm(samplike~edges+mutual)
cat('\nERGM normalized eta =', res[[1]]/sqrt(sum(res[[1]]^2)),'\n')

eta_k <- c(-2,-2)
reps <- 20
for(i in 1:reps) {
	eta_k <- eta_learning_compact(1000, 500, eta_k, g_obs)
}

cat('\nEta learning eta =', eta_k,'\n')