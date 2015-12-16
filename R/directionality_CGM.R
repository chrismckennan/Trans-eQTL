###The purpose of this script is to run a Gibbs sampler to determine which neighbors of a gene g with cis-eQTL s are unaffected and which are indirectly affected
###Created 151211

###Input: sufficient statistic, direct effect index, number of independent observations, sigma.a, simgma.a weights, Wishart Prior m, theta
#If beta = F, then the theta[1] is the prior probability a node is unaffected

###Output: Posterior probability that each index is an INDIRECT effect; The NA is the direct effect


Gibbs.dir <- function(n.iter, n.burn, suff.stat, D.gam, n.ind, sigma.a, weights.sigma, m, theta, beta=T) {			#P(Node is Unaffected) ~ Beta(theta[1], theta[2])
	n.nei <- dim(suff.stat$SYY)[1] - 1
	
	ind.infer <- rep(0, (n.nei+1))			##0 is unaffected, 1 is indirectly affected, 2 is directly affected
	ind.infer[D.gam] <- 2
	
	iterate.index <- which(ind.infer != 2)
	post.mean <- rep(0, n.nei+1); post.mean[D.gam] <- NA				#Posterior mean that index is indirectly affected
	
	for (i in 1:n.iter) {
		
		for (j in iterate.index) {
			tmp.0 <- ind.infer
			tmp.0[j] <- 0					##Make jth index 0, i.e. unaffected
			
			tmp.1 <- ind.infer
			tmp.1[j] <- 1					##Make jth index 1, i.e. indirectly affected
			
			log10bfU <- Log10BF.sing(suff.stat, D.gam, which(tmp.0 == 0), n.ind, sigma.a, weights.sigma, m)$log10BF
			log10bfI <- Log10BF.sing(suff.stat, D.gam, which(tmp.1 == 0), n.ind, sigma.a, weights.sigma, m)$log10BF
			
			if (beta) {
				pU.j <- ( theta[1] + length(which(tmp.1 == 0)) )/(sum(theta) + n.nei - 1)			#If we use a beta prior, find p(Node j is Unaffected | Rest of Node Assignments)
			} else {
				pU.j <- theta[1]
			}
			
			C.j <- max(log10bfU, log10bfI)
			post.probU.j <- pU.j*10^( log10bfU - C.j )/( pU.j*10^( log10bfU - C.j ) + (1-pU.j)*10^( log10bfI - C.j ) )
			
			is.I.j <- rbinom(1, 1, 1 - post.probU.j)			#Gibbs sampler
			ind.infer[j] <- is.I.j
			if (i > n.burn) {
				post.mean[j] <- post.mean[j] + is.I.j/(n.iter - n.burn)
			}		
		}
		
	}
	
	return(list(post.mean=post.mean, D.ind=D.gam))
}