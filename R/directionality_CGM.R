###The purpose of this script is to run a Gibbs sampler to determine which neighbors of a gene g with cis-eQTL s are unaffected and which are indirectly affected
###Created 151211

###Input: sufficient statistic, direct effect index, number of independent observations, sigma.a, simgma.a weights, Wishart Prior m, theta
#If direchlet = F, then the theta is proportional to the prior probability a node is unaffected, indirectly affected, directly affected by the cis-eQTL

###Output: Posterior probability that each index is UNAFFECTED by the cis-eQTL in the network; we do not care whether a gene is directly or indirectly affected by the cis-eQTL. The NA is the gene associated with the cis-eQTL that perturbs the network

##This was the intitial Gibbs Sampler, which does not assume there is a direct effect between the SNP and the source gene

Gibbs.dir.0 <- function(n.iter, n.burn, suff.stat, D.gam, n.ind, sigma.a, weights.sigma, m, theta=rep(0.5, 3), dirichlet=T) {			#P(Node Label) ~ Dirichlet(theta[1], theta[2], theta[3])
	n.nei <- dim(suff.stat$SYY)[1] - 1
	
	ind.infer <- rep(0, (n.nei+1))			##0 is unaffected, 1 is indirectly affected, 2 is directly affected
	ind.infer[D.gam] <- 2
	
	#iterate.index <- (1:(n.nei+1))[-D.gam]
	iterate.index <- (1:(n.nei+1))
	post.mean <- rep(0, n.nei+1); #post.mean[D.gam] <- NA				#Posterior mean that index is unaffected by cis-eQTL
	
	for (i in 1:n.iter) {
		
		for (j in iterate.index) {
			tmp.0 <- ind.infer
			tmp.0[j] <- 0					##Make jth index 0, i.e. unaffected
			
			tmp.1 <- ind.infer
			tmp.1[j] <- 1					##Make jth index 1, i.e. indirectly affected
			
			tmp.2 <- ind.infer
			tmp.2[j] <- 2					##Make jth index 2, i.e. directly affected
			
			log10bfU <- Log10BF.sing.all(suff.stat, which(tmp.0 == 2), which(tmp.0 == 0), n.ind, sigma.a, weights.sigma, m)$log10BF
			log10bfI <- Log10BF.sing.all(suff.stat, which(tmp.1 == 2), which(tmp.1 == 0), n.ind, sigma.a, weights.sigma, m)$log10BF
			log10bfD <- Log10BF.sing.all(suff.stat, which(tmp.2 == 2), which(tmp.2 == 0), n.ind, sigma.a, weights.sigma, m)$log10BF
			
			if (dirichlet) {
				pU.j <- ( theta[1] + length(which(ind.infer[-j] == 0)) )/(sum(theta) + n.nei)			#If we use a Dirichlet prior, find p(Node j is Unaffected | Rest of Node Assignments)
				pI.j <- ( theta[2] + length(which(ind.infer[-j] == 1)) )/(sum(theta) + n.nei)
				pD.j <- 1 - pU.j - pI.j
			} else {
				pU.j <- theta[1]
				pI.j <- theta[2]
				pD.j <- theta[3]
			}
			
			C.j <- max(log10bfU, log10bfI, log10bfD)
			post.probU.j <- pU.j*10^( log10bfU - C.j )/( pU.j*10^( log10bfU - C.j ) + pI.j*10^( log10bfI - C.j ) + pD.j*10^( log10bfD - C.j ) )
			post.probI.j <- pI.j*10^( log10bfI - C.j )/( pU.j*10^( log10bfU - C.j ) + pI.j*10^( log10bfI - C.j ) + pD.j*10^( log10bfD - C.j ) )
			
			v.gibbs <- which( rmultinom( 1, 1, c(post.probU.j, post.probI.j, 1 - (post.probI.j + post.probU.j)) ) == 1 ) - 1			#Gibbs sampler
			
			ind.infer[j] <- v.gibbs	
			if (i > n.burn) {
				post.mean[j] <- post.mean[j] + as.numeric(! v.gibbs)/(n.iter - n.burn)
			}		
		}
		
	}
	
	return(list(post.mean=post.mean, D.ind=D.gam))
}


##This version assumes there is a direct effect between the SNP and source gene (i.e. there is no need to include this in the Gibbs sampler, as it is assumed to be known)

Gibbs.dir.1 <- function(n.iter, n.burn, suff.stat, D.gam, n.ind, sigma.a, weights.sigma, m, theta=rep(0.5, 3), dirichlet=T) {			#P(Node Label) ~ Dirichlet(theta[1], theta[2], theta[3])
	n.nei <- dim(suff.stat$SYY)[1] - 1
	
	ind.infer <- rep(0, (n.nei+1))			##0 is unaffected, 1 is indirectly affected, 2 is directly affected
	ind.infer[D.gam] <- 2
	
	iterate.index <- (1:(n.nei+1))[-D.gam]
	
	post.mean <- rep(0, n.nei+1); post.mean[D.gam] <- NA				#Posterior mean that index is unaffected by cis-eQTL
	post.mean.I <- rep(0, n.nei+1); post.mean.I[D.gam] <- NA
	post.mean.D <- rep(0, n.nei+1); post.mean.D[D.gam] <- NA
	
	post.prob.nodes <- rep(0, 3)			#P( (p_U, p_I, p_D) | Data)
	
	for (i in 1:n.iter) {
		
		for (j in iterate.index) {
			tmp.0 <- ind.infer
			tmp.0[j] <- 0					##Make jth index 0, i.e. unaffected
			
			tmp.1 <- ind.infer
			tmp.1[j] <- 1					##Make jth index 1, i.e. indirectly affected
			
			tmp.2 <- ind.infer
			tmp.2[j] <- 2					##Make jth index 2, i.e. directly affected
			
			log10bfU <- Log10BF.sing.all(suff.stat, which(tmp.0 == 2), which(tmp.0 == 0), n.ind, sigma.a, weights.sigma, m)$log10BF
			log10bfI <- Log10BF.sing.all(suff.stat, which(tmp.1 == 2), which(tmp.1 == 0), n.ind, sigma.a, weights.sigma, m)$log10BF
			log10bfD <- Log10BF.sing.all(suff.stat, which(tmp.2 == 2), which(tmp.2 == 0), n.ind, sigma.a, weights.sigma, m)$log10BF
			
			if (dirichlet) {
				pU.j <- ( theta[1] + length(which(ind.infer[-j] == 0)) )/(sum(theta) + n.nei - 1)			#If we use a Dirichlet prior, find p(Node j is Unaffected | Rest of Node Assignments)
				pI.j <- ( theta[2] + length(which(ind.infer[-j] == 1)) )/(sum(theta) + n.nei - 1)
				pD.j <- 1 - pU.j - pI.j
			} else {
				pU.j <- theta[1]
				pI.j <- theta[2]
				pD.j <- theta[3]
			}
			
			C.j <- max(log10bfU, log10bfI, log10bfD)
			post.probU.j <- pU.j*10^( log10bfU - C.j )/( pU.j*10^( log10bfU - C.j ) + pI.j*10^( log10bfI - C.j ) + pD.j*10^( log10bfD - C.j ) )
			post.probI.j <- pI.j*10^( log10bfI - C.j )/( pU.j*10^( log10bfU - C.j ) + pI.j*10^( log10bfI - C.j ) + pD.j*10^( log10bfD - C.j ) )
			
			v.gibbs <- which( rmultinom( 1, 1, c(post.probU.j, post.probI.j, 1 - (post.probI.j + post.probU.j)) ) == 1 ) - 1			#Gibbs sampler
			
			ind.infer[j] <- v.gibbs	
			if (i > n.burn) {
				post.mean[j] <- post.mean[j] + as.numeric(! v.gibbs)/(n.iter - n.burn)
			    post.mean.I[j] <- post.mean.I[j] + as.numeric(v.gibbs == 1)/(n.iter - n.burn)
			    post.mean.D[j] <- post.mean.D[j] + as.numeric(v.gibbs == 2)/(n.iter - n.burn)
			}		
		}
		if (i > n.burn) {
			gibbs.probs <- rdirichlet(1, theta + c(length(which(ind.infer == 0)), length(which(ind.infer == 1)), length(which(ind.infer == 2)) - 1) )			#Gibbs sampler for (p_U, p_I, p_D)
			post.prob.nodes <- post.prob.nodes + gibbs.probs/(n.iter - n.burn)		
		}
	}
	
	return(list(post.mean.U=post.mean, post.mean.I=post.mean.I, post.mean.D=post.mean.D, D.ind=D.gam, post.probs=post.prob.nodes))
}