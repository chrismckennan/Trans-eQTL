###The purpose of this script is to run a Gibbs sampler to determine which neighbors of a gene g with cis-eQTL s are unaffected and which are indirectly affected
###Created 151211

###Input: sufficient statistic, direct effect index, number of independent observations, sigma.a, simgma.a weights, Wishart Prior m, theta
#If direchlet = F, then the theta is proportional to the prior probability a node is unaffected, indirectly affected, directly affected by the cis-eQTL

###Output: Posterior probability that each index is UNAFFECTED by the cis-eQTL in the network; we do not care whether a gene is directly or indirectly affected by the cis-eQTL. The NA is the gene associated with the cis-eQTL that perturbs the network

##This was the intitial Gibbs Sampler, which does not assume there is a direct effect between the SNP and the source gene

Gibbs.dir.0 <- function(n.iter, n.burn, suff.stat, D.gam, n.ind, sigma.a, weights.sigma, m, theta=rep(0.5, 3), dirichlet=T) {			#P(Node Label) ~ Dirichlet(theta[1], theta[2], theta[3])
	n.nei <- dim(suff.stat$SYY)[1] - 1
	
	ind.infer <- rbinom((n.nei+1), 1, 1/2)			##0 is unaffected, 1 is indirectly affected, 2 is directly affected
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
	if (! dirichlet) {
		theta <- theta/sum(theta)
		post.prob.nodes <- theta
	}
	
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
		if (i > n.burn && dirichlet) {
			gibbs.probs <- rdirichlet(1, theta + c(length(which(ind.infer == 0)), length(which(ind.infer == 1)), length(which(ind.infer == 2)) - 1) )			#Gibbs sampler for (p_U, p_I, p_D)
			post.prob.nodes <- post.prob.nodes + gibbs.probs/(n.iter - n.burn)		
		}
	}
	
	return(list(post.mean.U=post.mean, post.mean.I=post.mean.I, post.mean.D=post.mean.D, D.ind=D.gam, post.probs=post.prob.nodes))
}


##Return the indices of the n.corr most correlated genes with the source gene (this will include the source gene)
Corr.ind <- function(YtY, source.gene, n.corr) {
	Corr <- diag( 1/sqrt(diag(YtY)) ) %*% YtY %*% diag( 1/sqrt(diag(YtY)) )
	return(order(-abs(Corr[source.gene,]))[1:min((n.corr+1), nrow(YtY))])
}

##Function to create Gibbs sample output. This wil be the function to parallelize when running many Gibbs samplers
#The function writes an output file and will be parrellized
Gibbs.par <- function(file, directory.in, directory.out, max.genes) {
	file.path <- paste(directory.in, file, sep="/")
	all_data = readLines(file.path)
	tissue <- strsplit(all_data[1], split='\t', perl=T)[[1]][2]
	n.ind <- as.numeric(strsplit(all_data[2], split='\t', perl=T)[[1]][2])   #Number of independent measurements
	chr <- as.numeric(strsplit(all_data[3], split='\t', perl=T)[[1]][2])    #Chromosome number
	gene <- strsplit(all_data[4], split='\t', perl=T)[[1]][2]    #Gene of interest
	SNP <- strsplit(all_data[5], split='\t', perl=T)[[1]][2]    #eQTL of interest
	Gene.names <- strsplit(all_data[6], split='\t', perl=T)[[1]][2:length(strsplit(all_data[6], split='\t', perl=T)[[1]])]   #Column names for Y'Y
	n.genes <- length(Gene.names)    #Number of genes in the network
	source.ind.all.g <- which(Gene.names == gene)
	
	YtY <- array(NA, dim=c(n.genes, n.genes))   #Y'Y
	for (r in 1:n.genes) {
		YtY[r,] <- as.numeric(strsplit(all_data[7+r], split='\t', perl=T)[[1]])
	}
	ind.use.g <- Corr.ind(YtY, source.ind.all.g, max.genes)
	YtY.use.g <- YtY[ind.use.g, ind.use.g]
	
	sxx <- as.numeric(strsplit(all_data[8+n.genes], split='\t', perl=T)[[1]][2])   #X'X, a scalar
	YtX.use.g <- as.numeric(strsplit(all_data[10+n.genes], split='\t', perl=T)[[1]])[ind.use.g]    #Y'X, a vector
	
	suff.stat <- list(SYY = YtY.use.g/n.ind, sxx = sxx/n.ind, SYX = YtX.use.g/n.ind, SY1 = rep(0, length(ind.use.g)), mu.g = 0)
	
	n.genes.use <- length(ind.use.g)
	if (n.genes.use <= 15) {
		n.iter <- 3000
		n.burn <- 1000
	} else {
	if (n.genes.use <= 25) {
		n.iter <- 4000
		n.burn <- 2000
	} else {
		n.iter <- 5000
		n.burn <- 2500
		}
	}
	sigma.a <- c(0.1, 0.4)
	weights.sigma <- c(1, 1)
	
	out.file = file.path(directory.out, sub("(summary)", '\\1.gibbs', file, perl=T))
	
	gibbs <- Gibbs.dir.1(n.iter, n.burn, suff.stat, 1, n.ind, sigma.a, weights.sigma, n.genes.use-1)
	pup.all.g <- rep(NA, n.genes)
	pup.all.g[ind.use.g] <- gibbs$post.mean.U     #Posterior probability a gene is UNAFFECTED by the source gene
	
	gibbs.mat <- rbind(pup.all.g)
	colnames(gibbs.mat) <- Gene.names
	write.table(gibbs.mat, out.file, col.names=T, row.names=F, append=F, quote=F, sep="\t")
	pup.all.g[!is.na(pup.all.g)]
}