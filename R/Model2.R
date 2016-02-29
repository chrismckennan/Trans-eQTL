####The purpose of this script is to compute a Bayes Factor for the model where we condition on the network, Sigma
##It was created 160225
##The model is similar to Matthews. We do not integrate out Sigma. We can also compute the Bayes factor using only summary data

#Compute the Bayes factor with summary data input
BF.cgm <- function(suff.stat, I.index, D.index, U.index, n.ind, sigma.a, weights.sigma, lambda=1) {
	n.U <- length(U.index)
	n.D <- length(D.index)
	n.I <- length(I.index)
	n.total <- n.U + n.D + n.I        #Number of genes in the network (including the anchor gene)
	order.vec <- c(I.index, D.index, U.index)
	
	var.a <- sigma.a^2
	weights <- weights.sigma/sum(weights.sigma)
	
	##Unpack sufficient statistics in the order of I, D, U
	Sigma <- suff.stat$Sigma[order.vec, order.vec]    			
	if (n.I > 0) {
		if (n.I == 1) {
			Sigma.12 <- rbind(Sigma[1:n.I, (n.I+1):(n.I+n.D+n.U)])
		} else {
			if (n.I == n.total - 1) {
				Sigma.12 <- cbind(Sigma[1:n.I, (n.I+1):(n.I+n.D+n.U)])
			} else {
				Sigma.12 <- Sigma[1:n.I, (n.I+1):(n.I+n.D+n.U)]
			}
		}
		
		Sigma.22 <- as.matrix(Sigma[(n.I+1):(n.I+n.D+n.U), (n.I+1):(n.I+n.D+n.U)])	
	}
	Sigma.inv <- solve(Sigma)
	sxx <- suff.stat$sxx
	SYX <- suff.stat$SYX[order.vec]
	SY1 <- suff.stat$SY1[order.vec]
	mu.g <- suff.stat$mu.g
	
	##Create vectors and compute traces##
	if (n.I > 0) {
		f <- lambda * Sigma.12 %*% solve(Sigma.22, c(rep(1,n.D), rep(0, n.U)))     #an n.I vector
	} else {
		f <- c()
	}
	f1.vec <- cbind( c(f, rep(1, n.D), rep(0, n.U)) )
	tr.inner <- n.ind * (sxx - mu.g^2) * as.numeric( t(f1.vec) %*% Sigma.inv %*% f1.vec )
	tr.cross <- n.ind * as.numeric( t(f1.vec) %*% Sigma.inv %*% cbind( SYX - mu.g * SY1 ) )
	
	##Compute Bayes Factor##
	logbf.vec <- 1/2 * ( tr.cross^2/(tr.inner + 1/var.a) - log(var.a * tr.inner + 1) )
	C <- max(logbf.vec)
	return( (log( sum( exp(logbf.vec - C) * weights ) ) + C)/log(10) )     #Log10 BF
}

#Compute all possible Bayes Factors
BF.cgm.all <- function(suff.stat, anchor.index, n.ind, sigma.a, weights.sigma, lambda=1) {
	n.total <- length(suff.stat$SYX)
	n.nei <- n.total - 1
	
	BF.mat <- array(2, dim=c(3^n.nei, n.nei + 2))			#Columns correspond to genes in the network; 0: Unaffected, 1: Indirectly affected, 2: Directly affected; The gene corresponding to the cis-eQTL is assumed to be Directly affected by g
	tmp.list <- list()
	for (i in 1:n.nei) {
		tmp.list[[i]] <- c(0,1,2)
	}
	BF.mat[,(1:(n.nei+1))[-anchor.index]] <- as.matrix(expand.grid(tmp.list))		#BF.mat now contains all possible groupings, while keeping the direct effect fixed at 1
	
	for (i in 1:nrow(BF.mat)) {
		U.index <- which(BF.mat[i,1:n.total] == 0)
		D.index <- which(BF.mat[i,1:n.total] == 2)
		I.index <- which(BF.mat[i,1:n.total] == 1)
		BF.mat[i,n.total+1] <- BF.cgm(suff.stat, I.index, D.index, U.index, n.ind, sigma.a, weights.sigma, lambda)
	}
	return(list(BF.mat=BF.mat, D.gam=anchor.index))
}

##This version assumes there is a direct effect between the SNP and source gene (i.e. there is no need to include this in the Gibbs sampler, as it is assumed to be known)
#It is the same code as Gibbs.dir.1, except the Bayes factor calcualtion has now changed
Gibbs.dir.cgm <- function(n.iter, n.burn, suff.stat, D.gam, n.ind, sigma.a, weights.sigma, theta=rep(0.5, 3), dirichlet=T, update.lambda=F) {			#P(Node Label) ~ Dirichlet(theta[1], theta[2], theta[3])
	n.nei <- dim(suff.stat$SYY)[1] - 1
	lambda.old <- 1
	lambda.vec <- rep(0, n.iter - n.burn)
	
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
			
			log10bfU <- BF.cgm(suff.stat, which(tmp.0 == 1), which(tmp.0 == 2), which(tmp.0 == 0), n.ind, sigma.a, weights.sigma, lambda.old)
			log10bfI <- BF.cgm(suff.stat, which(tmp.1 == 1), which(tmp.1 == 2), which(tmp.1 == 0), n.ind, sigma.a, weights.sigma, lambda.old)
			log10bfD <- BF.cgm(suff.stat, which(tmp.2 == 1), which(tmp.2 == 2), which(tmp.2 == 0), n.ind, sigma.a, weights.sigma, lambda.old)
			
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
			if (abs(post.probU.j) < 1e-8) {
				post.probU.j = 0
			}
			post.probI.j <- pI.j*10^( log10bfI - C.j )/( pU.j*10^( log10bfU - C.j ) + pI.j*10^( log10bfI - C.j ) + pD.j*10^( log10bfD - C.j ) )
			if (abs(post.probI.j) < 1e-8) {
				post.probI.j = 0
			}			
			post.probD.j <- 1 - post.probU.j - post.probI.j
			if (abs(post.probD.j) < 1e-8) {
				post.probD.j = 0
			}	
			
			gibbs.sample <- rmultinom( 1, 1, c(post.probU.j, post.probI.j, post.probD.j) )
			v.gibbs <- which( gibbs.sample == 1 ) - 1			#Gibbs sampler
			log10bf.use <- c(log10bfU, log10bfI, log10bfD)[v.gibbs + 1]
			
			ind.infer[j] <- v.gibbs	
			if (i > n.burn) {
				post.mean[j] <- post.mean[j] + as.numeric(! v.gibbs)/(n.iter - n.burn)
			    post.mean.I[j] <- post.mean.I[j] + as.numeric(v.gibbs == 1)/(n.iter - n.burn)
			    post.mean.D[j] <- post.mean.D[j] + as.numeric(v.gibbs == 2)/(n.iter - n.burn)
			}		
		}
		##Update lambda##
		if (update.lambda) {
			lambda.new <- runif(1,0,1)
			log10bf.new <- BF.cgm(suff.stat, which(ind.infer == 1), which(ind.infer == 2), which(ind.infer == 0), n.ind, sigma.a, weights.sigma, lambda.new)
			if (log10bf.new - log10bf.use > log(runif(1), base=10)) {
				lambda.old <- lambda.new
			}
		}
		
		if (i > n.burn) {
			if (dirichlet) {
				gibbs.probs <- rdirichlet(1, theta + c(length(which(ind.infer == 0)), length(which(ind.infer == 1)), length(which(ind.infer == 2)) - 1) )            #Gibbs sampler for (p_U, p_I, p_D)
				post.prob.nodes <- post.prob.nodes + gibbs.probs/(n.iter - n.burn)	
			}			
			lambda.vec[i - n.burn] = lambda.old
		}
	}
	
	return(list(post.mean.U=post.mean, post.mean.I=post.mean.I, post.mean.D=post.mean.D, D.ind=D.gam, post.probs=post.prob.nodes, mean.lambda=lambda.vec))
}