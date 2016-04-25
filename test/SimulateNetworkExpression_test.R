####The purpose of this script is to simulate gene expression data given a gene g, a cis-eQTL s and neighbors of g, N_g
####Created 151130

###Input: A covariance matrix 'Sigma' for the expression of g and its neighbors, an effect size 'theta' for the cis-eQTL, a partition 'gamma' of genes in the network into Directly affected, Indirectly affected, and Unaffected and a genotype vector for 'n.ind' individuals at s

n.ind <- 1e4		#Number of independent individuals
n.nei <- 9			#Number of neighbors of g
f.s <- 0.1			#Minor Allele Frequency
Sigma <- create.sigma(n.nei+1)
lambda <- 1		##Shrinkage term; shrink mean of indirect effects
D.gam <- (4)
I.gam <- (1:3)			##Indices of indirectly affected neighbors; correspond to indices of Sigma
U.gam <- (5:(n.nei+1))		##Indices of unaffected neighbors; correspond to indices of Sigma
theta = 0.1
X.s <- rbinom(n=n.ind, size=1, prob=f.s) + rbinom(n=n.ind, size=1, prob=f.s)			#Genotype vector of n.ind individuals under HWE

sigma.a <- c(0.01, 0.4)				#Prior sd's on effect; used in BF calculation
weights.sigma <- c(1)					#Relative weights of each element of sigma.a
m = n.nei + 1					#degrees of freedom used in Wishart Prior; must be at least n.nei + 1

###Output is a matrix of log10 Bayes Factors for all 2^n.nei possible configurations



####Simulate expression data for n.ind individuals####

Y.gex <- Sim.gex(Sigma, D.gam, I.gam, U.gam, theta, X.s, lambda=lambda)			#n.ind x (n.nei + 1) gene expression matrix in the order of input Sigma

####Sufficient statistics necessary to calculate Bayes Factor####

suff.stat <- Suff.stat(Y.gex, X.s)

####Compute Bayes Factors####
#Return a 2^n.nei x (n.nei+2) matrix of log10 Bayes factors; the last column are the log10 BF's and the preceding columns give the partitions
#D.gam indicates the position of the direct effect

BF.obj <- Log10BF(suff.stat, D.gam, n.ind, sigma.a, weights.sigma, m)




create.sigma <- function(n) {			##Generates a random PD matrix 1/n*X'X, where X is an n x n random matrix
	X <- matrix(rnorm(n^2), n, n)
	return(1/n * t(X) %*% X)
}

#Simulate Expression for a single network, involving a gene g directly affected by cis-eQTL s, and gene g's neighbors whic have been partitioned into Indirectly affected genes (I.gam) and Unaffected genes (U.gam)
#theta is the cis-eQTL effect size and lambda is a shrinkage term to shrink indirect effects; in real data, indirect effects are generally not exactly as one would expect under statistical independence assumptions
#D.gam, I.gam, U.gam are the partition of the indices of Sigma into Directly affected, Indirectly affecte and Unaffected genes
#Output is a n.ind x (n.nei + 1) gene expression matrix
Sim.gex <- function(Sigma.in, D.gam, I.gam, U.gam, theta, X.s, lambda=1) {
	n.ind <- length(X.s)
	n.nei <- dim(Sigma.in)[1] - 1
	
	pos.vec <- c(I.gam, D.gam, U.gam)
	Sigma <- Sigma.in[pos.vec, pos.vec]			#Sigma = U Sigma.in U' for some permutation matrix U, given by I.gam, D.gam, U.gam
	
	perm.matrix <- array(0, dim=dim(Sigma.in))		#U in above equation
	for (r in 1:dim(Sigma.in)[1]) {
		perm.matrix[r, pos.vec[r]] = 1
	}
	
	I.gam.tmp <- (1:length(I.gam))
	D.gam.tmp <- length(I.gam) + 1
	U.gam.tmp <- ((D.gam.tmp+1):(length(I.gam) + length(U.gam) + 1))
	
	Sigma.12 <- Sigma[I.gam.tmp, c(D.gam.tmp, U.gam.tmp)]
	Sigma.22 <- Sigma[c(D.gam.tmp, U.gam.tmp), c(D.gam.tmp, U.gam.tmp)]
	
	tmp.vec <- lambda * Sigma.12 %*% chol2inv(chol(Sigma.22))[,1]
	mu.vec <- c(tmp.vec, 1, rep(0, length(U.gam)))
	
	R <- chol(Sigma)			#R'R = Sigma
	resids.iid <- matrix(rnorm((n.nei + 1)*n.ind), nrow=n.nei+1, ncol=n.ind)
	resids <- t(R) %*% resids.iid		#Matrix of gene expression residuals, a (n.nei + 1) x n matrix
	mu.mat <- cbind(mu.vec) %*% rbind(theta*X.s)		#Mean matrix of gene expression values, a (n.nei + 1) x n matrix
	
	Y.gex <- mu.mat + resids			#A (n.nei + 1) x n.ind expression matrix; The columns of Y are (Y_I, Y_D, Y_U)
	return(t(Y.gex) %*% perm.matrix)		#Rotate Y.gex into the order of Sigma.in; The output is a n.ind x (n.nei + 1) gene expression matrix
}


##Return sufficient statistics necessary to compute Bayes Factors in downstream analysis
##Y.gex is assumed to be a n.ind x (n.nei + 1) matrix, X.s a n.ind-vector of genetypes at s
##The output is a list, with names 'SYY', 'sxx', 'SYX', 'SY1', 'mu.g'
Suff.stat <- function(Y.gex, X.s) {
	n.ind <- length(X.s)
	SYY <- 1/n.ind * t(Y.gex) %*% Y.gex		#(n.nei + 1) x (n.nei + 1) matrix
	sxx <- 1/n.ind * sum(X.s*X.s)					#a scalar ~ 2f.s(1-f.s)
	SYX <- as.vector(1/n.ind * t(Y.gex) %*% cbind(X.s))		#a (n.nei + 1) vector
	SY1 <- apply(Y.gex, 2, mean)		#a (n.nei + 1) vector
	mu.g <- mean(X.s)			#a scalar
	
	return(list(SYY = SYY, sxx = sxx, SYX = SYX, SY1 = SY1, mu.g = mu.g))
}


##Compute log10 Bayes factors for all possible partitions gamma

Log10BF <- function(suff.stat, D.gam, n.ind, sigma.a, weights.sigma, m) {
	var.a <- sigma.a^2
	weights <- weights.sigma/sum(weights.sigma)
	
	SYY <- suff.stat$SYY
	sxx <- suff.stat$sxx
	SYX <- suff.stat$SYX
	SY1 <- suff.stat$SY1
	mu.g <- suff.stat$mu.g
	
	n.nei <- dim(SYY)[1] - 1	
	
	BF.mat <- array(1, dim=c(2^n.nei, n.nei + 2))			#Columns are c(c(Indicators if gene is affected by s), BF); assume that graph is such that there is only 1 possible direct effect, s -> g
	tmp.list <- list()
	for (i in 1:n.nei) {
		tmp.list[[i]] <- c(0,1)
	}
	BF.mat[,(1:(n.nei+1))[-D.gam]] <- as.matrix(expand.grid(tmp.list))		#BF.mat now contains all possible groupings, while keeping the direct effect fixed at 1
	

	for (i in 1:nrow(BF.mat)) {
		U.index <- which(BF.mat[i,] == 0)			#Indices of unaffected genes
		n.u <- length(U.index)
		
		##Build XX1 = 1/n.ind * [1 Y_u X_s]' [1 Y_u X_s] and XX0 = 1/n.ind * [1 Y_u]' [1 Y_u]##
		
		XX1 <- array(1, dim=c(n.u+2, n.u+2))
		
		if (n.u > 0) {
			XX1[1, 2:(n.u+2)] = c(SY1[U.index], mu.g)
			XX1[2:(n.u+2), 1] = c(SY1[U.index], mu.g)
			XX1[n.u+2, 2:(n.u+2)] = c(SYX[U.index], sxx)
			XX1[2:(n.u+1), n.u+2] = SYX[U.index]
			XX1[2:(n.u+1), 2:(n.u+1)] = SYY[U.index, U.index]
		
			XX0 <- XX1[1:(n.u+1), 1:(n.u+1)]	
		} else {
			XX1[1,2] <- mu.g
			XX1[2,1] <- mu.g
			XX1[2,2] <- sxx
			
			XX0 <- matrix(1, 1, 1)
		}
						
		##Build RSS1 = 1/n.ind * RSS(Y_D | [1 Y_u X_s], K) and RSS0 = 1/n.ind * RSS(Y_D | [1 Y_u], 0)##
		
		vec.i <- c(SY1[D.gam], SYY[D.gam, U.index], SYX[D.gam])
		
		if (n.u > 0) {
			RSS0 <- SYY[D.gam, D.gam] - sum( vec.i[-(n.u+2)] * solve(XX0, vec.i[-(n.u+2)]) )
		} else {
			RSS0 <- SYY[D.gam, D.gam] - vec.i[1]*vec.i[1]
		}		
		
		log.bf.i <- rep(0, length=length(var.a))
		for (s in 1:length(var.a)) {
			tmp.mat.i <- XX1
			tmp.mat.i[n.u+2, n.u+2] <- tmp.mat.i[n.u+2, n.u+2] + 1/var.a[s]/n.ind		#tmp.mat.i = XX1 + 1/n.ind * K
			RSS1 <- SYY[D.gam, D.gam] - sum(vec.i * solve(tmp.mat.i, vec.i))
			log.bf.i[s] <- length(D.gam)/2 * ( -log(var.a[s]) + log(det(XX0)) - log(n.ind) - log(det(tmp.mat.i)) ) + (n.ind + m - (n.nei + 1) + n.u + length(D.gam))/2 * ( log(RSS0) - log(RSS1) )
		}
		
		##Compute log Bayes Factor##
		
		BF.mat[i, n.nei + 2] = log( sum(weights*exp(log.bf.i - max(log.bf.i))) ) + max(log.bf.i)
	}
	BF.mat[,n.nei + 2] <- BF.mat[,n.nei + 2]/log(10)		#Log10 BF
	
	return(list(BF.mat = BF.mat, D.ind = D.gam))
}


#Compute log10 Bayes factor for a single partition gamma
#Inputs are the sufficient statistics, the indices corresponding to the directly affected and unaffected genes, #independent samples, vector of sigma.a's, weights for sigma.a's, m in Wisharg prio
#Output is the log10BF and the indices corresponding to the partitions

Log10BF.sing <- function(suff.stat, D.gam, U.gam, n.ind, sigma.a, weights.sigma, m) {
	var.a <- sigma.a^2
	weights <- weights.sigma/sum(weights.sigma)
	
	SYY <- suff.stat$SYY
	sxx <- suff.stat$sxx
	SYX <- suff.stat$SYX
	SY1 <- suff.stat$SY1
	mu.g <- suff.stat$mu.g
	
	n.nei <- dim(SYY)[1] - 1	
	n.u <- length(U.gam)
		
	##Build XX1 = 1/n.ind * [1 Y_u X_s]' [1 Y_u X_s] and XX0 = 1/n.ind * [1 Y_u]' [1 Y_u]##
		
	XX1 <- array(1, dim=c(n.u+2, n.u+2))
		
	if (n.u > 0) {
		XX1[1, 2:(n.u+2)] = c(SY1[U.gam], mu.g)
		XX1[2:(n.u+2), 1] = c(SY1[U.gam], mu.g)
		XX1[n.u+2, 2:(n.u+2)] = c(SYX[U.gam], sxx)
		XX1[2:(n.u+1), n.u+2] = SYX[U.gam]
		XX1[2:(n.u+1), 2:(n.u+1)] = SYY[U.gam, U.gam]
	
		XX0 <- XX1[1:(n.u+1), 1:(n.u+1)]	
	} else {
		XX1[1,2] <- mu.g
		XX1[2,1] <- mu.g
		XX1[2,2] <- sxx
		
		XX0 <- matrix(1, 1, 1)
	}
					
	##Build RSS1 = 1/n.ind * RSS(Y_D | [1 Y_u X_s], K) and RSS0 = 1/n.ind * RSS(Y_D | [1 Y_u], 0)##
	
	vec <- c(SY1[D.gam], SYY[D.gam, U.gam], SYX[D.gam])
	
	if (n.u > 0) {
		RSS0 <- SYY[D.gam, D.gam] - sum( vec[-(n.u+2)] * solve(XX0, vec[-(n.u+2)]) )
	} else {
		RSS0 <- SYY[D.gam, D.gam] - vec[1]*vec[1]
	}		
	
	log.bf <- rep(0, length=length(var.a))
	for (s in 1:length(var.a)) {
		tmp.mat <- XX1
		tmp.mat[n.u+2, n.u+2] <- tmp.mat[n.u+2, n.u+2] + 1/var.a[s]/n.ind		#tmp.mat = XX1 + 1/n.ind * K
		RSS1 <- SYY[D.gam, D.gam] - sum(vec * solve(tmp.mat, vec))
		log.bf[s] <- length(D.gam)/2 * ( -log(var.a[s]) + log(det(XX0)) - log(n.ind) - log(det(tmp.mat)) ) + (n.ind + m - (n.nei + 1) + n.u + length(D.gam))/2 * ( log(RSS0) - log(RSS1) )
	}
	
	##Compute log Bayes Factor##
	
	logBF = log( sum(weights*exp(log.bf - max(log.bf))) ) + max(log.bf)

	log10BF <- logBF/log(10)		#Log10 BF
	
	return(list(log10BF = log10BF, D.ind = D.gam, I.ind = (1:(n.nei+1))[ -c(D.gam,U.gam) ], U.ind=U.gam))
}

