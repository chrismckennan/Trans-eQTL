
library("GeneNet")


# shrinkage estimation of GGM
build_shrinkage_ggm <- function(gex, labels)
{
  nodelables = labels

  inferred.pcor <- ggm.estimate.pcor(gex) # assumes static gex data, rows are multivariate observation
  test.results <- ggm.test.edges(inferred.pcor, fdr=T, direct=F, plot=T)

  pruned.results <- extract.network(test.results, method.ggm=c("prob"), cutoff.ggm=0.8, verbose=TRUE)

  gr <- ggm.make.igraph(pruned.results, nodelables)
  #show.edge.weights(gr)
  return(gr)
}

plot_ggm <- function(gr)
{
  library("Rgraphviz")
  plot(gr)

}

network_features <- function(gn)
{
   ew <- edgeWeights(gn)
   # degree
   degree <- apply(ew, 1, function(x) { return(sum(abs(x)>0)) })
   print(sum(degree == 0))
   hist(degree, nclass=100)
}

# traditional relevance	network	(Butte et al., 2000)
# GEX data is n		rows (samples) by p columns (genes)

build_correlation_graph <- function(gex, gexnames)			#Edge weights are correlation p-value, threshholded at 2*0.95/(p^2 - p)
{
   n = nrow(gex)
   p = ncol(gex)
   corgraph = matrix(0,	p, p)
   edges = c()
   for(i in 1:(p-1)) {
     for(j in (i+1):p) {
      	   c2g = (cor.test(gex[,i], gex[,j])$p.value)
	   corgraph[i,j] <- c2g
	   corgraph[j,i] <- c2g
	   if(c2g < (0.95/(((p*p)-p)/2))) {
 	   	  edges = c(edges, (i-1), (j-1))
	   }
     }
   }
   g = graph(edges, n=p, directed=F)
   V(g)$color <- 'black'
   V(g)$label <- gexnames
   V(g)$name <- gexnames
   print(0.95/(((p*p)-p)/2))
   return(list(corpval = corgraph, g=g))
}

build_correlation_graph_fdr <- function(gex, gexnames)
{
   p = nrow(gex)
   n = ncol(gex)
   betahat = cor(t(gex))
   se.betahat = sqrt((1-betahat^2)/n)
   pval = 2*pnorm(-abs((betahat/se.betahat)[se.betahat>0]))
   pval = cbind(0,t(matrix((pval),p-1,p)))
   for(i in 2:p) {
      pval[i, 1:(i-1)] = pval[i, 2:i]
      pval[i,i] = 0
   } 
   return(pval)
}

compute_fdr <- function(pval_mat, cutoff) {
    p = nrow(pval_mat)
    total <- sum(apply(pval_mat, 1, function(x) { return(sum(x<cutoff)) } ))
    return(total/(p*cutoff))
}

# assumes n x p matrix of gex data
build_correlation_graph_fdr <- function(gex, gexnames)
{
   p = nrow(gex)
   n = ncol(gex)
   betahat = cor(t(gex))
   se.betahat = sqrt((1-betahat^2)/n)
   pval = 2*pnorm(-abs((betahat/se.betahat)[se.betahat>0]))
   
   qval <- qvalue_approx(as.vector(pval))
   qval <- cbind(0, matrix(qval<0.05,n,n-1))
   print(dim(qval))
   lapply(1:n, function(x) { t = qval[1,x]
                             qval[1,x] = qval[x,x]
                             qval[x,x] = t }) 
   return(list(edges=qval, corg=betahat))
}


# http://genomics.princeton.edu/storeylab/qvalue/results.R
qvalue <- function(p, alpha=NULL, lam=NULL, robust=F)
{
#This is a function for estimating the q-values for a given set of p-values. The
#methodology comes from a series of recent papers on false discovery rates by John
#D. Storey et al. See http://www.stat.berkeley.edu/~storey/ for references to these
#papers. This function was written by John D. Storey. Copyright 2002 by John D. Storey.
#All rights are reserved and no responsibility is assumed for mistakes in or caused by
#the program.
#
#Input
#=============================================================================
#p: a vector of p-values (only necessary input)
#alpha: a level at which to control the FDR (optional)
#lam: the value of the tuning parameter to estimate pi0 (optional)
#robust: an indicator of whether it is desired to make the estimate more robust 
#        for small p-values (optional)
#
#Output
#=============================================================================
#remarks: tells the user what options were used, and gives any relevant warnings
#pi0: an estimate of the proportion of null p-values
#qvalues: a vector of the estimated q-values (the main quantity of interest)
#pvalues: a vector of the original p-values
#significant: if alpha is specified, and indicator of whether the q-value fell below alpha 
#    (taking all such q-values to be significant controls FDR at level alpha)

#This is just some pre-processing
    m <- length(p)
#These next few functions are the various ways to automatically choose lam
#and estimate pi0
    if(!is.null(lam)) {
        pi0 <- mean(p>lam)/(1-lam)
        remark <- "The user prespecified lam in the calculation of pi0."
    }
    else{
        remark <- "A smoothing method was used in the calculation of pi0."
        library(stats)
        lam <- seq(0,0.95,0.01)
        pi0 <- rep(0,length(lam))
        for(i in 1:length(lam)) {
        pi0[i] <- mean(p>lam[i])/(1-lam[i])
        }
        spi0 <- smooth.spline(lam,pi0,df=3,w=(1-lam))
        pi0 <- predict(spi0,x=0.95)$y
    }
#The q-values are actually calculated here
    u <- order(p)
    v <- rank(p)
    qvalue <- pi0*m*p/v
    if(robust) {
        qvalue <- pi0*m*p/(v*(1-(1-p)^m))
        remark <- c(remark, "The robust version of the q-value was calculated. See Storey JD (2002) JRSS-B 64: 479-498.")
    }
    qvalue[u[m]] <- min(qvalue[u[m]],1)
    for(i in (m-1):1) {
    qvalue[u[i]] <- min(qvalue[u[i]],qvalue[u[i+1]],1)
    }
    #Here the results are returned
    if(!is.null(alpha)) {
        return(list(remarks=remark, pi0=pi0, qvalue=qvalue, significant=(qvalue <= alpha), pvalue=p))
    }
    else {
        return(list(remarks=remark, pi0=pi0, qvalue=qvalue, pvalue=p))
    }
}


# http://genomics.princeton.edu/storeylab/qvalue/results.R
#This is a function for estimating the q-values for a given set of p-values. The
#methodology comes from a series of recent papers on false discovery rates by John
#D. Storey et al. See http://www.stat.berkeley.edu/~storey/ for references to these
#papers. This function was written by John D. Storey. This is an approxate version.
#
#Input
#=============================================================================
#p: a vector of p-values (only necessary input)
#alpha: a level at which to control the FDR (optional)
#lam: the value of the tuning parameter to estimate pi0 (optional)
#robust: an indicator of whether it is desired to make the estimate more robust 
#        for small p-values (optional)
#
#Output
#=============================================================================
#remarks: tells the user what options were used, and gives any relevant warnings
#pi0: an estimate of the proportion of null p-values
#qvalues: a vector of the estimated q-values (the main quantity of interest)
#pvalues: a vector of the original p-values
#significant: if alpha is specified, and indicator of whether the q-value fell below alpha 
#    (taking all such q-values to be significant controls FDR at level alpha)

# assumes pi0=1
# does not take into account ties. use unique or p<=x for that (but will take longer).
qvalue_approx <- function(p)
{
    m <- length(p)
    u <- sort(p,index.return=TRUE)
    qval <- lapply(1:m, function(x) { return((u$x[x]*m)/x) } )
    return(qval[u$ix])
}

#sigma hat/root n
#sigma hat residual variance: 1 - beta hat squared
#beta hat/sigma hat/root n

#x = matrix(rnorm(48000),nrow=480, ncol=5)
#betahatx = cor(x)
#se.betahatx = sqrt((1-betahatx^2)/480)
#pval = (betahatx/(se.betahatx/sqrt(480)))


#se.betahatx = betahatx*sqrt(478/(1-betahatx^2))
#diag(se.betahatx) = .0
#qqnorm((betahatx/se.betahatx)[se.betahatx>0])
#qval = qvalues(2*pnorm(-abs((betahatx/se.betahatx)[se.betahatx>0])))

#versus

#cgfake <- build_correlation_graph(x, 1:100)

build_correlation_graph_differential_2not1 <- function(gex1, gex2, gexnames)
{
   p = nrow(gex1)
   n = ncol(gex1)
   gex = c()
   for(i in 1:p) {
      gex = cbind(gex, qqnorm(lm(gex2[i,]~gex1[i,])$residual,plot.it=FALSE)$x)		#Normal quantiles of regression residuals
   }
   return(build_correlation_graph_fdr(gex, gexnames))
}

build_correlation_graph_differential_1not2 <- function(gex1, gex2, gexnames)
{
   p = nrow(gex1)
   n = ncol(gex1)
   gex = c()
   for(i in 1:p) {
      gex = cbind(gex, qqnorm(lm(gex2[i,]~gex1[i,])$residual,plot.it=FALSE)$x)
   }
   return(build_correlation_graph_fdr(gex, gexnames))
}

build_correlation_graph_differential_different <- function(gex1, gex2, gexnames)
{
   p = nrow(gex1)
   n = ncol(gex1)
   gex = c()
   for(i in 1:p) {
      gex = cbind(gex, qqnorm(lm((gex2[i,]-gex1[i,])~(gex2[i,]+gex1[i,]))$residual,plot.it=FALSE)$x)
   }
   return(build_correlation_graph_fdr(gex, gexnames))
}

build_correlation_graph_differential_opposite <- function(gex1, gex2, gexnames)
{
   p = nrow(gex1)
   n = ncol(gex1)
   gex = c()
   for(i in 1:p) {
      gex = cbind(gex, qqnorm(lm((gex2[i,]-gex1[i,])~(gex2[i,]+gex1[i,]))$residual,plot.it=FALSE)$x)
   }
   return(build_correlation_graph_fdr(gex, gexnames))
}