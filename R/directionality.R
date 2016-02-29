
###Use a Gibbs sampler to estimate directionality of a graph; get posteriors of phenotypes lying in a particular category (I, U, D)

# estimate these parameters directly; this can be the prior
theta = c(0.4, 0.5, 0.1)			#What do sets of phenotypes do the indices correspond to? Assuming last one is Directly affected


# assumes genotypes: one per phenotype, same order, -1's if there's no eqtl.
# phenotypes: 
add_directionality <- function(ugraph, genotypes, phenotypes)		##Input is an undirected graph, genotypes, expression data; output is posterior probability of the arrow direction for each node; we can have bidirectional arrows if two genes are directly affected by a cis-eQTL and they directly affect one-another
{
   N <- length(V(ugraph))
   directions <- matrix(0, N, N)

   #V(ggmsg)$name # gene names

   #V(g) [ degree(g) > 5 ]

   iterations = 1000
   dirgibbs = matrix(0,N,N)
   for(iter in 1:iterations) {
    for(i in 1:N) {
     geno <- genotypes[,i]			#Genotypes of cis-eQTL for individuals
     if(geno[1] > -1) {					#Infer arrows if the SNP is an eQTL
        nei <- neighbors(ugraph,(i-1))			#Indices-1 of neighbors of gene i
	#print(nei+1)
	#print("running gibbs")
   	#print(i)
        directions <- gs_single_node(phenotypes, i, nei+1, directions, geno, theta)		#Draw arrows to and from neighbors of gene i
     }
    }
    dirgibbs = dirgibbs+directions
   }
   return(dirgibbs/iterations)
}

gs_single_node <- function(gex, index, nei, T, G, theta) {		#T is complete marix of current directions, G is genotype, theta is prior probability an arrow points one way or the other
   if(is.null(nei) | length(nei) == 0) {
      return(T)
   }
   N <- nrow(T)
   for(i in 1:length(nei)) {
      n = nei[i]
      if(T[index, n] == 1) {
         T[index,n] = 0      
      }
      # build set of affecteds
      afflist <- T[index,nei]*(nei)			#A 1 corresponds to an arrow pointing TOWARDS gene 'index'
      afflist <- afflist[afflist > 0]
      afflist <- c(afflist, n) #add on current gene as last
      aff <- cbind(gex[,index],gex[,afflist])			#Gene expression of gene index and genes that point TOWARDS index; the last column is the gene expression of the neighbor of interest
      znull <- rep(0, length(afflist)+1)
      bf0 <- MVassoctestfixedz(t(aff), G, c(rep(1,length(afflist)),0), znull)$lbf		#Neighbor i points away from index; this function is not in either of Barbara's other scripts
      bf1 <- MVassoctestfixedz(t(aff), G, c(rep(1,length(afflist)),1), znull)$lbf		#Neighbor i points towards index
      reversearrow <- T[n,index]
      bf0 <- (10^bf0)*theta[reversearrow+1]
      bf1 <- (10^bf1)*theta[reversearrow+2]
      pt <- bf0/(bf0+bf1)			#Posterior probability arrow points away from index
      print("Node")
      print(index)
      #print(nei)
      print(n)
      print("without direction")
      print(bf0)
      print("with direction")
      print(bf1)
      print(pt)
      # gibbs sampler
      u = rbinom(1,1,pt)
      print(u==1)
      if(u == 0) {
         T[index,n] = 1			#Record a 1 if neighbor i points towards index
      }      
   }
   return(T)
}

# assumes the first gene is the cis-eQTL of G
compare_trans_vs_null <- function(gex, G) 
{
   z <- c(1,1)
   znull <- c(1,0)
   return(MVassoctestfixedz(t(gex), G, z, znull)$lbf)
}
