
# estimate these parameters directly; this can be the prior
theta = c(0.4, 0.5, 0.1)


# assumes genotypes: one per phenotype, same order, -1's if there's no eqtl.
# phenotypes: 
add_directionality <- function(ugraph, genotypes, phenotypes)
{
   N <- length(V(ugraph))
   directions <- matrix(0, N, N)

   #V(ggmsg)$name # gene names

   #V(g) [ degree(g) > 5 ]

   iterations = 1000
   dirgibbs = matrix(0,N,N)
   for(iter in 1:iterations) {
    for(i in 1:N) {
     geno <- genotypes[,i]
     if(geno[1] > -1) {
        nei <- neighbors(ugraph,(i-1))
	#print(nei+1)
	#print("running gibbs")
   	#print(i)
        directions <- gs_single_node(phenotypes, i, nei+1, directions, geno, theta)
     }
    }
    dirgibbs = dirgibbs+directions
   }
   return(dirgibbs/iterations)
}

gs_single_node <- function(gex, index, nei, T, G, theta) {
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
      afflist <- T[index,nei]*(nei)
      afflist <- afflist[afflist > 0]
      afflist <- c(afflist, n) #add on current gene as last
      aff <- cbind(gex[,index],gex[,afflist])
      znull <- rep(0, length(afflist)+1)
      bf0 <- MVassoctestfixedz(t(aff), G, c(rep(1,length(afflist)),0), znull)$lbf
      bf1 <- MVassoctestfixedz(t(aff), G, c(rep(1,length(afflist)),1), znull)$lbf
      reversearrow <- T[n,index]
      bf0 <- (10^bf0)*theta[reversearrow+1]
      bf1 <- (10^bf1)*theta[reversearrow+2]
      pt <- bf0/(bf0+bf1)
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
         T[index,n] = 1
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
