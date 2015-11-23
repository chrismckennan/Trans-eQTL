#simple bf calculation for single SNP
log10BF = function(g,y,sigmaa,sigmad){
subset = complete.cases(y) & complete.cases(g)
y=y[subset]
g=g[subset]
n=length(g)
X = cbind(rep(1,n),g,g==1)
invnu = diag(c(0,1/sigmaa^2,1/sigmad^2))
invOmega = invnu + t(X) %*% X
B = solve(invOmega, t(X) %*% cbind(y))
invOmega0 = n
return(-0.5*log10(det(invOmega)) + 0.5*log10(invOmega0) - log10(sigmaa) - log10(sigmad) -(n/2) * (log10( t(y- X %*% B) %*% y) - log10(t(y) %*% y - n*mean(y)^2) ))  
}

# compute logZnd ratio in (52) of Minka
logZndratio=function(n0,d,n){
  return(sum(lgamma(0.5*(n+n0+1-(1:d)))-lgamma(0.5*(n0+1-(1:d)))))
}



#compute log marginal likelihood for multivariate regression
#from Minka paper eqn (52),
LogMarginalLikwithconst=function(Y,X,K,N0,S0){
  d = dim(Y)[1]
  N = dim(Y)[2]
  SXX = X %*% t(X) + K
  SYY = Y %*% t(Y)
  SYX = Y %*% t(X)
  SYgX = SYY - SYX %*% chol2inv(chol(SXX)) %*% t(SYX)
  return(logZndratio(N0,d,N) - (N*d/2)*log(pi) +  (d/2) * (determinant(K)$modulus - determinant(SXX)$modulus) + (N0/2) * determinant(S0)$modulus - (N+N0)/2 * determinant(SYgX + S0)$modulus)
}

#log like in limit S0 = 0, wiht no constants, and K=diag(0,0,...,sigmaa)
LogMarginalLiknoconst=function(Y,X,sigmaa,N0){
  d = dim(Y)[1]
  N = dim(Y)[2]
  
  SXX = X %*% t(X)
  SXX[2,2] = SXX[2,2]+sigmaa^{-2}
  SYY = Y %*% t(Y)
  SYX = Y %*% t(X)
  LXX = chol(SXX)
  SYgX = SYY - SYX %*% chol2inv(LXX) %*% t(SYX)
  return((d/2) * (-2*log(sigmaa) - 2*sum(log(diag(LXX)))) - (N+N0)/2 * determinant(SYgX)$modulus)
}

# Compute Bayes Factor for Y a d by N matrix, and X a N-vector of genotypes
# Added a parameter epsilon, with Psi = epsilon I, with the idea
# to check that if Y has exact linear dependencies then it does not effect the
# limiting BF
logBFall = function(g,Y,sigmaa,N0=0,epsilon=0){
subset = complete.cases(t(Y)) & complete.cases(g)
Y=Y[,subset]
g=g[subset]
d = dim(Y)[1]
N = dim(Y)[2]
if(N0==0){N0 = d-1}
X = rbind(rep(1,N),g)
SXX = X %*% t(X)
  SXX[2,2] = SXX[2,2]+sigmaa^{-2}
  SYY = Y %*% t(Y)
  centerY = Y - apply(Y,1,mean)
  SYY0 = centerY %*% t(centerY)
  SYX = Y %*% t(X)
  LXX = chol(SXX)
  SYgX = SYY - SYX %*% chol2inv(LXX) %*% t(SYX)
return((d/2) * (-2*log(sigmaa) - 2*sum(log(diag(LXX)))+ log(N)) + (N+N0)/2 * (determinant(SYY0+diag(epsilon,d))$modulus - determinant(SYgX+diag(epsilon,d))$modulus)[1])
}


# Compute smoothed Bayes Factor for Y a d by N matrix, and X a N-vector of genotypes
# Added a parameter epsilon, with Psi = epsilon I, with the idea
# to check that if Y has exact linear dependencies then it does not effect the
# limiting BF
logBFall.smoothed = function(g,Y,sigmaa){
subset = complete.cases(t(Y)) & complete.cases(g)
Y=Y[,subset]
g=g[subset]
d = dim(Y)[1]
N = dim(Y)[2]
N0 = d #N0 = m, set to >= d for smoothing
epsilon = N0
X = rbind(rep(1,N),g)
SXX = X %*% t(X)
  SXX[2,2] = SXX[2,2]+sigmaa^{-2}
  SYY = Y %*% t(Y)
  centerY = Y - apply(Y,1,mean)
  SYY0 = centerY %*% t(centerY)
  SYX = Y %*% t(X)
  LXX = chol(SXX)
  SYgX = SYY - SYX %*% chol2inv(LXX) %*% t(SYX)
return((d/2) * (-2*log(sigmaa) - 2*sum(log(diag(LXX)))+ log(N)) + (N+N0)/2 * (determinant(SYY0+diag(epsilon,d))$modulus - determinant(SYgX+diag(epsilon,d))$modulus)[1])
}


#this one, March 2010, written to test the new "efficient" BF calculation
#based on rank one update
#note that now Y is n by d
logBFall.rankone = function(g,Y,X0,sigmaa,N0=0){
subset = complete.cases(Y) & complete.cases(g)
Y=Y[subset,]
g=g[subset]
d = dim(Y)[2]
N = dim(Y)[1]
X1 = cbind(g)
if(N0==0){N0 = d-1}
X = cbind(X0,X1)
SXX0 = t(X0) %*% X0
invSXX0 = chol2inv(chol(SXX0))
phi1 = t(Y) %*% X0
phi2 = t(X0) %*% X1
C = t(X1) %*% X1 - t(phi2) %*% invSXX0 %*% phi2
invC = chol2inv(chol(C))
RSS0 = t(Y) %*% Y - phi1 %*% invSXX0 %*% t(phi1)
invRSS0 = chol2inv(chol(RSS0))
v = t(Y) %*% X1 - phi1 %*% invSXX0 %*% phi2
lambda = sigmaa^{-2} / (t(X1) %*% X1 - t(phi2) %*% invSXX0 %*% phi2)
k = as.numeric(1/(1 +lambda))
return((d/2) * log(1-k) - 0.5*(N+N0)*log(1-k*invC %*% t(v) %*% invRSS0 %*% v))
}

logBFall.rankone.new = function(g,Y,X0,sigmaa,m=0){
subset = complete.cases(Y) & complete.cases(g)
Y=Y[subset,]
g=g[subset]
d = dim(Y)[2]
n = dim(Y)[1]
X1 = cbind(g)
if(m==0){m = d-1}
X = cbind(X0,X1)
phi0 = t(Y) %*% X0
phi1 = t(Y) %*% X1

S00 = t(X0) %*% X0
S01 = t(X0) %*% X1
S11 = t(X1) %*% X1

invS00 = chol2inv(chol(S00))
RSS0 = t(Y) %*% Y - phi0 %*% invS00 %*% t(phi0)
L = chol(RSS0)

C= S11 - t(S01) %*% invS00 %*% S01
v = phi1 - phi0 %*% invS00 %*% S01
lambda = sigmaa^{-2} / C
k = as.numeric(1/(1+lambda))

a = backsolve(L,v)

return((d/2) * log(1-k) - 0.5*(n+m)*log(1-(k/C) * sum(a*a)))
}

#now G is an n by p matrix of SNPs
logBFall.rankone.matrix = function(G,Y,sigmaa,pi0=0.5,m=0){
subset = complete.cases(Y) & complete.cases(G)
Y=Y[subset,]
G=G[subset,]
n = dim(Y)[1]
d = dim(Y)[2]
p = dim(G)[2] #number of SNPs
if(m==0){m = d-1}

Y =scale(Y,center=T,scale=F)  #center Y and G to avoid needing intercept in regression
G = scale(G,center=T,scale=F)

phi = crossprod(Y,G) # this is t(Y) %*% G, a d by p matrix
SYY = crossprod(Y) # t(Y) %*% Y, a d by d matrix
S11 = colSums(G*G) # a p vector of ||g|| values

prior = rep(0,3^d)
z=matrix(0,nrow=3^d,ncol=d)
lbf = matrix(0,nrow=3^d, ncol=p)

for(i in 0:(3^d-1)){
  for(j in 1:d){
    z[i+1,j]= (i %% 3^j) %/% 3^{j-1}
  }
  prior[i+1] = computeprior(z[i+1,],pi0)
  #print(c(z[i+1,],prior[i+1]))
  U = (z[i+1,]==0)
  D = (z[i+1,]==1)
  if(prior[i+1]>0){
  lbf[i+1,] = logBF.fromsummaries(phi,SYY,S11,U,D,n,m,d,sigmaa) #note we just don't bother computing for models with prior = 0
  } else {lbf[i+1,] = 0}
}
prior[1] = pi0

return(list(z=z,lbf=lbf,prior=prior))

}

#note the "drop=F" commands below stop R collapsing matrices into vectors inappropriately
logBF.fromsummaries = function(phi,SYY,S11,U,D,n,m,d,sigmaa){
dd = sum(D)
du= sum(U)

if(du>0){
LUU = chol(SYY[U,U,drop=F]) # a du by du matrix
phi0 = SYY[U,D,drop=F]      #a du by da matrix
c = cbind(backsolve(LUU,phi[U,,drop=F]))#c solves LUU c = phiU, c is a du by p matrix
b = cbind(backsolve(LUU, phi0))  # b is du by da, and solves LUU b = phi0, so b'b = phi0 SYYU^-1 phi0
}
else{c=matrix(0,nrow=1,ncol=p); b=matrix(0,nrow=1,ncol=da);}

C = S11 - colSums(c*c)
u = phi[D,,drop=F] - crossprod(b,c)
RSS0 = SYY[D,D,drop=F] - crossprod(b)
L0 = chol(RSS0)
a = backsolve(L0,u)

lambda = sigmaa^{-2} / C
k = as.numeric(1/(1+lambda))

return((dd/2) * log(1-k) - 0.5*(n+m-(d-sum(D)-sum(U)))*log(1-(k/C) * colSums(a*a)))
}


#simple univariate bf calculation for single SNP, for comparison
log10BF = function(g,y,sigmaa,sigmad){
subset = complete.cases(y) & complete.cases(g)
y=y[subset]
g=g[subset]
n=length(g)
X = cbind(rep(1,n),g,g==1)
invnu = diag(c(0,1/sigmaa^2,1/sigmad^2))
invOmega = invnu + t(X) %*% X
B = solve(invOmega, t(X) %*% cbind(y))
invOmega0 = n
return(-0.5*log10(det(invOmega)) + 0.5*log10(invOmega0) - log10(sigmaa) - log10(sigmad) -(n/2) * (log10( t(y- X %*% B) %*% y) - log10(t(y) %*% y - n*mean(y)^2) ))  
}



#return log lik for model where Y[aset,]is directly affected by X; dset is downstream; uset is upstream
AssocLogLik=function(Y,X,uset,aset,dset,sigmaa,n0,epsilon=1e-6)
{
complete=complete.cases(cbind(t(t(X)),t(Y)))
Y=Y[,complete]
X=X[complete]
d =nrow(Y)
n=ncol(Y)
S0=diag(epsilon,d)
onevec=matrix(rep(1,n),nrow=1)
sigma0sq=100
Yu=matrix(Y[uset,],ncol=n)
Ya=matrix(Y[aset,],ncol=n)
Yd=matrix(Y[dset,],ncol=n)
du=length(uset)
da=length(aset)
dd=length(dset)
loglik=0
if(du>0){
S0u = matrix(S0[uset,uset],nrow=du)
loglik =loglik+ LogMarginalLikwithconst(Yu,onevec,matrix(1/sigma0sq,nrow=1,ncol=1),n0-da-dd,S0u)
}
if(da>0){
S0a = matrix(S0[aset,aset],nrow=da)
loglik =loglik+LogMarginalLikwithconst(Ya,rbind(onevec,Yu,X),diag(c(1/sigma0sq,diag(as.matrix(S0[uset,uset])),1/sigmaa^2)),n0-dd,S0a)
}
if(dd>0){
S0d = matrix(S0[dset,dset],nrow=dd)
loglik=loglik+LogMarginalLikwithconst(Yd,rbind(onevec,Yu,Ya),diag(c(1/sigma0sq,diag(as.matrix(S0[c(uset,aset),c(uset,aset)]))),nrow=d-dd+1),n0,S0d)
}
return(loglik)
}

allones = function(x){return(prod(x==1))}

#note on prior
# Let di be the number of the d variables in each of the groups i=0,1,2
# we assume that the null model (all 0) has probability pi
# the sum of all alternatives (d1>=1) has probability 1-pi
# conditional on d1>=1, we assume that A=d1+d2 is uniform on 1...d,
# and conditional on A, d1 is uniform on 1 to A.
#conditional on d0, d1 and d2 we assume all d!/(d0! d1! d2!) configurations are equally likely
# thus Pr(given configuration) = pi if configuration is null,
# and (1-pi) (1/d) (1/(d1+d2)) d0! d1! d2!/d!
computeprior = function(z,pi){
  dvec = tabulate(z+1,nbin=3)
  d = length(z)
  return(ifelse(dvec[2]>0,(1-pi)*(1/d) * (1/(dvec[2]+dvec[3])) * 1/choose(d, dvec[1]) * 1/choose(d-dvec[1],dvec[2]), 0))
  
}

#only difference here is that we allow that conditional on A, d1 is uniform
# on 0 to A instead of 1 to A. This way, all might be indirect effects.
#reason to allow this is that this function is used for marginal BFs when
#the variable whose marginal value is being computed may be the direct effect
#and want to allow that others may be indirect
computepriorvec = function(z,pi){
  dvec = t(apply(z+1,1,tabulate,nbin=3))
  d = dim(z)[2]
  p = (ifelse(dvec[,1]<d,(1-pi)*(1/d) * (1/(dvec[,2]+dvec[,3]+1)) * 1/choose(d, dvec[,1]) * 1/choose(d-dvec[,1],dvec[,2]), 0))
  nulls = (apply(z,1,sum)==0)
  p[nulls]=pi
#  return(rep(1/dim(z)[1],dim(z)[1]))
  return(p)
}


MVassoctest=function(Y,X,sigmaa=0.4,n0=0,pi=0.5){
complete=complete.cases(cbind(t(t(X)),t(Y)))
d = nrow(Y)
Y=matrix(Y[,complete],nrow=d)
Y=Y/apply(Y,1,sd)
Y=matrix(Y,nrow=d)
X=X[complete]

prior = rep(0,3^d)
if(n0==0){n0=d+1}
z=matrix(0,nrow=3^d,ncol=d)
al=rep(0,3^d)
for(i in 0:(3^d-1)){
  for(j in 1:d){
    z[i+1,j]= (i %% 3^j) %/% 3^{j-1}
  }
  prior[i+1] = computeprior(z[i+1,],pi)
  #print(c(z[i+1,],prior[i+1]))
  
  al[i+1] = ifelse(prior[i+1]>0,AssocLogLik(Y,X,(1:d)[z[i+1,]==0],(1:d)[z[i+1,]==1],(1:d)[z[i+1,]==2],sigmaa,n0),0) #note we just don't bother computing for models with prior = 0
  
}

prior[1] = pi # set prior on null model to be pi
al[1] = AssocLogLik(Y,X,1:d,NULL,NULL,sigmaa,n0) # compute assocloglik for null
al = ifelse(prior==0,al[1],al) # set all the ones with prior = 0 to have same loglik as null

BF = exp(al-al[1])
lBF=log10(BF)

alldirect = apply(z,1,allones)
lbfall = log10(BF[alldirect==1])

posterior = prior*BF/sum(prior*BF)

#compute marginal bf for each dimension, being direct (lmarginalbf1) or indirect (lmarginalbf2) effect, treating inclusion as indep of other variables
lmarginalbf1 = rep(0,d)
lmarginalbf2 = rep(0,d)
for(i in 1:d){
  condprior = computepriorvec(as.matrix(z[,-i]),pi)
  lmarginalbf1[i] = log10( sum((condprior*BF)[z[,i]==1]) / sum(condprior[z[,i]==1] )) - log10( sum((condprior*BF)[z[,i]==0]) / sum(condprior[z[,i]==0]))
  lmarginalbf2[i]  = log10( sum((condprior*BF)[z[,i]==2]) / sum(condprior[z[,i]==2] )) - log10( sum((condprior*BF)[z[,i]==0]) / sum(condprior[z[,i]==0]))
}

lunivariatebf = rev(lBF[apply(z,1,sum)==(2*d - 1)])

loverallbf = log10(sum(prior[-1]*BF[-1])/sum(prior[-1]))
p0=t(z==0) %*% posterior
p1=t(z==1) %*% posterior
p2=t(z==2) %*% posterior

 return(list(z=z,prior=prior,loglik=al,p=cbind(p0,p1,p2),lbf=lBF,loverallbf = loverallbf, lmarginalbf1= lmarginalbf1,lmarginalbf2=lmarginalbf2,lunivariatebf = lunivariatebf, lbfall=lbfall))
}

fastMVassoctest=function(Y,X,sigmaa=0.4,n0=0,pi=0.5){
complete=complete.cases(cbind(t(t(X)),t(Y)))
d = nrow(Y)
Y=matrix(Y[,complete],nrow=d)
Y=Y/apply(Y,1,sd)
Y=matrix(Y,nrow=d)
X=X[complete]
n = ncol(Y)

prior = rep(0,3^d)
if(n0==0){n0=d+1}

z=matrix(0,nrow=3^d,ncol=d)
al=rep(0,3^d)
unaff = rep(0,d) #indicators for whether each coordinate is unaffected
aff = rep(0,d)
iter = 0
for(i in 0:(2^d-1)){
  for(j in 1:d){
    unaff[j] = ((i %% 2^j) %/% 2^{j-1})
  }
  unaff = as.logical(unaff)
  aff.set = (1:d)[!unaff]
  naff = length(aff.set)

  #Form X1 and then Omega1 = X t(X) + K
  X1 = rbind(rep(1,n),Y[unaff,],X)
  Omega1 = X1 %*% t(X1)
  dd = nrow(Omega1)
  Omega1[dd,dd] = Omega1[dd,dd] + sigmaa^(-2)
  
  Omega0 = Omega1[-dd,-dd]

  L1 = chol(Omega1)
  L0 = chol(Omega0) #note could be more efficient here by exploiting nestd structure of Omega0 vs Omega1
  ldet1 = 2*sum(log(diag(L1)))  # log determinants of Omega1
  ldet0 = 2*sum(log(diag(L0))) # and Omega0

  Omega1inv = chol2inv(L1)
  Omega0inv = chol2inv(L0)
  
  Ya = matrix(Y[!unaff,],ncol=n)
  SYY = Ya %*% t(Ya)
  SYX1 = Ya %*% t(X1)
  SYX0 = matrix(SYX1[,-dd],ncol=dd-1)
  
  Sigma1 = SYY - SYX1 %*% Omega1inv %*% t(SYX1)
  Sigma0 = SYY - SYX0 %*% Omega0inv %*% t(SYX0) 

  loglikbase = -log(sigmaa) + (1/2) * ldet0 - (1/2) * ldet1
  
  d.aff = rep(0,naff) #indicator for which are directly affected
  for(k in 0:(2^naff-1)){
    iter = iter+1
    for(l in 1:naff)
      d.aff[l] = (k %% 2^l) %/% 2^{l-1}
    d.aff = as.logical(d.aff)
    n.indirect.aff = sum(!d.aff)
    m = n0 - n.indirect.aff
    z[iter,aff.set] = 2-d.aff
    prior[iter] = computeprior(z[iter,],pi)
    if(sum(d.aff)>0){
      al[iter] = sum(d.aff)*loglikbase + ((n+m)/2)*determinant(as.matrix(Sigma0[d.aff,d.aff]),logarithm=T)$modulus -((n+m)/2)*determinant(as.matrix(Sigma1[d.aff,d.aff]),logarithm=T)$modulus} else
      {al[iter] = 0}
   
  }  
}


prior[length(prior)] = pi # set prior on null model to be pi
BF = exp(al-al[1])
lBF=log10(BF)

alldirect = apply(z,1,allones)
lbfall = log10(BF[alldirect==1])

posterior = prior*BF/sum(prior*BF)

#compute marginal bf for each dimension, being direct (lmarginalbf1) or indirect (lmarginalbf2) effect, treating inclusion as indep of other variables
lmarginalbf1 = rep(0,d)
lmarginalbf2 = rep(0,d)
for(i in 1:d){
  condprior = computepriorvec(as.matrix(z[,-i]),pi)
  lmarginalbf1[i] = log10( sum((condprior*BF)[z[,i]==1]) / sum(condprior[z[,i]==1] )) - log10( sum((condprior*BF)[z[,i]==0]) / sum(condprior[z[,i]==0]))
  lmarginalbf2[i]  = log10( sum((condprior*BF)[z[,i]==2]) / sum(condprior[z[,i]==2] )) - log10( sum((condprior*BF)[z[,i]==0]) / sum(condprior[z[,i]==0]))
}

lunivariatebf = lBF[apply(z,1,sum)==(2*d - 1)]

loverallbf = log10(sum(prior[-1]*BF[-1])/sum(prior[-1]))
p0=t(z==0) %*% posterior
p1=t(z==1) %*% posterior
p2=t(z==2) %*% posterior

 return(list(z=z,prior=prior,loglik=al,p=cbind(p0,p1,p2),lbf=lBF,loverallbf = loverallbf, lmarginalbf1= lmarginalbf1,lmarginalbf2=lmarginalbf2,lunivariatebf = lunivariatebf, lbfall=lbfall))
}


bestanova = function(mantest){
  minp=1
  for(i in 1:length(summary.aov(mantest))){
    q = summary.aov(mantest)[[i]]$Pr[1]
    if(q < minp){minp = q}
  }
  return(minp)
}



#Generate multivariate normal with given mean and covariance
rmvnorm = function(n,mu,Sigma){
  r = length(mu)
  Z = matrix(rnorm(n*r),nrow=r)
  L = chol(Sigma)
  return(t(L) %*% Z)
}

# generate a random positive definite matrix
# not the inverse wishart.
rposdef <- function (n, ev = runif(n, 0, 10)) {
Z <- matrix(ncol=n, rnorm(n^2))
decomp <- qr(Z)
Q <- qr.Q(decomp) 
R <- qr.R(decomp)
d <- diag(R)
ph <- d / abs(d)
O <- Q %*% diag(ph)
Z <- t(O) %*% diag(ev) %*% O
return(Z)
}

#Generate random bivariate normal with unit variance and correlation rho
rbinorm = function(n,rho){
  return(rmvnorm(n,c(0,0),cbind(c(1,rho),c(rho,1))))
}
  


sim.2d.data = function(n=500,a=0.3,type){
  d=2
  X = rbinom(n,2,0.2)
  Z = matrix(rnorm(n*d), nrow=d)
  Y = Z
  if(type==1){ # negatively correlated, both affected
    Y[1,] = Z[1,] + a*X
    Y[2,] = (Z[2,] - Z[1,])/sqrt(2)+ a*X
  }
  if(type==2){ # uncorrelated, both affected
    Y[1,] =Z[1,]+a*X
    Y[2,] =Z[2,]+a*X
  }
  if(type==3){ #positive correlated,  both affected, cor = 0.7
    Y[1,] =(0.7*Z[1,]+0.3*Z[2,])/sqrt(0.7^2+0.3^2)+a*X
    Y[2,] =(0.3*Z[1,]+0.7*Z[2,])/sqrt(0.7^2+0.3^2)+a*X
  }
  if(type == 4){ # both affected, Y2 adds nothing (Y1 is intermediate)
    Y[1,] =Z[1,]+a*X
    Y[2,] =Y[1,]+Z[2,]
  }
  if(type == 5){ #Y1 unaffected, Y2 affected, Y1 correlated with Y2
    Y[1,] =Z[1,]
    Y[2,] =(Z[1,]+Z[2,])/sqrt(2)+a*X
  }
  #following types were chosen to illustrate potential "response" associations 
  if(type == 6){ #Y1 and Y2 affected the same
    Y[1,] = (Z[1,] + 0.2*Z[2,])/sqrt(1+0.2^2)+ a * X
    Y[2,] = (Z[1,] - 0.2*Z[2,])/sqrt(1+0.2^2) + a * X
   }
   if(type == 7){ #Y1 only affected
    Y[1,] = (Z[1,] + 0.2*Z[2,])/sqrt(1+0.2^2)+ a * X
    Y[2,] = (Z[1,] - 0.2*Z[2,])/sqrt(1+0.2^2) 
   }
   if(type == 8){ #Y2 only affected
    Y[1,] = (Z[1,] + 0.2*Z[2,])/sqrt(1+0.2^2)
    Y[2,] = (Z[1,] - 0.2*Z[2,])/sqrt(1+0.2^2) + a* X
   }
   if(type == 9){ #Y2 opposite affected than Y1
    Y[1,] = (Z[1,] + 0.2*Z[2,])/sqrt(1+0.2^2) - a * X
    Y[2,] = (Z[1,] - 0.2*Z[2,])/sqrt(1+0.2^2) + a * X
   }
   if(type == 10){ #Y2 differently affected than Y1
    Y[1,] = (Z[1,] + 0.2*Z[2,])/sqrt(1+0.2^2) + 0.5 * a *X
    Y[2,] = (Z[1,] - 0.2*Z[2,])/sqrt(1+0.2^2) + a* X
   }
   return(list(Y=Y,X=X))
}


#compute power at different sizes based on null and alternative simulated test
#statistics. Big values are evidence against null
power = function(null,alt,size){
qq = quantile(null,1-size)
return(1-ecdf(alt)(qq))
}

#needs library ggplot2
#powermatrix = function(s){
#  bf=data.frame(size=size,power=power(s$lbf.0.a1,s$lbf.a1,size),testtype="Multivariate")
#  bfall=data.frame(size=size, power=power(s$lbfall.0.a1,s$lbfall.a1,size),testtype="Multivariate")
#  bfuni = data.frame(size=size,power=power(s$lbfuni.0.a1,s$lbfuni.a1,size),testtype="Univariate")
#  man = data.frame(size=size,power=power(-log(s$p.man.0),-log(s$p.man),size),testtype="Multivariate")
#  aov = data.frame(size=size,power=power(-log(s$p.aov.0),-log(s$p.aov),size),testtype="Univariate")
#  return(melt(list(BFav=bf,BFall=bfall,BFuni=bfuni,ANOVA=aov),id.vars=c("size","power","testtype")))
#}

#library("ggplot2")

#pp1=data.frame(powermatrix(s1.d5.rep),scenario="Multivariate 1")
#pp2 = data.frame(powermatrix(s2.d5.rep),scenario="Independence")
#pp3=data.frame(powermatrix(s3.d5.rep),scenario="Rest Associated")
#pp4 = data.frame(powermatrix(s4.d5.rep),scenario="Latent Factor")
#pp5=data.frame(powermatrix(s5.d5.rep),scenario="Multivariate 2")
#pp6 = data.frame(powermatrix(s6.d5.rep),scenario="Rest Unassociated")

#pp.uni = rbind(pp2,pp6,pp3)
#pp.mv = rbind(pp1,pp5,pp4)
#names(pp.uni)[4] = "Method"
#names(pp.mv)[4] = "Method"

#pp.comb = data.frame(rbind(pp.uni,pp.mv))

#png("power_all6.png",height=360,width=540)
#p=ggplot(data=pp.comb,aes(size,power,colour=Method,lty=testtype))+geom_line() + facet_wrap(~scenario,nrow=2) +opts(axis.text.x  = theme_text(angle=90, colour="gray50",size=10))  + scale_color_manual(values=c("orange","red","blue","black")) 
##+geom_point(shape=16)
#print(p+scale_y_continuous(limits=c(0,1))+scale_x_continuous(limits=c(0,0.05)))

#dev.off()

#pdf("power_all6.pdf",height=6,width=8)
#p=ggplot(data=pp.comb,aes(size,power,colour=Method,lty=testtype))+geom_line() + facet_wrap(~scenario,nrow=2) +opts(axis.text.x  = theme_text(angle=90, colour="gray50",size=10))  + scale_color_manual(values=c("orange","red","blue","black")) 
##+geom_point(shape=16)
#print(p+scale_y_continuous(limits=c(0,1))+scale_x_continuous(limits=c(0,0.05)))

#dev.off()



#pdf("power_oneaffected.pdf",height=2.5,width=6)
#p=ggplot(data=pp.uni,aes(size,power,colour=Method,lty=testtype))+geom_line() + facet_grid(.~scenario) +opts(axis.text.x  = theme_text(angle=90, colour="gray50",size=10))  + scale_color_manual(values=c("orange","red","blue","black")) 
##+geom_point(shape=16)
#print(p+scale_y_continuous(limits=c(0,1))+scale_x_continuous(limits=c(0,0.05)))
#
#dev.off()
#
#pdf("power_multiaffected.pdf",height=2.5,width=6)
#p=ggplot(data=pp.mv,aes(size,power,colour=Method,lty=testtype))+geom_line() + facet_grid(.~scenario) +opts(axis.text.x  = theme_text(angle=90, colour="gray50", size=10))  + scale_color_manual(values=c("orange","red","blue","black")) 
##+geom_point(shape=16)
#print(p+scale_y_continuous(limits=c(0,1))+scale_x_continuous(limits=c(0,0.05)))
#dev.off()









powerplot = function(s,use=c("BF","BFuni","Manova/BFall","Anova"),usecol=c("red","blue","black","blue"),uselty=c(1,1,1,3),maxsize=0.05,...){

  nrep = length(s$lbf.0.a1)
truth = c(rep(0,nrep),rep(1,nrep))
  plot(1,1,type="n",xlim=c(0,maxsize),ylim=c(0,1),xlab="size",ylab="power",...)
for(i in 1:length(use)){
#if(use[i] == "all")
#   bf = c(s$lbfall.0.a1,s$lbfall.a1)
if(use[i] == "BFuni"){
  bf = c(s$lbfuni.0.a1,s$lbfuni.a1)
}
#if(use[i] == "ave")
#  bf = c(log10(10^s$lbfuni.0.a1+10^s$lbfall.0.a1+10^s$lbf.0.a1),log10(10^s$lbfuni.a1+10^s$lbfall.a1 + 10^s$lbf.a1))
if(use[i] == "BF"){
  use[i]="BFav"
  bf = c(s$lbf.0.a1,s$lbf.a1)}
if(use[i] == "BF2"){
  bf = c(10^s$lbf.0.a1+10^s$lbf.0.a2,10^s$lbf.a1+10^s$lbf.a2)}
if(use[i] == "Manova/BFall"){
  use[i]="BFall"
  bf = c(s$lbfall.0.a1,s$lbfall.a1)
#  bf=1-c(s$p.man.0,s$p.man)
}
if(use[i] == "Anova"){
  bf=1-c(s$p.aov.0,s$p.aov)
}
o = order(bf,decreasing=T)
ntrue=cumsum(truth[o])
lines(((1:(2*nrep))-ntrue)/nrep,ntrue/nrep,xlim=c(0,maxsize),lwd=2,col=usecol[i],lty=uselty[i],...)
}

legend(0.03,0.25,cex=0.7,lty=uselty,col=usecol,legend=use,lwd=2)
}






cdfplot=function(s){
  plot(rev(sort(s$lbf.a1)),seq(0,1,length=10000),col=2,type="l",ylim=c(0,0.4))
lines(rev(sort(s$lbfall.a1)),seq(0,1,length=10000),col=1,type="l",ylim=c(0,0.4))
lines(rev(sort(s$lbfuni.a1)),seq(0,1,length=10000),col=4,type="l",ylim=c(0,0.4))
}

#par(mfcol=c(2,3))
#cdfplot(s1.d5.rep)
#cdfplot(s2.d5.rep)
#cdfplot(s3.d5.rep)
#cdfplot(s4.d5.rep)
#cdfplot(s5.d5.rep)
#cdfplot(s6.d5.rep)


#pdf("powerplot_rep_nobf.pdf",width=10,height=6)
#par(mfcol=c(1,2))
#powerplot(s2.d5.rep,use=c("Manova/BFall","Anova"),usecol=c("black","blue"),uselty=c(1,1),main="Univariate effect (Independence)")
#powerplot(s1.d5.rep,use=c("Manova/BFall","Anova"),usecol=c("black","blue"),uselty=c(1,1,1),main="Multivariate effect")
#dev.off()



#pdf("powerplot_rep.pdf",width=10,height=6)
#par(mfcol=c(1,2))
#powerplot(s2.d5.rep,use=c("BF","Manova/BFall","Anova"),usecol=c("red","black","blue"),uselty=c(1,1,1),main="Univariate effect (Independence)")
#powerplot(s1.d5.rep,use=c("BF","Manova/BFall","Anova"),usecol=c("red","black","blue"),uselty=c(1,1,1),main="Multivariate effect")
#dev.off()


#pdf("powerplot_rep2.pdf",width=10,height=4)
#par(mfcol=c(1,3))
#powerplot(s3.d5.rep,use=c("BF","Manova/BFall","Anova"),usecol=c("red","black","blue"),uselty=c(1,1,1),main="Univariate effect (rest downstream)")
#powerplot(s6.d5.rep,use=c("BF","Manova/BFall","Anova"),usecol=c("red","black","blue"),uselty=c(1,1,1),main="Univariate effect (rest upstream)")
#powerplot(s4.d5.rep,use=c("BF","Manova/BFall","Anova"),usecol=c("red","black","blue"),uselty=c(1,1,1),main="Multivariate effect (latent factor)")
#dev.off()


#pdf("powerplot_ashg.pdf",width=10,height=4)
#par(mfcol=c(1,3))
#powerplot(s2.d5.rep,use=c("Manova/BFall","BFuni","BF"),usecol=c("red","black","blue"),uselty=c(1,1,1),main="Univariate effect (Independence)")
#powerplot(s1.d5.rep,use=c("Manova/BFall","BFuni","BF"),usecol=c("red","black","blue"),uselty=c(1,1,1),main="Multivariate effect")
#powerplot(s4.d5.rep,use=c("Manova/BFall","BFuni","BF"),usecol=c("red","black","blue"),uselty=c(1,1,1),main="Multivariate effect (latent factor)")
#dev.off()

#pdf("powerplot.rep.pdf",width=10,height=6)
#par(mfcol=c(2,3))
#powerplot(s2.d5.rep,main="Univariate effect (Independence)")
#powerplot(s1.d5.rep,main="Multivariate effect")
#powerplot(s3.d5.rep,main="Univariate effect (rest downstream)")
#powerplot(s4.d5.rep,main="Multivariate (latent factor) effect")
#powerplot(s5.d5.rep,main="2-comp effect (rest upstream)")
#powerplot(s6.d5.rep,main="Univariate effect (rest upstream)")
#dev.off()






#par(mfcol=c(2,3))
#powerplot(s2.d5.rep,main="Univariate effect, Independence")
#powerplot(s1.d5.rep,main="Multivariate effect")
#powerplot(s3.d5.rep,main="Univariate effect, rest downstream")
#powerplot(s4.d5.rep,main="Multivariate factor effect")
#powerplot(s5.d5.rep,main="2 components affected")
#powerplot(s6.d5.rep,main="1 component affected")

#some numerics, in case it is interesting
#> mean(s1$lbf.a1> quantile(s1$lbf.0.a1,0.95))
#[1] 0.9422
#> mean(s1$lbf.a4> quantile(s1$lbf.0.a4,0.95))
#[1] 0.9315
#> mean(s1$p.aov < quantile(s1$p.aov.0,0.05))
#[1] 0.7715
#> mean(s1$p.man < quantile(s1$p.man.0,0.05))
#[1] 0.9453
#
#> mean(s2$lbf.a1> quantile(s2$lbf.0.a1,0.95))
#[1] 0.5005
#> mean(s2$lbf.a4> quantile(s2$lbf.0.a4,0.95))
#[1] 0.5083
#> mean(s2$p.aov < quantile(s2$p.aov.0,0.05))
#[1] 0.5059
#> mean(s2$p.man < quantile(s2$p.man.0,0.05))
#[1] 0.4581





fplot=function(Y,X,...){
  plot(Y[1,],Y[2,],col=X+1,ylab="Y2",xlab="Y1",...)
points( lapply(split(Y[1,],X),mean),lapply(split(Y[2,],X),mean),col=c(1,2,3),cex=4,pch=20)
}

#s1.2d = sim.2d.data(type=1,a=1,n=500)
#s2.2d = sim.2d.data(type=2,a=1,n=500)
#s3.2d = sim.2d.data(type=3,a=1,n=500)
#s4.2d = sim.2d.data(type=4,a=1,n=500)
#s5.2d = sim.2d.data(type=5,a=1,n=500)
#s6.2d = sim.2d.data(type=6,a=1,n=500)

#pdf("bothassoc_negcorr.pdf",width=7,height=7)
#fplot(s1.2d$Y,s1.2d$X,main="Both affected, negative correlation")
#dev.off()

#pdf("bothassoc_nocorr.pdf",width=7,height=7)
#fplot(s2.2d$Y,s2.2d$X,main="Both affected, no correlation")
#dev.off()

#pdf("bothassoc_poscorr.pdf",width=7,height=7)
#fplot(s3.2d$Y,s3.2d$X,main="Both affected, positive correlation")
#dev.off()
#
#pdf("condindep.pdf",width=7,height=7)
#fplot(s4.2d$Y,s4.2d$X, main="Both affected, Y2 adds nothing (Y2=Y1+noise)")
#dev.off()

#pdf("oneassoc.pdf", width=7,height=7)
#fplot(s5.2d$Y,s5.2d$X, main="Y2 only affected")
#dev.off()

#pdf("sameassoc.pdf",width=7,height=7)
#fplot(s6.2d$Y,s6.2d$X, main="Y1 and Y2 affected the same")
#dev.off()

#pdf("bothassoc.pdf",width=14,height=7)
#par(mfcol=c(1,2))
#fplot(s1.2d$Y,s1.2d$X,main="Both affected, negative correlation")
#fplot(s2.2d$Y,s2.2d$X, main="Both affected, no correlation")
#dev.off()



#sim summary
#  par(mfcol=c(2,5),mar=c(2,2,2,2))

#does 100 2d simulations of given type and applies MVassoctest to them
simset=function(type){
s = list()
t.04=list()
t.01=list()
for(i in 1:100){
  s[[i]] = sim.2d.data(type=type)
  t.04[[i]] = MVassoctest(s[[i]]$Y,s[[i]]$X,sigmaa=0.4)
  t.01[[i]] = MVassoctest(s[[i]]$Y,s[[i]]$X,sigmaa=0.1)
}
return(list(s=s,t.01=t.01,t.04=t.04))
}

mv.vs.uni.compplot=function(sim,...){
t=sim$t.04
nrep = length(t)
mvbf = rep(0,nrep)
unibf = rep(0,nrep)
for(j in 1:100){mvbf[j]=t[[j]]$loverallbf;
                unibf[j] = log10(mean(10^t[[j]]$lunivariatebf));}
xmax = max(c(mvbf,unibf))
xmin = min(c(mvbf,unibf))
if(xmax<12){xmax=12}
if(xmin>-1){xmin=-1}
plot(unibf,mvbf,xlim=c(xmin,xmax),ylim=c(xmin,xmax),xlab="log10(mean univariate BF)",ylab="log10(multivarate BF)",...)
abline(a=0,b=1)
abline(v=4,h=4,lty=2)
}

mv.vs.all.compplot=function(sim,...){
t=sim$t.04
nrep = length(t)
mvbf = rep(0,nrep)
unibf = rep(0,nrep)
for(j in 1:100){mvbf[j]=t[[j]]$loverallbf;
                unibf[j] = t[[j]]$lbfall;}
xmax = max(c(mvbf,unibf))
xmin = min(c(mvbf,unibf))
if(xmax<12){xmax=12}
if(xmin>-1){xmin=-1}
plot(unibf,mvbf,xlim=c(xmin,xmax),ylim=c(xmin,xmax),xlab="log10(BFall)",ylab="log10(multivarate BF)",...)
abline(a=0,b=1)
abline(v=4,h=4,lty=2)
}


#sim1 = simset(1)
#sim2=simset(2)
#sim3=simset(3)
#sim4=simset(4)
#sim5=simset(5)

#pdf("mv_vs_uni_1.pdf",width=7,height=7)
#mv.vs.uni.compplot(sim1,main="Both affected, negative correlation")
#dev.off()
#pdf("mv_vs_uni_2.pdf",width=7,height=7)
#mv.vs.uni.compplot(sim2,main="Both affected, no correlation")
#dev.off()
#pdf("mv_vs_uni_3.pdf",width=7,height=7)
#mv.vs.uni.compplot(sim3,main="Both affected, positive correlation")
#dev.off()
#pdf("mv_vs_uni_4.pdf",width=7,height=7)
#mv.vs.uni.compplot(sim4,main="Both affected, Y2 adds nothing")
#dev.off()
#pdf("mv_vs_uni_5.pdf",width=7,height=7)
#mv.vs.uni.compplot(sim5,main="Y2 only affected, Y1 correlated with Y2")
#dev.off()


#pdf("mv_vs_all_1.pdf",width=7,height=7)
#mv.vs.all.compplot(sim1,main="Both affected, negative correlation")
#dev.off()
#pdf("mv_vs_all_2.pdf",width=7,height=7)
#mv.vs.all.compplot(sim2,main="Both affected, no correlation")
#dev.off()
#pdf("mv_vs_all_3.pdf",width=7,height=7)
#mv.vs.all.compplot(sim3,main="Both affected, positive correlation")
#dev.off()
#pdf("mv_vs_all_4.pdf",width=7,height=7)
#mv.vs.all.compplot(sim4,main="Both affected, Y2 adds nothing")
#dev.off()
#pdf("mv_vs_all_5.pdf",width=7,height=7)
#mv.vs.all.compplot(sim5,main="Y2 only affected, Y1 correlated with Y2")
#dev.off()




#x1 = rep(0,100)
#x2=x1
#y=x1
#for(j in 1:100){x1[j]=t[[j]]$lmarginalbf1[2]; x2[j] = t[[j]]$lmarginalbf2[2]; y[j]=t[[j]]$lunivariatebf[2]}
#hist(y-x1)
#hist(y-x2)

#t1=MVassoctest(Y1,X)
#t2=MVassoctest(Y2,X)
#t3=MVassoctest(Y3,X)
#t4=MVassoctest(Y4,X)
#t5=MVassoctest(Y5,X)

simple.5d.simdata = function(a,n,type=1){
X = rbinom(n,2,0.2)
d=5
Z = matrix(rnorm(n*d), nrow=d)
Y=Z
f=NULL
if(type==1){ 
Y[1,] =Z[1,]+a*X
Y[2,] =(Z[2,]+Z[1,])/sqrt(2)
Y[3,]=Z[3,]+a*X
Y[4,]=(Z[4,]+Z[2,]+Z[3,])/sqrt(3)
Y[5,]=(Y[3,]+Z[5,])/sqrt(2)
}
if(type == 2){ #Ys are uncorrelated, only Y1 is affected
Y[1,] =Z[1,]+a*X
Y[2,] =Z[2,]
Y[3,]=Z[3,]
Y[4,]=Z[4,]
Y[5,]= Z[5,]
}
if(type==3){ #here Y1 is affected by X, and rest are all downstream
Y[1,] =Z[1,]+a*X
Y[2,] =(Z[2,]+Y[1,])/sqrt(2)
Y[3,]=(Z[3,]+0.5*Y[1,])/sqrt(1+0.5^2)
Y[4,]=(Z[4,]+0.2*Y[1,])/sqrt(1+0.2^2)
Y[5,]= (Z[5,]+0.1*Y[1,])/sqrt(1+0.1^2)
}
if(type==4){ # here a single factor is really responsible for association
f = a*X + 0.5*rnorm(n)
Y[1,] = 0.5*Z[1,] + 0.3*f
Y[2,]= 0.5*Z[2,] + 0.2*f
Y[3,] = 0.5*Z[3,] - 0.3*f
Y[4,] = 0.5*Z[4,] + 0.5*f
Y[5,] = 0.5*Z[5,]+ 0.2*f
}
if(type==5){ # here it is a 2-component model, with other 3 being unaffected
Y[1,] = Z[1,]
Y[2,]=  Z[2,]
Y[3,] = Z[3,]
Y[4,] = (Z[4,] + Y[1,] +Y[2,])/sqrt(3)+ a*X
Y[5,] = (Z[5,] + Y[3,])/sqrt(2) + a*X
}
if(type==6){ # here it is a 1-component model, with others being unaffected
Y[1,] = Z[1,]
Y[2,]=  Z[2,]
Y[3,] = Z[3,]
Y[4,] = Z[4,]
Y[5,] = (Z[5,] + Y[1,] + Y[2,] + Y[3,] + Y[4,])/sqrt(5) + a*X
}

return(list(Y=Y,X=X,f=f))
}


intermediate.5d.simdata = function(a,n){
X = rbinom(n,2,0.2)
d=5
Z = matrix(rnorm(n*d), nrow=d)
Y=Z
Y[1,] =Z[1,]+a*X
Y[2,] =Z[2,]+a*X
Y[3,]=Z[3,]+a*X
Y[4,]=Z[4,]+a*X
Y[5,]=Z[5,] + 0.25*(Y[1,]+Y[2,]+Y[3,]+Y[4,])
return(list(Y=Y,X=X))
}


simple.5d.simdata2 = function(a,n){
X = rbinom(n,2,0.2)
d=5
Z = matrix(rnorm(n*d), nrow=d)
Y=Z
Y[1,] =Z[1,]+a*X
Y[2,] =Z[2,]
Y[3,]=Z[3,]
Y[4,]=Z[4,]
Y[5,]= Z[5,]
return(list(Y=Y,X=X))
}



#sbf = rep(0,100)
#bf=sbf
#prob = list()
#for(j in 1:100){
#sim = simple.5d.simdata(0.3,500)
#Y = sim$Y
#X = sim$X
#for(i in  1:d){Y[d,]=Y[d,]/sd(Y[d,])}
#t=MVassoctest(Y,X)
## single SNP BFs
#sbf[j] = max(log10(exp(t$loglik[which(apply(t$z,1,sum)==9)] - t$loglik[1])))
#bf[j] = t$lbf
#prob[[j]] = t$p
#}
  
#AssocLogLik(Y,X,1:5,NULL,NULL,0.4,5)

#X=matrix(rep(1,n),nrow=1)
#K = matrix(1/sigma0sq,nrow=1,ncol=1)
#S0= diag(n0,d)
#
#LogMarginalLikwithconst(Y,X,K,n0,S0)
#uset=(1:3)
#Yu=Y[uset,]
#Ya=Y[-uset,]
#
#LogMarginalLikwithconst(Yu,X,K,n0-nrow(Ya),S0[uset,uset]) + LogMarginalLikwithconst(Ya,rbind(X,Yu),diag(c(1/sigma0sq,diag(S0[uset,uset]))),n0,S0[-uset,-uset])


#n0=10


#uset=c(1,2)
#aset=c(3,4)
#dset=c(5,6)
#n=100
#d=6
#Y = matrix(rnorm(n*d), nrow=d)

                        
#compute log marginal likelihood for multivariate regression
#from Minka paper, up to a constant
LogMarginalLik=function(Y,X,K,N0){
  d = dim(Y)[1]
  N = dim(Y)[2]
  S0 = N0 * diag(d)
  SXX = X %*% t(X) + K
  SYY = Y %*% t(Y)
  SYX = Y %*% t(X)
  SYgX = SYY - SYX %*% chol2inv(chol(SXX)) %*% t(SYX)
  return((d/2) * (determinant(K)$modulus - determinant(SXX)$modulus) + (N0/2) * determinant(S0)$modulus - (N+N0)/2 * determinant(SYgX + S0)$modulus)
}


#testing marginal likelihood calculation
#N=1000
#K1 = diag(c(0.00001,4))
#K0 = diag(c(0.00001,1e10))
#N0 = 0.001
#
#N=100
#g = rbinom(N,2,0.5)
#a = 1
#y1 = rnorm(N)
#y2 = 0.1*a*g
#Y = rbind(y1,y1+y2)
#X = rbind(rep(1,N),g)
#LogMarginalLik(Y,X,K0,N0)-LogMarginalLik(Y,X,K1,N0)

#check transformation invariant
#LogMarginalLik(rbind(y1,y2),X,K0,N0)-LogMarginalLik(rbind(y1,y2),X,K1,N0)

#check loss of information from including irrelevant info
#LogMarginalLik(rbind(y2),X,K0,N0)-LogMarginalLik(rbind(y2),X,K1,N0)
#
## try some null simulations
#temp = rnorm(N)
#YNULL = rbind(temp,temp+rnorm(N))
#LogMarginalLik(YNULL,X,K0,N0)-LogMarginalLik(YNULL,X,K1,N0)

# Examine X+Y and X-Y transform

#N=1000
#g = rbinom(N,2,0.5)
#a=0.3
#z = rnorm(N)
#y1 = z + rnorm(N,0,0.7)
#y2 = z + rnorm(N,0,0.7) + a*g
#Y = rbind(y1,y2)
#X = rbind(rep(1,N),g)
#LogMarginalLik(Y,X,K0,N0)-LogMarginalLik(Y,X,K1,N0)

#Y = rbind(y1+y2)
#LogMarginalLik(Y,X,K0,N0)-LogMarginalLik(Y,X,K1,N0)
#Y = rbind(y1-y2)
#LogMarginalLik(Y,X,K0,N0)-LogMarginalLik(Y,X,K1,N0)




#try repeated measures approach
#N=100
#g = rbinom(N,2,0.5)
#X = cbind(rep(1,2*N),rbind(diag(N),diag(N)),c(g,rep(0,N)),c(g,g))
#mu = rnorm(N)
#Y1 = mu1+


#PARC analysis examples (not subfractions)

#parc.snps = read.table("PARCmvphenotypes/rs.MS.Sep09.txt",header=F)
#postmean.geno=function(x){
#y = as.numeric(x[-(1:3)])
#ind = as.vector( t(cbind(1:(length(y)/2),1:(length(y)/2) )) )
#return(as.vector( by(y,ind, function(x) {x[2] + 2*(1-x[1] - x[2]) } ) ))
#}
#parc.postmean = matrix(nrow= dim(parc.snps)[1],ncol=(dim(parc.snps)[2]-3)/2)
#for(i in 1:dim(parc.snps)[1]){
#parc.postmean[i,] = postmean.geno(parc.snps[i,])
#}

# parc.crp = scan("PARCmvphenotypes/log_CRP_rqav_noID.gd.txt")
# parc.ldlc = scan("PARCmvphenotypes/log_LDLC_rqav_noID.gd.txt")
# parc.hdlc = scan("PARCmvphenotypes/log_HDLC_rqav_noID.gd.txt")
# parc.tg = scan("PARCmvphenotypes/log_Tg_rqav_noID.gd.txt")
# parc.tc = scan("PARCmvphenotypes/log_TC_rqav_noID.gd.txt")

#parc.Y = rbind(parc.crp,parc.ldlc,parc.hdlc,parc.tg)
#parc.t = list()
#for(i in 1:5){
#parc.t[[i]] = MVassoctest(parc.Y,parc.postmean[i,],sigmaa=0.2)}

#parc.Y2 = rbind(parc.crp,parc.ldlc,parc.hdlc,parc.tg,parc.tc)
#parc.t2 = list()
#for(i in 1:5){
#parc.t2[[i]] = MVassoctest(parc.Y2,parc.postmean[i,],sigmaa=0.2)}

#pimage=function(t,thresh=0.1,...){
#p = t$p
#image(c(1,2,3),1:dim(p)[1],t(p),col=gray(seq(1,0,length=10)),tck=0,xlab="",ylab="",labels=FALSE,...)
#axis(1,at=1:3,labels=c("Unaffected","Direct","Indirect"))
#labs=c("CRP","LDLC","HDLC","Tg","TC")
#if(dim(p)[1]==4){labs = c("CRP","LDLC","HDLC","Tg")}
#axis(2,at=1:dim(p)[1],labels=labs)
#for(i in 2:3){
#  for(j in 1:dim(p)[1]){
#    if(p[j,1]<thresh){
#text(i,j,t(round(p[j,i],2)),col=2,cex=2)
#}}}
#}

#pdf("pplot_cetp.pdf",width=4,height=7)
#pimage(parc.t[[1]],main="CETP")
#dev.off()
#pdf("pplot_crp.pdf",width=4,height=7)
#pimage(parc.t[[2]],main="CRP")
#dev.off()
#pdf("pplot_apoe.pdf",width=4,height=7)
#pimage(parc.t[[3]],main="APOE")
#dev.off()
#pdf("pplot_gckr.pdf",width=4,height=7)
#pimage(parc.t[[4]],main="GCKR")
#dev.off()
#pdf("pplot_tcf1.pdf",width=4,height=7)
#pimage(parc.t[[5]],main="TCF1")
#dev.off()




#pdf("pplot2_cetp.pdf",width=4,height=7)
#pimage(parc.t2[[1]],main="CETP")
#dev.off()
#pdf("pplot2_crp.pdf",width=4,height=7)
#pimage(parc.t2[[2]],main="CRP")
#dev.off()
#pdf("pplot2_apoe.pdf",width=4,height=7)
#pimage(parc.t2[[3]],main="APOE")
#dev.off()
#pdf("pplot2_gckr.pdf",width=4,height=7)
#pimage(parc.t2[[4]],main="GCKR")
#dev.off()
#pdf("pplot2_tcf1.pdf",width=4,height=7)
#pimage(parc.t2[[5]],main="TCF1")
#dev.off()


#pdf("parc_pairs.pdf",width=7,height=7)
#pairs(t(parc.Y2),labels=c("CRP","LDLC","HDLC","Tg","TC"))
#dev.off()


#PARC analysis examples 2nd try (not subfractions)

#parc.snps2 = read.table("PARCmvphenotypes/rs.MS.Sep09II.txt",header=F)
postmean.geno=function(x){
y = as.numeric(x[-(1:3)])
ind = as.vector( t(cbind(1:(length(y)/2),1:(length(y)/2) )) )
return(as.vector( by(y,ind, function(x) {x[2] + 2*(1-x[1] - x[2]) } ) ))
}
#parc.postmean2 = matrix(nrow= dim(parc.snps2)[1],ncol=(dim(parc.snps2)[2]-3)/2)
#for(i in 1:dim(parc.snps2)[1]){
#parc.postmean2[i,] = postmean.geno(parc.snps2[i,])
#}

# parc.crp = scan("PARCmvphenotypes/log_CRP_rqav_noID.gd.txt")
# parc.ldlc = scan("PARCmvphenotypes/log_LDLC_rqav_noID.gd.txt")
# parc.hdlc = scan("PARCmvphenotypes/log_HDLC_rqav_noID.gd.txt")
# parc.tg = scan("PARCmvphenotypes/log_Tg_rqav_noID.gd.txt")
# parc.tc = scan("PARCmvphenotypes/log_TC_rqav_noID.gd.txt")

#parc.Y2 = rbind(parc.crp,parc.ldlc,parc.hdlc,parc.tg,parc.tc)
#parc.t2 = list()
#for(i in 1:9){
#parc.t2[[i]] = MVassoctest(parc.Y2,parc.postmean2[i,],sigmaa=0.2)}


#pdf("pplot_apob_II.pdf",width=4,height=7)
#pimage(parc.t2[[1]],main="APOB")
#dev.off()
#pdf("pplot_gckr_II.pdf",width=4,height=7)
#pimage(parc.t2[[2]],main="GCKR")
#dev.off()
#pdf("pplot_ldlr_II.pdf",width=4,height=7)
#pimage(parc.t2[[3]],main="LDLR")
#dev.off()
#pdf("pplot_apoe_II.pdf",width=4,height=7)
#pimage(parc.t2[[4]],main="APOE")
#dev.off()
#pdf("pplot_crp_II.pdf",width=4,height=7)
#pimage(parc.t2[[5]],main="CRP")
#dev.off()
#pdf("pplot_tcf1_II.pdf",width=4,height=7)
#pimage(parc.t2[[6]],main="TCF1")
#dev.off()
#pdf("pplot_cetp_II.pdf",width=4,height=7)
#pimage(parc.t2[[7]],main="CETP")
#dev.off()
#pdf("pplot_lpl_II.pdf",width=4,height=7)
#pimage(parc.t2[[8]],main="LPL")
#dev.off()
#pdf("pplot_lipc_II.pdf",width=4,height=7)
#pimage(parc.t2[[9]],main="LIPC")
#dev.off()






#pdf("pplot2_cetp.pdf",width=4,height=7)
#pimage(parc.t2[[1]],main="CETP")
#dev.off()
#pdf("pplot2_crp.pdf",width=4,height=7)
#pimage(parc.t2[[2]],main="CRP")
#dev.off()
#pdf("pplot2_apoe.pdf",width=4,height=7)
#pimage(parc.t2[[3]],main="APOE")
#dev.off()
#pdf("pplot2_gckr.pdf",width=4,height=7)
#pimage(parc.t2[[4]],main="GCKR")
#dev.off()
#pdf("pplot2_tcf1.pdf",width=4,height=7)
#pimage(parc.t2[[5]],main="TCF1")
#dev.off()


#pdf("parc_pairs.pdf",width=7,height=7)
#pairs(t(parc.Y2),labels=c("CRP","LDLC","HDLC","Tg","TC"))
#dev.off()
