f.rbing=function(n,lam){
  ## lam contains the q-1 non negative eigenvalues
  lam=sort(lam,decreasing=TRUE) ## sort the eigenvalues in desceding order
  nsamp=0
  X=NULL
  lam.full=c(lam,0)
  q=length(lam.full)
  A=diag(lam.full)
  SigACG.inv=diag(q)+2*A
  SigACG=solve(SigACG.inv)
  Ntry=0
  while(nsamp < n) {
  x.samp=FALSE
  while(x.samp==FALSE) {
  yp=MASS::mvrnorm(n=1,mu=rep(0,q),Sig=SigACG)
  y=yp/sqrt(t(yp)%*%yp)
  lratio=-t(y)%*%A%*%y -q/2*log(q) + 0.5*(q-1) + q/2*log(t(y)%*%SigACG.inv%*%y)
  if(log(runif(1)) < lratio) {
  X=c(X,y)
  x.samp=TRUE
  nsamp = nsamp+1
  }
  Ntry=Ntry+1 }
  }
  if(n>1) X=matrix(X,byrow=T,ncol=q)
  ## the X contains the simulated values
  ## the avtry is the estimate of the M in rejection sampling
  ## 1/M is the probability of acceptance
  list(X=X,avtry=Ntry/n) 
}

rbingham=function(n,A){
  p = 3
  lam=eigen(A)$values
  V=eigen(A)$vectors
  lam=lam-lam[p]
  lam=lam[-p]
  x=f.rbing(5*n,lam)$X
  y=x%*%t(V)
  y 
}

fb.sim=function(n,k,m,A) {
  ## n is the required sample size
  ## k is the concentration parameter, the Fisher part
  ## m is the mean vector, the Fisher part
  ## A is the symmetric matrix, the Bingham part
  q=length(m)
  A1=A+k/2*( diag(q)-m%*%t(m) )
  lam=eigen(A1)$values
  V=eigen(A1)$vectors
  lam=lam-lam[q]
  lam=lam[-q]
  x=f.rbing(5*n,lam)$X
  x=x%*%t(V)
  u=log(runif(5*n))
  ffb=k*x%*%m-diag(x%*%A%*%t(x))
  fb=k-diag(x%*%A1%*%t(x))
  ina=1:c(5*n)
  keep=ina[u<=c(ffb-fb)]
  ind=sample(keep,n)
  y=x[ind,]
  y
}

n <- 100
k <- 100
m <- c(0,0,1)
A <- matrix(c(-47.5,0,0,0,47.5,0,0,0,0),nrow=3,ncol=3)
y <- fb.sim(n,k,m,A)
write(t(y),"rdata",sep="\t",ncolumns=3)

