# list of functions which are used by R codes.
assig <- function(n_args){
  cargs <- vector("list", length(n_args))
  for(i in 1:length(n_args)) cargs[[i]] <- 1:n_args[i]
  t(expand.grid(cargs))
}

##--------------produce the B-spline functions----------------------##
bsbasefun <- function(X,K,degr){
  n = dim(X)[1]
  p = dim(X)[2]
  nk = K - degr
  u.k = seq(0, 1, length=nk+2)[-c(1,nk+2)]
  BS = NULL
  for(j in 1:p){
    Knots = as.numeric(quantile(X[,j], u.k))  
    BS0 = bs(X[,j], knots=Knots, intercept=TRUE, degree=degr)
    BS = cbind(BS,BS0[,-1])
  }
  BS = scale(BS,center = T, scale = F)
  id = seq(1,p*K,K)
  Z = NULL
  for(j in 1:K){
    Z = cbind(Z,BS[,id+j-1])
  }
  return(Z)
}

##--------------produce p*q estimated functions----------------------##
trans <- function(X,D3,p,q,K,degr,s0){
  Z = bsbasefun(X,K,degr)
  id = seq(1,p*K,p) 
  fjl = assig(c(s0,q)) 
  funhat = NULL
  for(j0 in 1:(q*s0)){ 
    qj1 = fjl[1,j0]
    qj = fjl[2,j0]
    funhat = cbind(funhat,Z[,id+qj1-1]%*%D3[qj,id+qj1-1]) 
  }
  return(funhat)
}

truefuns <- function(n,q,p,s,D2,seed_id){
  set.seed(seed_id)
  X <- matrix(runif(n*p),n,p)
  X1 <- X[,1:s]
  basefuns1 <- sin(2*pi*X1)
  basefuns2 <- cos(pi*X1)
  f0 <- matrix(rep(basefuns1,q),nrow=n)*matrix(rep(D2[1,],each=n),n) + matrix(rep(basefuns2,q),nrow=n)*matrix(rep(D2[2,],each=n),n)
  return(list(X=X,f0=f0))
  
}
##--------------generate data----------------------##
generateData <- function(n,q,p,s,D2,group,sigma2=NULL,seed_id=NULL, rho=0.0){
  if(n<2) stop("n must be not smaller than 2")
  if(q<1) stop("q must be not smaller than 1")
  if(p<1) stop("p must be not smaller than 1")
  if(s<1) stop("s must be not smaller than 1")
  pq <- dim(D2)
  if(is.null(pq)|(pq[1]<=1)|(pq[2]<=1))stop("D2 should be a matrix with 2 or more rows and columns")
  if(is.null(sigma2)) sigma2 = 0.1
  if(sigma2<=0){
    warning("sigma2 <= 0; set to 0.1")
    sigma2=0.1
  }
  if(is.null(seed_id)) seed_id=1000
  set.seed(seed_id)
  X <- matrix(runif(n*p), nrow = n)

  gunique <- unique(group)
  ng = length(gunique)
  fl = matrix(0,n,q)
  f0 = matrix(0,n,q)
  
  for(j in 1:ng){
    #Xj <- X[,((j-1)*len+1):(j*len)]
    Xj <- X[,group==gunique[j]]
    X1 <- Xj[,1:s]
    basefuns1 <- sin(2*pi*X1)
    basefuns2 <- cos(pi*X1)
    fl <- fl + basefuns1%*%matrix(D2[1,],nrow=s)+basefuns2%*%matrix(D2[2,],nrow=s)
    #f0 <- f0 + matrix(rep(basefuns1,q),nrow=n)*matrix(rep(D2[1,],each=n),n) + matrix(rep(basefuns2,q),nrow=n)*matrix(rep(D2[2,],each=n),n)
  }
  eps <- matrix(rnorm(n*q),n,q)
  Y <- fl + sigma2*eps
  return(list(Y=Y,X=X,f0=f0))
}

##--------------generate data----------------------##
generateData2 <- function(n,q,p,s,D3,K,degr,group,sigma2,seed_id, t=0){
  if(n<2) stop("n must be not smaller than 2")
  if(q<1) stop("q must be not smaller than 1")
  if(p<1) stop("p must be not smaller than 1")
  if(s<1) stop("s must be not smaller than 1")
  pq <- dim(D3)
  if(is.null(pq)|(pq[1]<=1)|(pq[2]<=1))stop("D3 should be a matrix with 2 or more rows and columns")
  if(is.null(sigma2)) sigma2 = 0.1
  if(sigma2<=0){
    warning("sigma2 <= 0; set to 0.1")
    sigma2=0.1
  }
  if(is.null(seed_id)) seed_id=1000
  set.seed(seed_id)
  W <- matrix(runif(n*p), nrow = n)
  U <- runif(n)
  X <- (W+t*matrix(rep(U,p),nrow = n))/(1+t)
  
  gunique <- unique(group)
  ng = length(gunique)
  fl = matrix(0,n,q)
  for(j in 1:ng){
    Xj <- X[,group==gunique[j]]
    X1 <- Xj[,1:s]
    fl = fl + bsbasefun(X1,K,degr)%*%t(D3)
  }
  f0 = fl
  eps <- matrix(rnorm(n*q),n,q)
  Y <- fl + eps*sigma2
  return(list(Y=Y,X=X,f0=f0))
}

##--------------generate data----------------------##
generateDataT4 <- function(n,q,p,s,D42,group,sigma2=NULL,seed_id=NULL, rho=0.0){
  if(n<2) stop("n must be not smaller than 2")
  if(q<1) stop("q must be not smaller than 1")
  if(p<1) stop("p must be not smaller than 1")
  if(s<1) stop("s must be not smaller than 1")
  pq <- dim(D42)  #???
  if(is.null(pq)|(pq[1]<=1)|(pq[2]<=1))stop("D2 should be a matrix with 2 or more rows and columns")
  if(is.null(sigma2)) sigma2 = 0.1
  if(sigma2<=0){
    warning("sigma2 <= 0; set to 0.1")
    sigma2=0.1
  }
  if(is.null(seed_id)) seed_id=1000
  set.seed(seed_id)
  X <- matrix(runif(n*p), nrow = n)
  
  gunique <- unique(group)
  ng = length(gunique)
  len = p/ng
  fl = matrix(0,n,q)
  f0 = matrix(0,n,q)
  
  id = NULL
  for(j in 1:q) id = c(id, c(1:len)+(j-1)*len*ng)
  for(j in 1:ng){
    D2 = D42[,id+(j-1)*len]
    Xj <- X[,group==gunique[j]]
    X1 <- Xj[,1:s]
    basefuns1 <- sin(2*pi*X1)
    basefuns2 <- cos(pi*X1)
    fl <- fl + basefuns1%*%matrix(D2[1,],nrow=s)+basefuns2%*%matrix(D2[2,],nrow=s)
  }
  eps <- matrix(rnorm(n*q),n,q)
  Y <- fl + sigma2*eps
  return(list(Y=Y,X=X,f0=f0))
}

##--------------generate data----------------------##
generateData2T4 <- function(n,q,p,s,D44,K,degr,group,sigma2,seed_id, t=0){
  if(n<2) stop("n must be not smaller than 2")
  if(q<1) stop("q must be not smaller than 1")
  if(p<1) stop("p must be not smaller than 1")
  if(s<1) stop("s must be not smaller than 1")
  pq <- dim(D44)
  if(is.null(pq)|(pq[1]<=1)|(pq[2]<=1))stop("D3 should be a matrix with 2 or more rows and columns")
  if(is.null(sigma2)) sigma2 = 0.1
  if(sigma2<=0){
    warning("sigma2 <= 0; set to 0.1")
    sigma2=0.1
  }
  if(is.null(seed_id)) seed_id=1000
  set.seed(seed_id)
  W <- matrix(runif(n*p), nrow = n)
  U <- runif(n)
  X <- (W+t*matrix(rep(U,p),nrow = n))/(1+t)
  
  gunique <- unique(group)
  G = length(gunique)
  len = s
  fl = matrix(0,n,q)
  for(j in 1:G){
    D3 = D44[,((j-1)*len*K+1):(j*len*K)]
    Xj <- X[,group==gunique[j]]
    X1 <- Xj[,1:s]
    fl = fl + bsbasefun(X1,K,degr)%*%t(D3)
  }
  f0 = fl
  eps <- matrix(rnorm(n*q),n,q)
  Y <- fl + eps*sigma2
  return(list(Y=Y,X=X,f0=f0))
}
##--------------plot curve of function f_{jl} ----------------------##
plotfuns <- function(fit,funTrueID,true.curve=FALSE){
  # funTrueID = c(j,l) is the index of the f_{jl}th function
  qj1 = funTrueID[1]
  qj = funTrueID[2]
  j0 = s0*(qj-1) + qj1
  
  n = nrow(fit$Y)
  p = ncol(fit$X)
  q = ncol(fit$Y)
  s0 = fit$s0
  X = fit$X0
  n0 = nrow(X)
  funhat = trans(X,fit$Dnew,p,q,fit$K,fit$degr,s0)
  xs = sort(X[,qj1],index.return = T)
  w = xs$x
  iw = xs$ix
  
  if(true.curve){
    D2 = fit$D2
    
    X1 <- X[,1:s0]
    basefuns1 <- sin(2*pi*X1)
    basefuns2 <- cos(pi*X1)
    f0 <- matrix(rep(basefuns1,q),nrow=n0)*matrix(rep(D2[1,],each=n0),n0)+matrix(rep(basefuns2,q),nrow=n0)*matrix(rep(D2[2,],each=n0),n0)

    if(qj1>s0)  f11 = rep(0,n)
    if(qj1<=s0) f11 = f0[,j0]
    
    f11hat = funhat[,j0]
    plot(w,f11[iw],type = "l",col = "blue",ylim=c(min(f11[iw])-0.2,max(f11[iw])+0.2),ylab="f(x)",xlab="x",
         main="True: blue-dot, Est: red-solid",pch=1,lty=3,lwd=3,cex=3)
    lines(w,f11hat[iw],col="red",pch=3,lty=1,lwd=3)
  }
  else{
    f11hat = funhat[,j0]
    plot(w,f11hat[iw],type = "l",col = "red",ylab="f(x)",xlab="x",main="Estimated curve",pch=1,lty=1,lwd=3,cex=3)
  }
}



