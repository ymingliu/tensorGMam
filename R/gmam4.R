
##--------------without sparsity----------------------##
gmam4 <- function(Y,X,group=NULL,K=6,r1=NULL,r2=NULL,r3=NULL,r4=NULL,D0=NULL,intercept=TRUE,degr=3,eps=1e-4,max_step=20){
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  nx <- ncol(X)
  if(is.null(group)) group = rep(1,nx)
  gunique = unique(group)
  G = length(unique(group))
  p <- nx/G
  if(is.null(r1)) r1 <- 2 
  if(is.null(r2)) r2 <- 2
  if(is.null(r3)) r3 <- 2
  if(is.null(r4)) r4 <- 2
  if(degr>K-1) stop("K must be larger than degree+1 !")

  opts = list(eps=eps,eps1=eps,max_step=max_step,max_step1=max_step,
              n=n,r1=r1,r2=r2,r3=r3,r4=r4,p=p,q=q,G=G,K=K,degr=degr) 
  # initial A,B,C,D,S
  if(is.null(D0)){
    set.seed(1)
    A = rbind(diag(r1), matrix(0,p-r1,r1))
    B = rbind(diag(r2), matrix(0,K-r2,r2))
    C = rbind(diag(r3), matrix(0,G-r3,r3))
    D = rbind(diag(r4), matrix(0,q-r4,r4))
    S = matrix(rnorm(r1*r2*r3*r4),r4,r1*r2*r3)
    D0 = list(S=S,A=A,B=B,C=C,D=D)
  }
  else{
    A = D0$A
    B = D0$B
    C = D0$C
    D = D0$D
    S = D0$S
  }
  
  gunique <- unique(group)
  Z = NULL; for(i in 1:G)  Z = cbind(Z, bsbasefun(X[,group==gunique[i]],K,degr))
  Zbar = colMeans(Z)
  Ybar = colMeans(Y)
  Z = Z - matrix(rep(Zbar,each=n),n)
  Y1 = Y - matrix(rep(Ybar,each=n),n)
  fit = EstimationT4(Y1,Z,S,A,B,C,D,opts)
  Dn = fit$Dnew %*% fit$Snew %*%t(kronecker(fit$Cnew, kronecker(fit$Bnew,fit$Anew)))
  if(intercept)  mu = Ybar-Dn%*%Zbar
  else mu = rep(0,q)
  return(list(rss=fit$likhd,
              Anew = fit$Anew,
              Bnew = fit$Bnew,
              Cnew = fit$Cnew,
              Dnew = fit$Dnew,
              Snew = fit$Snew,
              mu = mu,
              degr = degr,
              K = K
              )
         )
}