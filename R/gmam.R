
##--------------without sparsity----------------------##
gmam <- function(Y,X,group=NULL,K=6,r1=NULL,r2=NULL,r3=NULL,D0=NULL,intercept=TRUE,degr=3,eps=1e-4,max_step=20){
  n <- nrow(Y)
  q <- ncol(Y)
  nx <- ncol(X)
  if(is.null(group)) group = rep(1,nx)
  gunique = unique(group)
  G = length(gunique)
  p = rep(0,G)
  for(g in 1:G) p[g] = sum(group==gunique[g])
  if(is.null(r1)) r1 <- 2 
  if(is.null(r2)) r2 <- 2
  if(is.null(r3)) r3 <- 2
  if(degr>K-1) stop("K must be larger than degree+1 !")
  # initial A,B,C,S
  if(is.null(D0)){
    set.seed(1)
    B = rbind(diag(r2), matrix(0,K-r2,r2))
    C = rbind(diag(r3), matrix(0,q-r3,r3))
    S = matrix(rnorm(r1*r2*r3),r3,r1*r2)
    D0 = list();    
    for(j in 1:G){
      A = rbind(diag(r1), matrix(0,p[j]-r1,r1))
      SABC = list(S=S,A=A,B=B,C=C)
      D0[[j]] = SABC
    }
  }
  opts = list(eps=eps,eps1=eps,max_step=max_step,max_step1=max_step,n=n,r1=r1,r2=r2,r3=r3,p=p,q=q,degr=degr,K=K,G=G,nx=nx)
  Sinit = list()
  Ainit = list()
  Binit = list()
  Cinit = list()
  Z = list()
  Zbar = NULL
  for(i in 1:G){
    Sinit[[i]] = as.matrix(D0[[i]]$S[1:r3,1:(r1*r2)])
    Ainit[[i]] = as.matrix(D0[[i]]$A[,1:r1])
    Binit[[i]] = as.matrix(D0[[i]]$B[,1:r2])
    Cinit[[i]] = as.matrix(D0[[i]]$C[,1:r3])
    Z[[i]] = as.matrix(bsbasefun(X[,group==gunique[i]],K,degr))
    Zbar1 = colMeans(Z[[i]])
    Z[[i]] = Z[[i]] - matrix(rep(Zbar1,each=n),n)
    Zbar = cbind(Zbar, Zbar1)
  }
  Ybar = colMeans(Y)
  Y1 = Y - matrix(rep(Ybar,each=n),n)
  fit = Estimation(Y1,Z,Sinit,Ainit,Binit,Cinit,opts)
  Dn = NULL
  for (g in 1:G) Dn = cbind(Dn, fit$Cnew[[g]] %*% fit$Snew[[g]] %*%t(kronecker(fit$Bnew[[g]], fit$Anew[[g]])))
  if(intercept)  mu = Ybar-Dn%*%as.vector(Zbar)
  else mu = rep(0,q)
  return(list(Snew=fit$Snew,
              Anew=fit$Anew, 
              Bnew=fit$Bnew, 
              Cnew=fit$Cnew, 
              rss=fit$likhd,
              mu = mu,
              opts = opts
              )
         )
}