
##--------------main by BIC without sparsity----------------------##
gmam_dr <- function(Y,X,group,method="BIC",ncv=10,K_index=NULL,r1_index=NULL,r2_index=NULL,r3_index=NULL,D0=NULL,
                    intercept=TRUE,degr=3,eps=1e-4,max_step=20,eps1=1e-4,max_step1=20){
  n <- nrow(Y)
  q <- ncol(Y)
  nx <- ncol(X)
  if(is.null(group)) group = rep(1,nx)
  gunique = unique(group)
  G = length(unique(group))
  p = rep(0,G)
  for(g in 1:G) p[g] = sum(group==gunique[g])
  K1 <- 6
  if(is.null(K_index)) K_index = min(6,K1-1):max(8,K1+1)
  if(is.null(r1_index)) r1_index = 1:min(ceiling(log(n)),min(p))
  if(is.null(r2_index)) r2_index = 1:min(K_index)
  if(is.null(r3_index)) r3_index = 1:min(ceiling(log(n)),q)
  if(degr>min(6,K1-1)-1) stop("K must be larger than degree+1 !")
  #---------------- The selection by BIC  ---------------------#  
  if(is.null(D0)){
    set.seed(1)
    r1_max = max(r1_index) 
    r2_max = max(r2_index) 
    r3_max = max(r3_index) 
    K_max = max(K_index)
    B = rbind(diag(r2_max), matrix(0,K_max-r2_max,r2_max))
    C = rbind(diag(r3_max), matrix(0,q-r3_max,r3_max))
    S = matrix(rnorm(r1_max*r2_max*r3_max),r3_max,r1_max*r2_max)
    D0 = list(); 
    for(j in 1:G){
      A = rbind(diag(r1_max), matrix(0,p[j]-r1_max,r1_max))
      SABC = list(S=S,A=A,B=B,C=C)
      D0[[j]] = SABC
    }
  }
  opts = list(eps=eps,eps1=eps1,max_step=max_step,max_step1=max_step1,n=n,r1=2,r2=2,r3=2,p=p,q=q,degr=degr,K=max(K_index),G=G,nx=nx)
  if(method=="BIC") fit_dr = gmam_bic(Y,X,group,K_index,r1_index,r2_index,r3_index,D0,intercept,opts)
  if(method=="CV") fit_dr = gmam_cv(Y,X,group,ncv,K_index,r1_index,r2_index,r3_index,D0,intercept,opts)
  return(fit_dr)
}