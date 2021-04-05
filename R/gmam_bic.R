
##--------------main by BIC without sparsity----------------------##
gmam_bic <- function(Y,X,group,K_index,r1_index,r2_index,r3_index,D0,intercept,opts){
  n = opts$n
  p = opts$p
  q = opts$q
  G = opts$G
  degr = opts$degr
  gunique <- unique(group)
  Ybar = colMeans(Y)
  Y1 = Y - matrix(rep(Ybar,each=n),n)
  Sinit = list()
  Ainit = list()
  Binit = list()
  Cinit = list()
  Z = list()
  RSS = NULL
  for(K in K_index){
    for(i in 1:G){
      Z[[i]] = as.matrix(bsbasefun(X[,group==gunique[i]],K,degr))
      Zbar1 = colMeans(Z[[i]])
      Z[[i]] = Z[[i]] - matrix(rep(Zbar1,each=n),n)
    }
    opts$K = K
    for(r3 in r3_index){
      opts$r3 = r3
      for(r2 in r2_index){
        opts$r2 = r2
        for(r1 in r1_index){
          opts$r1 = r1
          for(g in 1:G){
            Sinit[[g]] = as.matrix(D0[[g]]$S[1:r3,1:(r1*r2)])
            Ainit[[g]] = as.matrix(D0[[g]]$A[,1:r1])
            Binit[[g]] = as.matrix(D0[[g]]$B[,1:r2])
            Cinit[[g]] = as.matrix(D0[[g]]$C[,1:r3])
          }
          fit = Estimation(Y1,Z,Sinit,Ainit,Binit,Cinit,opts)
          df = 0; for(g in 1:G) df = df + r1*r2*r3 + p[g]*r1 + K*r2 + q*r3 - r1^2 - r2^2 - r3^2
          RSS = c(RSS,fit$likhd+log(n)*df)
          #stop("stop here!")
        }
      }
    }
  }
  selected = which.min(RSS)
  opt = assig(c(length(r1_index),length(r2_index),length(r3_index),length(K_index)))[,selected]
  r1_opt = r1_index[opt[1]]
  r2_opt = r2_index[opt[2]]
  r3_opt = r3_index[opt[3]]
  K_opt = K_index[opt[4]]
  #---------------- The estimation after selection ---------------------#
  opts$r1 = r1_opt
  opts$r2 = r2_opt
  opts$r3 = r3_opt
  opts$K = K_opt
  Zbar = NULL
  for(i in 1:G){
    Sinit[[i]] = as.matrix(D0[[i]]$S[1:r3_opt,1:(r1_opt*r2_opt)])
    Ainit[[i]] = as.matrix(D0[[i]]$A[,1:r1_opt])
    Binit[[i]] = as.matrix(D0[[i]]$B[,1:r2_opt])
    Cinit[[i]] = as.matrix(D0[[i]]$C[,1:r3_opt])
    Z[[i]] = as.matrix(bsbasefun(X[,group==gunique[i]],K_opt,degr))
    Zbar1 = colMeans(Z[[i]])
    Z[[i]] = Z[[i]] - matrix(rep(Zbar1,each=n),n)
    Zbar = cbind(Zbar, Zbar1)
  }
  fit = Estimation(Y1,Z,Sinit,Ainit,Binit,Cinit,opts) 
  Dn = NULL
  for (g in 1:G) Dn = cbind(Dn, fit$Cnew[[g]] %*% fit$Snew[[g]] %*%t(kronecker(fit$Bnew[[g]], fit$Anew[[g]])))
  if(intercept)  mu = Ybar-Dn%*%as.vector(Zbar)
  else mu = rep(0,q)
  return(list(Snew=fit$Snew,
              Anew=fit$Anew, 
              Bnew=fit$Bnew, 
              Cnew=fit$Cnew, 
              mu = mu,
              rss=fit$likhd,
              rk_opt=c(r1_opt,r2_opt,r3_opt,K_opt),
              selected=selected,
              opts = opts
              )
         )
}