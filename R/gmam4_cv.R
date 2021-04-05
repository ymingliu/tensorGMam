
##--------------Estimation without Penalty----------------------##
gmam4_cv <- function(Y,X,group,ncv,K_index,r1_index,r2_index,r3_index,D0,intercept,opts){
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
  
  len_cv = ceiling(n/ncv)  
  RSS = rep(0,length(r1_index)*length(r2_index)*length(r3_index)*length(K_index))
  for(jj in 1:ncv){ # start CV
    cv.id = ((jj-1)*len_cv+1):(jj*len_cv)
    if(jj==ncv) cv.id = ((jj-1)*len_cv+1):n
    Ytrain = Y1[-cv.id,]
    Xtrain = X[-cv.id,]
    Ytest = Y1[cv.id,]
    Xtest = X[cv.id,]
    
    RSS0 = NULL
    for(K in K_index){
      opts$K = K
      # Ztrain = list()
      # Ztest = list()
      # for(i in 1:G){
      #   Ztrain[[i]] = as.matrix(bsbasefun(Xtrain[,group==gunique[i]],K,degr))
      #   Ztest[[i]] = as.matrix(bsbasefun(Xtest[,group==gunique[i]],K,degr))
      # } 
      Ztrain = NULL
      Ztest = NULL
      for(i in 1:G) {
        Ztrain = cbind(Ztrain, bsbasefun(X[,group==gunique[i]],K,degr))
        Ztest = cbind(Ztest, bsbasefun(X[,group==gunique[i]],K,degr))
      } 
      for(r3 in r3_index){
        opts$r3 = r3
        for(r2 in r2_index){
          opts$r2 = r2
          for(r1 in r1_index){
            opts$r1 = r1
            S = as.matrix(D0$S[1:r4,1:(r1*r2*r3)])
            A = as.matrix(D0$A[,1:r1])
            B = as.matrix(D0$B[,1:r2])
            C = as.matrix(D0$C[,1:r3])
            D = as.matrix(D0$D[,1:r4])
            # fit = EstimationT4(Y,Z,S,A,B,C,D,opts)  
            # Dn  = fit$Dn
            # temp = Ytest
            # for(i in 1:G){
            #    Dn = fit$Dnew[[i]] %*% fit$Snew[[i]] %*%t(kronecker(fit$Cnew[[i]], kronecker(fit$Bnew[[i]],fit$Anew[[i]])))
            #    temp = temp - Ztest[[i]] %*% t(Dn)
            #  }
            fit = EstimationT4(Ytrain,Ztrain,S,A,B,C,D,opts) 
            Dn = fit$Dnew %*% fit$Snew %*%t(kronecker(fit$Cnew, kronecker(fit$Bnew,fit$Anew)))
            temp = Ytest - Ztest %*% t(Dn)
            RSS0 = c(RSS0,sum(temp^2))
          }
        }
      }
    }
    RSS = RSS + RSS0
  } # end of CV
  selected = which.min(RSS)
  opt = assig(c(length(r1_index),length(r2_index),length(r3_index),length(r4_index),length(K_index)))[,selected]
  r1_opt = r1_index[opt[1]]
  r2_opt = r2_index[opt[2]]
  r3_opt = r3_index[opt[3]]
  r4_opt = r4_index[opt[4]]
  K_opt = K_index[opt[5]]
  #---------------- The estimation after selection ---------------------#
  opts$r1 = r1_opt
  opts$r2 = r2_opt
  opts$r3 = r3_opt
  opts$r4 = r4_opt
  opts$K = K_opt
  
  A = as.matrix(D0$A[,1:r1_opt])
  B = as.matrix(D0$B[,1:r2_opt])
  C = as.matrix(D0$C[,1:r3_opt])
  D = as.matrix(D0$D[,1:r4_opt])
  S = as.matrix(D0$S[1:r4_opt,1:(r1_opt*r2_opt*r3_opt)])
  Z = NULL;  for(i in 1:G)  Z = cbind(Z, bsbasefun(X[,group==gunique[i]],K_opt,degr))
  Zbar = colMeans(Z)
  Z = Z - matrix(rep(Zbar,each=n),n)
  fit = EstimationT4(Y1,Z,S,A,B,C,D,opts)
  Dn = fit$Dnew %*% fit$Snew %*%t(kronecker(fit$Cnew, kronecker(fit$Bnew,fit$Anew)))
  if(intercept)  mu = Ybar-Dn%*%Zbar
  else mu = rep(0,q)
  return(list(Snew=fit$Snew,
              Anew=fit$Anew, 
              Bnew=fit$Bnew, 
              Cnew=fit$Cnew, 
              Dnew=fit$Dnew,
              rss=fit$likhd,
              mu = mu,
              rk_opt=c(r1_opt,r2_opt,r3_opt,r4_opt,K_opt),
              selected=selected,
              opts = opts
              )
         )
}
