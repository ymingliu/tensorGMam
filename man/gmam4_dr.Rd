\name{gmam4_dr}
\alias{gmam4_dr-function}
\alias{gmam4_dr}
\docType{package}
\title{
  Fit GMAM without sparsity assumption, and with ranks selected by \code{BIC} or \code{CV}.
}
\description{
  Fit a low-dimensional grouped multivariate additive model using B-splines, without aparsity assumptions, and with ranks \eqn{r_1, r_2, r_3, r_4} selected by \code{BIC} or \code{CV}.
} 

\usage{gmam4_dr(Y, X, group=NULL, method="BIC", ncv=10, K_index=NULL, r1_index=NULL, r2_index=NULL, r3_index=NULL,
                r4_index=NULL, D0=NULL, intercept=TRUE, degr=3, eps=1e-4, max_step=20, eps1=1e-4, max_step1=20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A \eqn{n\times q} numeric matrix of responses.}
  
  \item{X}{A \eqn{n\times p} numeric design matrix for the model.}
  
  \item{group}{The group index of covariates. It is a p-vector and varies from a integer of 1 to group size G}
  
  \item{method}{The method to be applied to select parameters.  Either "BIC"
                (the default), or "CV".}
    
  \item{ncv}{The number of cross-validation folds.  Default is 10. If \code{method} is  
             \code{"BIC"}, \code{ncv} is useless.}
  
  \item{K_index}{A user-specified sequence of \code{K} values, where \code{K} is he number of B-spline base function. Default is \code{k_index=6}.}
    
  \item{r1_index}{A user-specified sequence of \eqn{r_1} values, where 
                   \eqn{r_1} is the first dimension of single value matrix of the 
                   tensor. Default is 
                   \code{r1_index}\eqn{=1,\cdots,\min(\log(n)\rceil,p)}.
                   }
   
  \item{r2_index}{A user-specified sequence of \eqn{r_2} values, where 
                  \eqn{r_2} is the second dimension of single value matrix of the 
                  tensor. Default is \code{r2_index}\eqn{=1,\cdots},\code{max{K_index}}.
                  }
                  
  \item{r3_index}{A user-specified sequence of \eqn{r_3} values, where 
                  \eqn{r_3} is the third dimension of single value matrix of the 
                  tensor. Default is \code{r3_index}\eqn{=1,\cdots},\code{max{ng}}.
                  }
   
  \item{r4_index}{A user-specified sequence of \eqn{r_4} values, where 
                  \eqn{r_4} is the fourth dimension of single value matrix of the  
                  tensor. Default is 
                  \code{r4_index}\eqn{=1,\cdots,\min(\log(n)\rceil,q)}.
                  }
                   
  \item{D0}{A user-specified list of initial coefficient matrices of \eqn{S}, \eqn{A}, \eqn{B}, \eqn{C}, \eqn{D}. By default,             initial matrices are provided by random.}
  
  \item{Intercept}{A logical value indicating whether to estimate the intercept \eqn{\mu}. Default is \code{TRUE}.}
               
  \item{degr}{ the number of knots of B-spline base function. Default is \code{3}.}
  
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
             relative change in any coefficient is less than \code{eps}.  Default
             is \code{1e-4}.}
             
  \item{max_step}{Maximum number of iterations.  Default is 20.}
  
  \item{eps1}{Convergence threshhold in updating \eqn{A} with sparsity assumption. The algorithm of updating \eqn{A} iterates               until the relative change in the row of \eqn{A} is less than \code{eps1}. Default is \code{1e-4}.}
  
  \item{max_step1}{Maximum number of iterations in updating \eqn{A} with sparsity assumption.  Default is 20.}

}

\details{This function gives \code{pxq} functional coefficients' estimators of GMAM. The singular value matrix of 
  tensor is a \eqn{r_1\times r_2\times r_3\times r_4 }-tensor. We choose \eqn{r_1}, \eqn{r_2}, \eqn{r_3}
  and \eqn{r_4}  by \code{BIC} or \code{CV}.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Snew}{Estimator of \eqn{S_{(3)}}, where \eqn{S} is the core tensor of \eqn{D}.}
  
  \item{Anew}{Estimator of \eqn{A}, which is 1-mode slice of \eqn{D}.}
  
  \item{Bnew}{Estimator of \eqn{B}, which is 2-mode slice of \eqn{D}.}
  
  \item{Cnew}{Estimator of \eqn{C}, which is 3-mode slice of \eqn{D}.}
  
  \item{Dnew}{Estimator of \eqn{D}, which is 4-mode slice of \eqn{D}.}
  
  \item{mu}{Estimator of intercept.}
  
  \item{rss}{Residual sum of squares (RSS).}
  
  \item{rk_opt}{The optimal parametres that slected by \code{BIC} or \code{CV}. It is a vector with length 5, which are selected \eqn{r_1},  \eqn{r_2}, \eqn{r_3} and \eqn{r_4}, and \eqn{K}.}
  
  \item{selected}{Which \eqn{\lambda} is selection.}
  
  \item{opts}{The parameter settings, optimal ranks and \eqn{K}}
}

%\author{
%Your Name, email optional.
%Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
%}
\references{A tensor estimation approach to integrative mulit-view multivariate additive models.
}
\seealso{
  gmam, gmam_sparse_dr
}
\examples{ 
  n <- 500
  p <- 20
  q <- 10
  G <- 4
  s <- p/G
  K <- 6
  r10 = r20 = r30 = r40 = 2
  S4 <- matrix(runif(r10*r20*r30*r40,3,7), nrow = r40)
  T1 <- matrix(rnorm(s*r10), nrow = s)
  U1 <- qr.Q(qr(T1))
  T1 <- matrix(rnorm(K*r20), nrow = K)
  U2 <- qr.Q(qr(T1))
  T1 <- matrix(rnorm(G*r30), nrow = G)
  U3 <- qr.Q(qr(T1))
  T1 <- matrix(rnorm(q*r40), nrow = q)
  U4 <- qr.Q(qr(T1))
  D4 <- U4 \%*\% S4 \%*\% t(kronecker(U3, kronecker(U2,U1)))
  D42 <- TransferModalUnfoldingsT(D4,4,2,c(s,K,G,q))
  mydata <- generateDataT4(n,q,p,s,D42)
  Y <- mydata$Y  
  X <- mydata$X
  group <- rep(1:G,each=p/G)
  
  fit <- gmam4_dr(Y, X, group)
  Dhat = fit$Dnew \%*\% fit$Snew \%*\% t(kronecker(fit$Cnew, kronecker(fit$Bnew, fit$Anew)))
  mu <- fit$mu
  opt <- fit$rk_opt
}