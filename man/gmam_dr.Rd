\name{gmam_dr}
\alias{gmam_dr-function}
\alias{gmam_dr}
\docType{package}
\title{
  Fit GMAM without sparsity assumption and ranks selected by \code{BIC} or \code{CV}.
}
\description{
  Fit a low-dimensional grouped multivariate additive model using B-splines, without aparsity assumptions, and with ranks \eqn{r_1, r_2, r_3} selected by \code{BIC} or \code{CV}.
} 

\usage{gmam_dr(Y, X, group, method="BIC", ncv=10, K_index=NULL, r1_index=NULL, r2_index=NULL, r3_index=NULL, D0=NULL,
              intercept=TRUE, degr=3, eps=1e-4, max_step=20, eps1=1e-4, max_step1=20)
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
                  tensor. Default is 
                  \code{r3_index}\eqn{=1,\cdots,\min(\log(n)\rceil,q)}.
                  }
                   
  \item{D0}{A user-specified list of initial coefficient tensor of all groups. Each element in \eqn{D0} is still a list,                 which contains initial values of coefficient matrix \eqn{S}, \eqn{A}, \eqn{B}, \eqn{C}. By default, initial                  matrices are provided by random.}
  
  \item{Intercept}{A logical value indicating whether to estimate the intercept \eqn{\mu}. Default is \code{TRUE}.}
               
  \item{degr}{ the number of knots of B-spline base function. Default is \code{3}.}
  
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
             relative change in any coefficient is less than \code{eps}.  Default
             is \code{1e-4}.}
             
  \item{max_step}{Maximum number of iterations.  Default is 20.}
  
  \item{eps1}{Convergence threshhold in updating \eqn{A} with sparsity assumption. The algorithm of updating \eqn{A} iterates               until the relative change in the row of \eqn{A} is less than \code{eps1}. Default is \code{1e-4}.}
  
  \item{max_step1}{Maximum number of iterations in updating \eqn{A} with sparsity assumption.  Default is 20.}
  
}

\details{This function gives \code{pxq} functional coefficients' estimators of GMAM. Each singular value matrix of 
  tensor is a \eqn{r_1\times r_2\times r_3}-tensor. We choose \eqn{r_1}, \eqn{r_2} 
  and \eqn{r_3}  by \code{BIC} or \code{CV}.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Snew}{An estimator list of \eqn{S_{(3)}} for all groups, where each \eqn{S} is the core tensor of \eqn{D}.}
  
  \item{Anew}{An estimator list of \eqn{A} for all groups, where each \eqn{A} is 1-mode slice of \eqn{D}.}
  
  \item{Bnew}{An estimator list of \eqn{B} for all groups, where each \eqn{B} is 2-mode slice of \eqn{D}.}
  
  \item{Cnew}{An estimator list of \eqn{C} for all groups, where each \eqn{C} is 3-mode slice of \eqn{D}.}
  
  \item{mu}{Estimator of intercept.}
  
  \item{rss}{Residual sum of squares (RSS).}
  
  \item{rk_opt}{The optimal parametres that slected by \code{BIC} or \code{CV}. It is a vector with length 4, which are selected \eqn{r_1},  \eqn{r_2}, \eqn{r_3}, and \eqn{K}.}
  
  \item{selected}{Which \eqn{\lambda} is selection.}
  
  \item{opts}{Optimal ranks and \eqn{K}.}
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
  r10 = r20 = r30 = 2
  S3 <- matrix(runif(r10*r20*r30,3,7), nrow = r30)
  T1 <- matrix(rnorm(s*r10), nrow = s)
  U1 <- qr.Q(qr(T1))
  T1 <- matrix(rnorm(K*r20), nrow = K)
  U2 <- qr.Q(qr(T1))
  T1 <- matrix(rnorm(q*r30), nrow = q)
  U3 <- qr.Q(qr(T1))
  D3 <- U3 \%*\% S3 \%*\% t(kronecker(U2,U1))
  D2 <- TransferModalUnfoldingsT(D3,3,2,c(s,K,q))
  mydata <- generateData(n,q,p,s,D2)
  Y <- mydata$Y  
  X <- mydata$X
  group <- rep(1:G,each=p/G)
  
  fit <- gmam_dr(Y, X, group)
  Dhat = NULL 
  for (g in 1:G) Dhat <- cbind(Dhat, fit$Cnew[[g]] \%*\% fit$Snew[[g]] \%*\% t(kronecker(fit$Bnew[[g]], fit$Anew[[g]]))) 
  mu <- fit$mu
  opt <- fit$rk_opt
}