\name{gmam_sparse_dr}
\alias{gmam_sparse_dr-function}
\alias{gmam_sparse_dr}
\docType{package}
\title{
  Fit GMAM with sparsity assumption and ranks selected by \code{BIC} or \code{CV}.
}
\description{
  Fit a high-dimensional grouped multivariate additive model using B-splines, with aparsity assumptions and ranks selected by \code{BIC} or \code{CV}. The multivariate sparse group \code{lasso} (\code{mcp} or \code{scad}) and the coordinate descent algorithm are used to estimate
  functions for sparsity situation. The tuning parameter is selected by \code{BIC} or \code{CV}, which matchs the method of rank selection.
}

\usage{
gmam_sparse_dr(Y, X, group=NULL, K_index=NULL, r1_index=NULL, r2_index=NULL, r3_index=NULL,
               method="BIC", ncv=10, penalty="LASSO", lambda=NULL, D0=NULL, intercept=TRUE,
               nlam=50, degr=3, lam_min=0.01, eps=1e-4, max_step=20, eps1=1e-4, max_step1=20, 
               gamma=2, dfmax=NULL, alpha=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A \eqn{n\times q} numeric matrix of responses.}
  
  \item{X}{A \eqn{n\times p} numeric design matrix for the model.}
  
  \item{group}{The group index of covariates. It is a p-vector and varies from a integer of 1 to group size G}
  
  \item{K_index}{A user-specified sequence of \code{K} values, where \code{K} is he number 
                 of B-spline base function. Default is \code{k_index=6}.}
                 
  \item{r1_index}{A user-specified sequence of \eqn{r_1} values, where 
                  \eqn{r_1} is the first dimension of single value matrix of the  
                  tensor. Default is 
                  \code{r1_index}\eqn{=1,\cdots,\min(\lceil\log(n)\rceil,p)}.
                  }
            
  \item{r2_index}{A user-specified sequence of \eqn{r_2} values, where 
                  \eqn{r_2} is the second dimension of single value matrix of the
                  tensor. Default is \code{r2_index}\eqn{=1,\cdots},\code{max{K_index}}.
                  }
                  
  \item{r3_index}{A user-specified sequence of \eqn{r_3} values, where 
                  \eqn{r_3} is the third dimension of single value matrix of the  
                  tensor. Default is 
                  \code{r3_index}\eqn{=1,\cdots,\min(\lceil\log(n)\rceil,q)}.
                  }
  
  \item{method}{The method to be applied to select parameters.  Either "BIC"
                (the default), or "CV".}
    
  \item{ncv}{The number of cross-validation folds.  Default is 10. If \code{method} is \code{BIC}, \code{ncv} is useless.}
  
  \item{penalty}{The penalty to be applied to the model. Either "LASSO" (the default), 
                 "SCAD", or "MCP".}
                 
  \item{lambda}{A user-specified sequence of lambda values.  By default,
                a sequence of values of length \code{nlam} is computed, equally
                spaced on the log scale.}
                
  \item{D0}{A user-specified list of initial coefficient tensor of all groups. Each element in \eqn{D0} is still a list,                 which contains initial values of coefficient matrix \eqn{S}, \eqn{A}, \eqn{B}, \eqn{C}. By default, initial                  matrices are provided by random.}
  
  \item{Intercept}{A logical value indicating whether to estimate the intercept \eqn{\mu}. Default is \code{TRUE}.}
  
  \item{nlam}{The number of lambda values.  Default is 50.}
  
  \item{degr}{The number of knots of B-spline base function. Default is \code{degr = 3}.}
  
  \item{lam_min}{The smallest value for lambda, as a fraction of
                 lambda.max.  Default is 1e-2.}
                 
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
             relative change in any coefficient is less than \code{eps}.  Default
             is \code{1e-4}.}
             
  \item{max_step}{Maximum number of iterations.  Default is 20.}
  
  \item{eps1}{Convergence threshhold in updating \eqn{A} with sparsity assumption. The algorithm of updating \eqn{A} iterates               until the relative change in the row of \eqn{A} is less than \code{eps1}. Default is \code{1e-4}.}
  
  \item{max_step1}{Maximum number of iterations in updating \eqn{A} with sparsity assumption.  Default is 20.}
  
  \item{gamma}{The tuning parameter of the MCP/SCAD penalty (see details).}
  
  \item{dfmax}{Upper bound for the number of nonzero coefficients.
               Default is no upper bound.  However, for large data sets,
               computational burden may be heavy for models with a large number of
               nonzero coefficients.}
    
  \item{alpha}{Tuning parameter for the Mnet estimator which controls
               the relative contributions from the LASSO, MCP/SCAD penalty and the ridge,
               or L2 penalty.  \code{alpha=1} is equivalent to LASSO, MCP/SCAD penalty,
               while \code{alpha=0} would be equivalent to ridge regression.
               However, \code{alpha=0} is not supported; \code{alpha} may be
               arbitrarily small, but not exactly 0.}

}

\details{
  This function gives \code{pxq} functional coefficients' estimators of GMAM. Each singular value matrix of 
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
  
  \item{df}{Degrees of freedom.}
  
  \item{activeX}{The active set.}
  
  \item{lambda}{The sequence of regularization parameter values in the path.}
  
  \item{selectedID }{The index of \code{lambda} corresponding to \code{lambda_opt}.}
  
  \item{lambda_opt }{The value of \code{lambda} with the minimum \code{BIC} or \code{CV} value.}
  
  \item{RSS}{The values of \code{BIC} or \code{CV}, which is a vector.}
  
  \item{rk_opt}{The optimal parametres that slected by \code{BIC} or \code{CV}. It is a vector with length 4, which are selected \eqn{r_1},  \eqn{r_2}, \eqn{r_3}, and \eqn{K}.}
  
  \item{opts}{Optimal ranks and \eqn{K}.}
  
  \item{opts_pen}{The parameter settings related to penalty item.}
  
  %\item{...}{ Other options for CompositeQuantile.}
}

%\author{
%Your Name, email optional.
%Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
%}
\references{
  A tensor estimation approach to integrative mulit-view multivariate additive models.
}
\seealso{
  gmam_dr, gmam_sparse
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
  
  fit <- gmam_sparse_dr(Y, X, group)
  Dhat = NULL 
  for (g in 1:G) Dhat <- cbind(Dhat, fit$Cnew[[g]] \%*\% fit$Snew[[g]] \%*\% t(kronecker(fit$Bnew[[g]], fit$Anew[[g]]))) 
  mu <- fit$mu
  opt <- fit$rk_opt
}