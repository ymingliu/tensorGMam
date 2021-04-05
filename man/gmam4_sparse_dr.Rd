\name{gmam4_sparse_dr}
\alias{gmam4_sparse_dr-function}
\alias{gmam4_sparse_dr}
\docType{package}
\title{
  Fit GMAM with sparsity assumption and ranks selected by \code{BIC} or \code{CV}.
}
\description{
  Fit a high-dimensional grouped multivariate additive model using B-splines, with aparsity assumptions and ranks selected by \code{BIC} or \code{CV}. The multivariate sparse group \code{lasso} (\code{mcp} or \code{scad}) and the coordinate descent algorithm are used to estimate
  functions for sparsity situation. The tuning parameter is selected by \code{BIC} or \code{CV}, which matchs the method of rank selection.
}

\usage{
gmam4_sparse_dr(Y, X, group, K_index=NULL, r1_index=NULL, r2_index=NULL, r3_index=NULL, r4_index=NULL, 
                method="BIC", ncv=10, penalty="LASSO", lambda=NULL, D0=NULL, intercept=TRUE, nlam=20,
                degr=3, lam_min=0.01, eps=1e-4, max_step=10, eps1=1e-4, max_step1=10, gamma=2, dfmax=NULL, alpha=1)
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
                  \code{r3_index}\eqn{=1,\cdots,\min(\lceil\log(n)\rceil,ng)}.
                  }
                  
  \item{r4_index}{A user-specified sequence of \eqn{r_4} values, where 
                  \eqn{r_4} is the fourth dimension of single value matrix of the  
                  tensor. Default is 
                  \code{r4_index}\eqn{=1,\cdots,\min(\lceil\log(n)\rceil,q)}.
                  }                

  \item{method}{The method to be applied to select parameters.  Either "BIC"
                (the default), or "CV".}
    
  \item{ncv}{The number of cross-validation folds.  Default is 10. If \code{method} is \code{BIC}, \code{ncv} is useless.}
  
  \item{penalty}{The penalty to be applied to the model. Either "LASSO" (the default), 
                 "SCAD", or "MCP".}
                  
  \item{lambda}{A user-specified sequence of lambda values.  By default,
                a sequence of values of length \code{nlam} is computed, equally
                spaced on the log scale.}

  \item{D0}{A user-specified list of initial coefficient matrices of \eqn{S}, \eqn{A}, \eqn{B}, \eqn{C}, \eqn{D}. By default,             initial matrices are provided by random.}
  
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
  This function gives \code{pxq} functional coefficients' estimators of GMAM. The singular value matrix of 
  tensor is a \eqn{r_1\times r_2\times r_3\times r_4}-tensor. We choose \eqn{r_1}, \eqn{r_2}, \eqn{r_3}
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
  
  \item{df}{Degrees of freedom.}
  
  \item{activeX}{The active set.}
  
  \item{lambda}{The sequence of regularization parameter values in the path.}
  
  \item{selectedID }{The index of \code{lambda} corresponding to \code{lambda_opt}.}
  
  \item{lambda_opt }{The value of \code{lambda} with the minimum \code{BIC} or \code{CV} value.}
  
  \item{RSS}{The values of \code{BIC} or \code{CV}, which is a vector.}
  
  \item{rk_opt}{The optimal parametres that slected by \code{BIC} or \code{CV}. It is a vector with length 4, which are selected \eqn{r_1},  \eqn{r_2}, \eqn{r_3}, and \eqn{K}.}
  
  \item{opts}{The parameter settings, optimal ranks and \eqn{K}}
  
  \item{opts_pen}{The parameter settings related to penalty item}
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
  gmam4_dr, gmam4_sparse
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
  
  fit <- gmam4_sparse_dr(Y, X, group)
  Dhat = fit$Dnew \%*\% fit$Snew \%*\% t(kronecker(fit$Cnew, kronecker(fit$Bnew, fit$Anew)))
  mu <- fit$mu
  opt <- fit$rk_opt

}