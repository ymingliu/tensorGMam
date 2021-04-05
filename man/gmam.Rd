\name{gmam}
\alias{gmam-function}
\alias{gmam}
\docType{package}
\title{
  Fit GMAM without sparsity assumption and with fixed ranks.
}
\description{
  Fit a low-dimensional grouped multivariate additive model using B-splines, without aparsity assumptions, and given ranks given ranks \eqn{r_1, r_2, r_3}.
}

\usage{gmam(Y, X, group=NULL, K=6, r1=NULL, r2=NULL, r3=NULL, D0=NULL, intercept=TRUE, degr=3, eps=1e-4, max_step=20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A \eqn{n\times q} numeric matrix of responses.}
  
  \item{X}{A \eqn{n\times p} numeric design matrix for the model.}
  
  \item{group}{The group index of covariates. It is a p-vector and varies from a integer of 1 to group size G}
  
  \item{K}{The number of B-spline base function, that is the plus of both degrees of base functioin and the number of knots. Default is \code{6}.}
  
  \item{r1}{The first dimension of single value matrix of each tensor. Default is 2.}
  
  \item{r2}{The second dimension of single value matrix of each tensor. Default is 2.}
  
  \item{r3}{The third dimension of single value matrix of each tensor. Default is 2.}
  
  \item{D0}{A user-specified list of initial coefficient tensor of all groups. Each element in \eqn{D0} is still a list,                 which contains initial values of coefficient matrix \eqn{S}, \eqn{A}, \eqn{B}, \eqn{C}. By default, initial                  matrices are provided by random.}
              
  \item{degr}{ the number of knots of B-spline base function. Default is \code{3}.}
  
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
    relative change in any coefficient is less than \code{eps}.  Default
    is \code{1e-4}.}
    
  \item{max_step}{Maximum number of iterations.  Default is 20.}
}

\details{This function gives \code{pxq} functional coefficients' estimators of GMAM. Each singular value matrix of 
  tensor is a \eqn{r_1\times r_2\times r_3}-tensor and \eqn{r_1}, \eqn{r_2} 
  and \eqn{r_3} are fixed.
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
  
  \item{opts}{The parameter settings to obtain estimator of \eqn{D}}
  %\item{...}{ Other options for CompositeQuantile.}
}

%\author{
%Your Name, email optional.
%Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
%}
\references{
A tensor estimation approach to integrative mulit-view multivariate additive models.
}
\keyword{High-dimensional; Sparse models; Tensor estimation; Tucker decomposition. }
\seealso{
  gmam_sparse
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
  
  fit <- gmam(Y, X, group, K, r10, r20, r30)
  Dhat = NULL 
  for (g in 1:G) Dhat <- cbind(Dhat, fit$Cnew[[g]] \%*\% fit$Snew[[g]] \%*\% t(kronecker(fit$Bnew[[g]], fit$Anew[[g]]))) 
  mu <- fit$mu

}