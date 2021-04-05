\name{gmam4}
\alias{gmam4-function}
\alias{gmam4}
\docType{package}
\title{
  Fit GMAM without sparsity assumption and with fixed ranks.
}
\description{
  Fit a low-dimensional grouped multivariate additive model using B-splines, without aparsity assumptions, and given ranks given ranks \eqn{r_1, r_2, r_3, r_4}.
}

\usage{gmam4(Y, X, group=NULL, K=6, r1=NULL, r2=NULL, r3=NULL, r4=NULL, D0=NULL, 
             intercept=TRUE, degr=3, eps=1e-4, max_step=20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A \eqn{n\times q} numeric matrix of responses.}
  
  \item{X}{A \eqn{n\times p} numeric design matrix for the model.}
  
  \item{group}{The group index of covariates. It is a p-vector and varies from a integer of 1 to group size G}
  
  \item{K}{The number of B-spline base function, that is the plus of both degrees of base functioin and the number of knots. Default is \code{6}.}
  
  \item{r1}{The first dimension of single value matrix of the tensor. Default is 2.}
  
  \item{r2}{The second dimension of single value matrix of the tensor. Default is 2.}
  
  \item{r3}{The third dimension of single value matrix of the tensor. Default is 2.}
  
  \item{r4}{The fourth dimension of single value matrix of the tensor. Default is 2.}
  
  \item{D0}{A user-specified list of initial coefficient matrices of \eqn{S}, \eqn{A}, \eqn{B}, \eqn{C}, \eqn{D}. By default,             initial matrices are provided by random.}
  
  \item{Intercept}{A logical value indicating whether to estimate the intercept \eqn{\mu}. Default is \code{TRUE}.}
  
  \item{degr}{The number of knots of B-spline base function. Default is \code{3}.}

  \item{eps}{Convergence threshhold.  The algorithm iterates until the
    relative change in any coefficient is less than \code{eps}.  Default
    is \code{1e-4}.}
    
  \item{max_step}{Maximum number of iterations.  Default is 20.}
}

\details{This function gives \code{pxq} functional coefficients' estimators of MAM. The singular value matrix of 
  tensor is a \eqn{r_1\times r_2\times r_3\times r_4}-tensor and \eqn{r_1}, \eqn{r_2}, \eqn{r_3} 
  and \eqn{r_4}  are fixed.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Anew}{Estimator of \eqn{A}, which is 1-mode slice of \eqn{D}.}
  
  \item{Bnew}{Estimator of \eqn{B}, which is 2-mode slice of \eqn{D}.}
  
  \item{Cnew}{Estimator of \eqn{C}, which is 3-mode slice of \eqn{D}.}
  
  \item{Dnew}{Estimator of \eqn{D}, which is 4-mode slice of \eqn{D}.}
  
  \item{mu}{Estimator of intercept.}
  
  \item{rss }{Residual sum of squares (RSS).}
  
  \item{degr}{The number of knots of B-spline base function.}
  
  \item{K}{The number of B-spline base function, that is the plus of both degrees of base functioin and the number of knots.}
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
  gmam4_sparse
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
  
  fit <- gmam4(Y, X, group, K, r10, r20, r30, r40)
  Dhat = fit$Dnew \%*\% fit$Snew \%*\% t(kronecker(fit$Cnew, kronecker(fit$Bnew, fit$Anew)))
  mu <- fit$mu
}