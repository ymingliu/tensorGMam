\name{generateData}
\alias{generateData}
\docType{package}
\title{
  Generate data from GMAM model.
}
\description{
  Generate data for a high-dimensional multivariate additive model with grouped covariates, under sparsity assumptions on groups and within each group. Each group has the same size.
}

\usage{
generateData(n, q, p, s, D42, sigma2=NULL, indexF=NULL, seed_id=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Sample size.}
  
  \item{q}{The number of responses, \eqn{q \ge q1}.}
  
  \item{p}{The number of covariates, \eqn{p \ge q1}.}
  
  \item{s}{The true covariates associating to response, \eqn{s \ge q1}.}
  
  \item{D42}{The mode of unfolding \eqn{D_{(2)}}, where \eqn{D} is the 4th order tensor.}
  
  \item{sigma2}{err variance. Default is 0.1.}
  
  \item{indexF}{A \eqn{qxs} matrix. The index of significant predictors corresponding to response \eqn{y_l}. 
  Default is the matrix with each row being \eqn{(1,2,\cdots,s)}.}
  
  \item{seed_id}{A positive integer, the seed for generating the random numbers.}
}

\details{This function gives \eqn{p x q} functional coefficients' estimators of GMAM. The singular value matrix of 
  tensor is a \eqn{r_1xr_2xr_3}-tensor. We choose \eqn{r_1}, \eqn{r_2} 
  and \eqn{r_3}  by \code{BIC} or \code{CV}.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Y}{Response, a \eqn{nxq}-matrix.}
  
  \item{X}{Design matrix, a \eqn{nxp}-matrix.}
  
  \item{f0}{True functions.}
}

%\author{
%Your Name, email optional.
%Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
%}
\references{
  A tensor estimation approach to integrative mulit-view multivariate additive models.
}
\keyword{ High-dimensional; Sparse models; Tensor estimation; Tucker decomposition. }
\seealso{
  generateData
}
\examples{
  # Example 1
  n <- 500
  p <- 20
  q <- 10
  G <- 4
  s <- p/G
  D42 <- matrix(runif(p*q,3,7),p,q)
  mydata <- generateDataT4(n,q,p,s,D42)
  Y <- mydata$Y  
  X <- mydata$X
    
    
  # Example 2
  n <- 500
  p <- 20
  q <- 10
  G <- 4
  s <- p/G
  K <- 6
  r10 = r20 = r30 = r40 = 2
  S3 <- matrix(runif(r10*r20*r30*r40,3,7), nrow = r40)
  T1 <- matrix(rnorm(s*r10), nrow = s)
  U1 <- qr.Q(qr(T1))
  T1 <- matrix(rnorm(K*r20), nrow = K)
  U2 <- qr.Q(qr(T1))
  T1 <- matrix(rnorm(G*r30), nrow = G)
  U3 <- qr.Q(qr(T1))
  T1 <- matrix(rnorm(q*r40), nrow = q)
  U4 <- qr.Q(qr(T1))
  D4 <- U4 \%*\% S3 \%*\% t(kronecker(U3,kronecker(U2,U1)))
  D42 <- TransferModalUnfoldingsT(D4,4,2,c(s,K,G,q))
  mydata <- generateDataT4(n,q,p,s,D42)
  Y <- mydata$Y  
  X <- mydata$X
}
