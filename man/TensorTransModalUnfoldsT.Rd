\name{TransferModalUnfoldingsT}
\alias{TransferModalUnfoldingsT-function}
\alias{TransferModalUnfoldingsT}
\docType{package}
\title{Transfer a tensor's modal unfoldings to another.}
\description{
  Transfer a tensor's modal unfoldings to another.
}

\usage{
TransferModalUnfoldingsT(S, d1, d2, dim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
  A mode-d1-unfolding of a tensor, input unfolding
}
  \item{d1}{
  An integer, the mode of unfolding \eqn{S_{(d_1)}}
}
  \item{d2}{
  An integer, the mode of output unfolding \eqn{S_{(d_2)}}
}
  \item{dim}{
  The dimension vector of tensor \eqn{S} and \eqn{length(dim) \ge 3}. e.g. for a tensor with size \eqn{r_1xr_2xr_3}, \eqn{dim = c(r_1, r_2, r_3)}
}
}

\details{
  This function transfers an input mode-d1-unfolding \eqn{S_{(d_1)}} 
  to mode-d2-unfolding \eqn{S_{(d_2)}}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{D}{the output mode-d2-unfolding, \eqn{S_{(d_2)}}}
}

%\author{
%Your Name, email optional.
%Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
%}
\references{
  A tensor estimation approach to integrative mulit-view multivariate additive models.
}
\keyword{ High-dimensional, Sparse models; Tensor estimation; Tucker decomposition. }
\examples{
  D1 <- matrix(1:48,nrow = 4) # A tensor unfolding with size 4*12
  D2 <- TransferModalUnfoldingsT(D1,1,2,c(4,3,2,2))
}
